#include "GanonClassify.hpp"
#include "hierarchical_interleaved_bloom_filter.hpp"

#include <robin_hood.h>

#include <utils/IBFConfig.hpp>
#include <utils/LCA.hpp>
#include <utils/SafeQueue.hpp>
#include <utils/StopClock.hpp>
#include <utils/adjust_seed.hpp>
#include <utils/dna4_traits.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/tuple.hpp>
#include <cereal/types/vector.hpp>

#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/alphabet/views/char_to.hpp> // para la conversion 

#include <cinttypes>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include<filesystem> //necesario para calcular el tamano del archivo

#include <mpi.h>

// librerias para el mapeo del archivo a maemoria ram 
#include <sys/mman.h> // para mmap 
#include <sys/stat.h> // para struct stat y fstat
#include <fcntl.h> // para open 
#include <unistd.h> // para close 
#include <string_view> // para optimizar la conversion sin copias
		       //#include <cereal/types/utility.hpp> // para serializar std::pair con cereal



		       //para la serializacion con cereal de std::pair
namespace cereal {
	template<class Archive, class F, class S> 
		void serialize(Archive & archive, std::pair<F, S> & p)
		{
			archive(p.first, p.second); 
		}
}

namespace GanonClassify
{

	namespace detail
	{

#ifdef LONGREADS
		typedef uint32_t TIntCount;
#else
		typedef uint16_t TIntCount;
#endif

		typedef raptor::hierarchical_interleaved_bloom_filter< seqan3::data_layout::uncompressed > THIBF;
		typedef seqan3::interleaved_bloom_filter< seqan3::data_layout::uncompressed >              TIBF;
		typedef robin_hood::unordered_map< std::string, std::tuple< size_t, double > >             TMatches;
		typedef std::vector< std::tuple< size_t, std::string > >                                   TBinMap;
		typedef robin_hood::unordered_map< std::string, std::vector< size_t > >                    TMap;
		typedef robin_hood::unordered_map< std::string, double >                                   TTargetFpr;

		struct Node
		{
			std::string parent;
			std::string rank;
			std::string name;
		};


		struct ReadBatches
		{

			ReadBatches()
			{
			}

			ReadBatches( bool _paired )
			{
				paired = _paired;
			}

			ReadBatches( bool _paired, std::vector< std::string > _ids, std::vector< std::vector< seqan3::dna4 > > _seqs )
			{
				paired = _paired;
				ids    = _ids;
				seqs   = _seqs;
			}

			ReadBatches( bool                                       _paired,
					std::vector< std::string >                 _ids,
					std::vector< std::vector< seqan3::dna4 > > _seqs,
					std::vector< std::vector< seqan3::dna4 > > _seqs2 )
			{
				paired = _paired;
				ids    = _ids;
				seqs   = _seqs;
				seqs2  = _seqs2;
			}

			bool                                       paired = false;
			std::vector< std::string >                 ids;
			std::vector< std::vector< seqan3::dna4 > > seqs;
			std::vector< std::vector< seqan3::dna4 > > seqs2{};
		};

		struct ReadMatch
		{
			ReadMatch()
			{
			}

			ReadMatch( std::string _target, size_t _kmer_count )
			{
				target     = _target;
				kmer_count = _kmer_count;
			}

			std::string target;
			size_t      kmer_count;
		};

		struct ReadOut
		{
			ReadOut()
			{
			}

			ReadOut( std::string _readID )
			{
				readID = _readID;
			}

			std::string              readID;
			std::vector< ReadMatch > matches;
		};

		struct Rep
		{
			// Report with counts of matches and reads assigned (unique or lca) for each target
			size_t matches      = 0;
			size_t lca_reads    = 0;
			size_t unique_reads = 0;
			//para que cereal pueda serializar el struct
			template <class Archive>
				void serialize( Archive& ar ){
					ar (matches, lca_reads, unique_reads);
				}
		};

		typedef robin_hood::unordered_map< std::string, Rep >  TRep;
		typedef robin_hood::unordered_map< std::string, Node > TTax;

		struct Total
		{
			size_t reads_processed  = 0;
			size_t length_processed = 0;
			size_t reads_classified = 0;
			size_t matches          = 0;
			size_t unique_matches   = 0;
		};

		struct Stats
		{
			Total total;
			// number of reads in the input files
			size_t input_reads = 0;
			// Total for each hierarchy
			std::map< std::string, Total > hierarchy_total;

			void add_totals( std::string hierarchy_label, std::vector< Total > const& totals )
			{
				// add several totals (from threads) into the stats
				for ( auto const& t : totals )
				{
					total.reads_processed += t.reads_processed;
					total.length_processed += t.length_processed;
					total.reads_classified += t.reads_classified;
					hierarchy_total[hierarchy_label].reads_processed += t.reads_processed;
					hierarchy_total[hierarchy_label].reads_classified += t.reads_classified;
					hierarchy_total[hierarchy_label].length_processed += t.reads_classified;
				}
			}

			void add_reports( std::string hierarchy_label, TRep const& report )
			{
				// add values from reports to stats
				for ( auto const& [target, rep] : report )
				{
					total.matches += rep.matches;
					total.unique_matches += rep.unique_reads;
					hierarchy_total[hierarchy_label].matches += rep.matches;
					hierarchy_total[hierarchy_label].unique_matches += rep.unique_reads;
				}
			}
		};

		struct FilterConfig
		{
			FilterConfig()
			{
			}

			FilterConfig( std::string _ibf_file, double _rel_cutoff )
			{
				ibf_file   = _ibf_file;
				rel_cutoff = _rel_cutoff;
			}

			std::string ibf_file;
			std::string tax_file = "";
			double      rel_cutoff;
			IBFConfig   ibf_config;
			TTargetFpr  target_fpr;
		};

		struct HierarchyConfig
		{
			std::vector< FilterConfig > filters;
			uint8_t                     kmer_size;
			uint32_t                    window_size;
			double                      rel_filter;
			double                      fpr_query;
			std::string                 output_file_lca;
			std::string                 output_file_all;
		};

		template < typename TFilter >
			struct Filter
			{
				TFilter      ibf;
				TMap         map;
				TTax         tax;
				size_t       bin_count = 0;
				FilterConfig filter_config;
			};


		//funcion auxiliar para comprobar si estamos en una cabecera valida o es info de calidad,
		// miramos si la linea empieza por @ y si la tercera linea empieza por +.
		// Mas info de como es el formato FASTQ en https://es.wikipedia.org/wiki/Formato_FASTQ

		bool is_valid_fastq_start( std::ifstream& file ) {
			std::streampos pos_actual = file.tellg(); //guarda la pos actual del puntero de lectura 
			std::string l1, l2, l3; // 3 lineas del bloque fastq 

			if ( !std::getline(file, l1) ){
				file.clear(); //limpiar flags de error antes de volver
				return false; 
			}
			if (l1.empty() || l1[0] != '@'){
				file.clear();
				file.seekg(pos_actual); //volvemos a la pos actual y q no es comienzo
				return false;
			}

			//leer las dos lineas mas para ver q es un + y comprobar q de verdad es el inicio 
			if ( !std::getline(file, l2) ){
				file.clear();
				file.seekg(pos_actual); //volvemos a la pos actual y q no es comienzo
				return false;
			}
			if ( !std::getline(file, l3) ){
				file.clear();
				file.seekg(pos_actual); //volvemos a la pos actual y q no es comienzo
				return false;
			}

			bool valid = (!l3.empty() && l3[0] == '+'); //la tercera linea debe empezar por + para ser un bloque fastq valido
			file.clear();
			file.seekg(pos_actual); //volvemos a la pos actual para q el proceso de lectura normal pueda leer el bloque completo
			return valid;
		}

		// parser paralelo de las lecturas, cada rank abrira se encargara de una parte del archivo, para ello claro, el archvivo debe estar descomprimido

		//para leer su trozo cada rank hago uso de seekg que me permite mover el puntero de lectura a la posicion que yo quiera.


		// funcion auxiliar pa buscar saltos de linea de una manera rapida en memoria ram (inline, pa decirle al compilador q en vez de la funcion meta su codigo en la llamada)
		inline char* find_next_newline( char* cursor, char* end ){
			return (char*) memchr( cursor, '\n', end - cursor ); //memchr busca el siguiente carater '\n' en el bloque de memoria entre cursor y end, me devuelve un puntero a ese caracter o nullptr si no lo encuentra
		}

		void parse_reads_parallel( SafeQueue< ReadBatches >& queue1, Stats& stats, Config const& config, uint64_t start1, uint64_t end1, uint64_t start2, uint64_t end2, int rank, bool interleaved, bool two_files){
			try{
				if (two_files) {
					// si es paired reads es q hay dos archivos 
					std::string filename1 = config.paired_reads[0];
					std::string filename2 = config.paired_reads[1];

					int fd1 = open(filename1.c_str(), O_RDONLY);
					int fd2 = open(filename2.c_str(), O_RDONLY);
					if (fd1 == -1 || fd2 == -1) {
						std::cerr << "[RANK " << rank << "] Error abriendo archivos Paired-End." << std::endl;
						if (fd1 != -1) close(fd1); if (fd2 != -1) close(fd2);
						return;
					}

					struct stat st1, st2; //sacar el tamano de los dos archivos 
					fstat(fd1, &st1); fstat(fd2, &st2);

					long page_size = sysconf(_SC_PAGESIZE); //tamano de pagina del sistema 
					size_t margin = 64 * 1024; //margen de seguridad para no quedar en medio de una lectura y no poder leer mas 

					// Mapear Archivo 1
					uint64_t map_offset1 = (start1 / page_size) * page_size;
					uint64_t local_skip1 = start1 - map_offset1;
					size_t map_length1 = (end1 - map_offset1) + margin;
					if (map_offset1 + map_length1 > st1.st_size) map_length1 = st1.st_size - map_offset1;
					char* map_ptr1 = static_cast<char*>(mmap(nullptr, map_length1, PROT_READ, MAP_PRIVATE, fd1, map_offset1));

					// Mapear Archivo 2
					uint64_t map_offset2 = (start2 / page_size) * page_size;
					uint64_t local_skip2 = start2 - map_offset2;
					size_t map_length2 = (end2 - map_offset2) + margin;
					if (map_offset2 + map_length2 > st2.st_size) map_length2 = st2.st_size - map_offset2;
					char* map_ptr2 = static_cast<char*>(mmap(nullptr, map_length2, PROT_READ, MAP_PRIVATE, fd2, map_offset2));

					madvise(map_ptr1, map_length1, MADV_SEQUENTIAL | MADV_WILLNEED);
					madvise(map_ptr2, map_length2, MADV_SEQUENTIAL | MADV_WILLNEED);

					char* cursor1 = map_ptr1 + local_skip1;
					char* end_ptr1 = map_ptr1 + map_length1;
					size_t logical_size1 = end1 - map_offset1;
					char* logical_end1 = map_ptr1 + std::min(logical_size1, map_length1);

					char* cursor2 = map_ptr2 + local_skip2;
					char* end_ptr2 = map_ptr2 + map_length2;
					size_t logical_size2 = end2 - map_offset2;
					char* logical_end2 = map_ptr2 + std::min(logical_size2, map_length2);

					// Sincronizar cursor 1
					if (start1 > 0) {
						char* nl = find_next_newline(cursor1, end_ptr1);
						if (nl) cursor1 = nl + 1;
						while (cursor1 < logical_end1) {
							char* at = (char*)memchr(cursor1, '@', end_ptr1 - cursor1);
							if (!at) break;
							char* l1 = find_next_newline(at, end_ptr1);
							if (l1) {
								char* l2 = find_next_newline(l1 + 1, end_ptr1);
								if (l2 && l2 + 1 < end_ptr1 && *(l2 + 1) == '+') { cursor1 = at; break; }
							}
							cursor1 = at + 1;
						}
					}

					// Sincronizar cursor 2
					if (start2 > 0) {
						char* nl = find_next_newline(cursor2, end_ptr2);
						if (nl) cursor2 = nl + 1;
						while (cursor2 < logical_end2) {
							char* at = (char*)memchr(cursor2, '@', end_ptr2 - cursor2);
							if (!at) break;
							char* l1 = find_next_newline(at, end_ptr2);
							if (l1) {
								char* l2 = find_next_newline(l1 + 1, end_ptr2);
								if (l2 && l2 + 1 < end_ptr2 && *(l2 + 1) == '+') { cursor2 = at; break; }
							}
							cursor2 = at + 1;
						}
					}

					size_t batch_size = config.n_reads;
					ReadBatches rb{ true }; // true porque es paired
					rb.ids.reserve(batch_size); rb.seqs.reserve(batch_size); rb.seqs2.reserve(batch_size);

					while (cursor1 < logical_end1 && cursor2 < logical_end2) {
						if (cursor1 >= end_ptr1 || cursor2 >= end_ptr2) break;

						// Extraer ID 1
						char* id1_start = cursor1 + 1;
						char* id1_end = find_next_newline(id1_start, end_ptr1);
						if (!id1_end) break;
						char* space1 = (char*)memchr(id1_start, ' ', id1_end - id1_start);
						size_t id1_len = space1 ? space1 - id1_start : id1_end - id1_start;

						// Extraer ID 2
						char* id2_start = cursor2 + 1;
						char* id2_end = find_next_newline(id2_start, end_ptr2);
						if (!id2_end) break;
						char* space2 = (char*)memchr(id2_start, ' ', id2_end - id2_start);
						size_t id2_len = space2 ? space2 - id2_start : id2_end - id2_start;

						rb.ids.emplace_back(id1_start, id1_len);

						// Secuencia 1
						char* seq1_start = id1_end + 1;
						char* seq1_end = find_next_newline(seq1_start, end_ptr1);
						std::string_view seq1_v(seq1_start, seq1_end - seq1_start);
						auto dna1 = seq1_v | seqan3::views::char_to<seqan3::dna4>;
						rb.seqs.emplace_back(dna1.begin(), dna1.end());

						char* qual1 = find_next_newline(find_next_newline(seq1_end + 1, end_ptr1) + 1, end_ptr1);
						cursor1 = qual1 + 1;

						// Secuencia 2
						char* seq2_start = id2_end + 1;
						char* seq2_end = find_next_newline(seq2_start, end_ptr2);
						std::string_view seq2_v(seq2_start, seq2_end - seq2_start);
						auto dna2 = seq2_v | seqan3::views::char_to<seqan3::dna4>;
						rb.seqs2.emplace_back(dna2.begin(), dna2.end());

						char* qual2 = find_next_newline(find_next_newline(seq2_end + 1, end_ptr2) + 1, end_ptr2);
						cursor2 = qual2 + 1;

						if (rb.ids.size() >= batch_size) {
							stats.input_reads += rb.ids.size();
							queue1.push(std::move(rb));
							rb = ReadBatches{ true };
							rb.ids.reserve(batch_size); rb.seqs.reserve(batch_size); rb.seqs2.reserve(batch_size);
						}
					}

					if (!rb.ids.empty()) { stats.input_reads += rb.ids.size(); queue1.push(std::move(rb)); }

					munmap(map_ptr1, map_length1); close(fd1);
					munmap(map_ptr2, map_length2); close(fd2);
					queue1.notify_push_over();

				} else {
					if (config.single_reads.empty()) return;
					std::string filename = config.single_reads[0];

					int fd = open(filename.c_str(), O_RDONLY);
					if (fd == -1) return;

					struct stat st;
					if (fstat(fd, &st) == -1) { close(fd); return; }

					uint64_t file_size = st.st_size;
					long page_size = sysconf(_SC_PAGESIZE);
					uint64_t map_offset = (start1 / page_size) * page_size;
					uint64_t local_skip = start1 - map_offset;
					size_t margin = 64 * 1024;
					size_t map_length = (end1 - map_offset) + margin;
					if (map_offset + map_length > file_size) map_length = file_size - map_offset;

					char* map_ptr = static_cast<char*>(mmap(nullptr, map_length, PROT_READ, MAP_PRIVATE, fd, map_offset));
					if (map_ptr == MAP_FAILED) { close(fd); return; }

					madvise(map_ptr, map_length, MADV_SEQUENTIAL | MADV_WILLNEED);

					char* cursor = map_ptr + local_skip;
					char* end_ptr = map_ptr + map_length;
					size_t logical_size = end1 - map_offset;
					char* logical_end = map_ptr + std::min(logical_size, map_length);

					if (start1 > 0) {
						char* newline = find_next_newline(cursor, end_ptr);
						if (newline) cursor = newline + 1;
						bool valido = false;
						while (cursor < logical_end) {
							char* at_ptr = (char*)memchr(cursor, '@', end_ptr - cursor);
							if (!at_ptr) break;
							char* l1_end = find_next_newline(at_ptr, end_ptr);
							if (l1_end) {
								char* l2_end = find_next_newline(l1_end + 1, end_ptr);
								if (l2_end && (l2_end + 1 < end_ptr) && *(l2_end + 1) == '+') {
									if (interleaved) {
										std::string_view id_view(at_ptr, l1_end - at_ptr);
										bool is_first = (id_view.find("/1") != std::string_view::npos || id_view.find(" 1:") != std::string_view::npos || id_view.find(" 1/") != std::string_view::npos);
										if (!is_first && (id_view.find("/2") != std::string_view::npos || id_view.find(" 2:") != std::string_view::npos || id_view.find(" 2/") != std::string_view::npos)) {
											cursor = l2_end + 1; continue;
										}
									}
									cursor = at_ptr; valido = true; break;
								}
							}
							cursor = at_ptr + 1;
						}
						if (!valido) { munmap(map_ptr, map_length); close(fd); queue1.notify_push_over(); return; }
					}

					size_t batch_size = config.n_reads;
					ReadBatches rb{ interleaved };
					rb.ids.reserve(batch_size); rb.seqs.reserve(batch_size);
					if (interleaved) rb.seqs2.reserve(batch_size);

					while (cursor < logical_end) {
						if (cursor >= end_ptr) break;
						char* id_start = cursor + 1;
						char* id_end = find_next_newline(id_start, end_ptr);
						if (!id_end) break;

						char* space_pos = (char*)memchr(id_start, ' ', id_end - id_start);
						size_t id_len = (space_pos) ? (space_pos - id_start) : (id_end - id_start);
						rb.ids.emplace_back(id_start, id_len);

						char* seq_start = id_end + 1;
						char* seq_end = find_next_newline(seq_start, end_ptr);
						if (!seq_end) break;
						std::string_view seq_view(seq_start, seq_end - seq_start);
						auto dna_view = seq_view | seqan3::views::char_to<seqan3::dna4>;
						rb.seqs.emplace_back(dna_view.begin(), dna_view.end());

						char* plus_line_end = find_next_newline(seq_end + 1, end_ptr);
						if (!plus_line_end) break;
						char* qual_line_end = find_next_newline(plus_line_end + 1, end_ptr);
						if (!qual_line_end) break;
						cursor = qual_line_end + 1;

						if (interleaved) {
							if (cursor >= end_ptr) break;
							char* id2_end = find_next_newline(cursor, end_ptr);
							if (!id2_end) break;
							char* seq2_start = id2_end + 1;
							char* seq2_end = find_next_newline(seq2_start, end_ptr);
							if (!seq2_end) break;

							std::string_view seq2_view(seq2_start, seq2_end - seq2_start);
							auto dna2_view = seq2_view | seqan3::views::char_to<seqan3::dna4>;
							rb.seqs2.emplace_back(dna2_view.begin(), dna2_view.end());

							char* plus2_line_end = find_next_newline(seq2_end + 1, end_ptr);
							if (!plus2_line_end) break;
							char* qual2_line_end = find_next_newline(plus2_line_end + 1, end_ptr);
							if (!qual2_line_end) break;
							cursor = qual2_line_end + 1;
						}

						if (rb.ids.size() >= batch_size) {
							stats.input_reads += rb.ids.size();
							queue1.push(std::move(rb));
							rb = ReadBatches{ interleaved };
							rb.ids.reserve(batch_size); rb.seqs.reserve(batch_size);
							if (interleaved) rb.seqs2.reserve(batch_size);
						}
					}
					if (!rb.ids.empty()) { stats.input_reads += rb.ids.size(); queue1.push(std::move(rb)); }

					queue1.notify_push_over();
					munmap(map_ptr, map_length); close(fd);
				}

			} catch (const std::exception& e){
				std::cerr << "\n[RANK " << rank << " FATAL ERROR EN HILO]: " << e.what() << "\n" << std::endl;
				queue1.notify_push_over();
			} catch(...){
				std::cerr << "\n[RANK " << rank << " FATAL ERROR EN HILO]: Excepción desconocida.\n" << std::endl;
				queue1.notify_push_over();
			}
		}


		/*
		   void parse_reads_parallel( SafeQueue< ReadBatches >& queue1, Stats& stats, Config const& config, uint64_t start, uint64_t end, int rank, bool interleaved){
		// funcion commo la otra parse reads pero si se le pasa un archivo con paired reads (interleaved a true) lo lee como paired_reads y si no como single reads. 
		// por que un archivo con las lecturas entremezcladas y no dos archivos con ellas separadas -> pues por q si le meto dos archivos y estoy usando bytes para inciios y finales aproximados , como la calidad de cada secuencia puede ser distinto, seria un rollo y con mal rendimiento ir sincronizando ambos archivos tratando de encontrar las secuencias que deberia ir juntas, es mas facil antes de ejecutar el programa pasarle un programa estilo seqtk y que las entremezcle una vez a uno con ambas lecturas juntas y listo. 

		try{

		if (config.single_reads.empty()){
		return; // si no hay reads pa fuera, lo pongo single por q antes de llamar a la funcion lo metere en single reads y ponde el interleaved a true o falsse
		}

		std::string filename = config.single_reads[0]; //nombre del archivo 

		//abrir el archivo y mapearlo 
		int fd = open(filename.c_str(), O_RDONLY); //abro el archivo en modo lectura, me da un descriptor de fichero,  .c_str() es necesario para convertir el string a un puntero a char que es lo q espera open

		if (fd == -1){
		std::cerr << "[RANK " << rank << "] Error opening file, open() fail: " << filename << std::endl; //mensaje de error si no de pude abirir el archivo 
		return;
		}

		struct stat st; //struct para almacenar la info del archivo 
		if (fstat ( fd, &st ) == -1){
		std::cerr << "[RANK " << rank << "] Error getting file size with fstat(): " << filename << std::endl; //mensaje de error si no de pude obtener el tamano del archivo 
		close(fd); //cerrar el descriptor de fichero antes de salir
		return;
		}

		uint64_t file_size = st.st_size; //tamanao del archivo en bytes
		long page_size = sysconf(_SC_PAGESIZE); //obtenemos el tamano de pagina del sistema con sysconf, esto es necesario para calcular la posicion de inicio alineada a pagina
		uint64_t map_offset = (start / page_size) * page_size; // calcular la posicion de incio de mapeao alineada, redondear hacia abajo claro, tiene q ser multiplo del tamano de pagina del sistema 
		uint64_t local_skip = start - map_offset; //calcular cuanto hay q saltar desde el inicio pa leer la parte que nos toca, saltar desde el inicio de la parte mapeada
		size_t margin = 64 * 1024; // 64kB de margen de seguridad para no leer fuera del archivo, si terminamos en medio de una secuencia poder acabar de leerla aunque nos saltemos de la parte q nos tocaba para no dejarla a medio leer
		size_t map_length = (end - map_offset) + margin; //calcular el tamano a mapear del archivo a memoria ram 
		if ( map_offset + map_length > file_size ){
		map_length = file_size - map_offset; //ajustar el tamano a mapear para no leer fuera del archivo 
		}

		char* map_ptr = static_cast<char*>( mmap( nullptr, map_length, PROT_READ, MAP_PRIVATE, fd, map_offset ) ); //mmap devuelve un puntero a la memoria mapeada, el primer parametro es nullptr para que el sistema elija la direccion de memoria, el segundo es el tamano a mapear, el tercero son las protecciones (lectura), el cuarto son las banderas (MAP_PRIVATE para que los cambios no se reflejen en el archivo), el quinto es el descriptor de fichero y el sexto es el offset del archivo a mapear 
		if ( map_ptr == MAP_FAILED ){
		std::cerr << "[RANK " << rank << "] Error mapping file to memory with mmap(): " << filename << std::endl; //mensaje de error si no de pude mapear el archivo a memoria ram 
		close(fd); //cerrar el descriptor de fichero antes de salir
		return;
		}

		madvise( map_ptr, map_length, MADV_SEQUENTIAL | MADV_WILLNEED ); //decirle al kernel que vamos a leer el archivo de manera secuencial y que se prepare para eso, esto puede mejorar el rendimiento al evitar lecturas innecesarias y optimizar el uso de la cache del sistema operativo 


		char* cursor = map_ptr + local_skip; //cursor va a estar donde comenzamos a leer q es el comienzo de la memoria mapeada mas lo q hay q saltarse
		char* end_ptr = map_ptr + map_length; //donde terminamos de leer incluyendo el margen
		size_t logical_size = end - map_offset; // tamano logico de la parque q nos toca leer 
		char* logical_end = map_ptr + std::min(logical_size, map_length); //puntero al final logico de nuestra parte, para evitar problemas de apuntar fuera del archivo si end - map_offset es mayor que map_length 
										  //char* logical_end = map_ptr + (end - map_offset);

										  //mirar el archivo 

										  if (start > 0){
		//saltar primera linea q estara cortada 
		char* newline = find_next_newline( cursor, end_ptr );
		if ( newline ) cursor = newline +1; //avanazr al inicio de la sig linea 
						    //buscar inicio de bloque valido
						    bool valido = false; // marcar si hay o no un bloque valido

						    while ( cursor < logical_end ){
		// buscar el @ 
		char* at_ptr = (char*) memchr( cursor, '@', end_ptr - cursor); // buscar el sig @ a partir del cursor 
		if ( !at_ptr ) break; // si no hay mas @ pa fuera 

		//verificar q no sea un @ aleatorio, comprobar q la tercera linea sea un + 
		char* l1_end = find_next_newline( at_ptr, end_ptr ); //buscar el final de la linea del @ 
		if(l1_end){
		char* l2_end = find_next_newline( l1_end + 1, end_ptr ); //buscar el final de la linea de la secuencia 
		if (l2_end){
			if ((l2_end +1 < end_ptr) && *(l2_end +1) == '+'){
				//es una lectura valida 

				//ahroa hay q mirar si es interleaved, si lo es, estamos en la primera lectura o en la segunda (la pareja)

				if ( interleaved ){ //herramienta pa meter 2 archivos de lecturas separadas -> https://github.com/linsalrob/fastq-pair
						    //
						    //si es interleaved hay q mirar si es la primera lectura o la segunda , la pareja vamos 
						    //mirar el ID a ver si va marcado que sea 1 o 2 
					std::string_view id_view ( at_ptr, l1_end - at_ptr); // string view del id 
											     //parsearlo a ver si tiene un 1 o un 2 
					bool is_first = false; // de momento no espero q sea la primera 
					if(id_view.find("/1") != std::string_view::npos || id_view.find(" 1:") != std::string_view::npos || id_view.find(" 1/") != std::string_view::npos){
						is_first = true; //es la primera lectura de la pareja 
					} 
					//si parece la segunda pues la salto y busco la siguiente 
					if (!is_first && ((id_view.find("/2") != std::string_view::npos) || (id_view.find(" 2:") != std::string_view::npos) || (id_view.find(" 2/") != std::string_view::npos))){
						//es la segunda lectura de la pareja, la salto y busco la siguiente 
						cursor = l2_end +1; //avanzar al inicio de la siguiente linea 
						continue; //saltar esta lectura y buscar la siguiente 
					}

				}
				// ahora el inicio es valido 
				cursor = at_ptr; //avanzar el cursor al inicio del bloque valido
				valido = true; //marcar q es inicio valido 
				break; //salir del bucle de busqueda de bloque valido
			}
		}
	}
	//avanzar cursor 
	cursor = at_ptr +1; //avanzar al siguiente caracter despues del @ para seguir buscando el siguiente bloque valido
	}
	if ( !valido ){
		//si no hay nada valido en nuestro trozo, pues terminar 
		munmap(	map_ptr, map_length ); //desmapear la memoria mapeada, esto es importante para liberar los recursos del sistema y evitar fugas de maemoria
		close(fd); //cerrar el descriptor de fichero, esto es importante para liberar los recursos del sistema y evitar fugas de descriptores de ficihero
		queue1.notify_push_over(); //notificar que ya no se van a enviar mas batches, esto es importante para que el hilo de procesamiento sepa cuando terminar y no se quede esperando indefinidamente en la cola 
		return; //terminar la funcion, no hay nada que procesar en nuestro bloque
	}

	} 

	// ahora toca parsear el bloque que dependera de si es single o paired claro 
	size_t batch_size = config.n_reads; // numero de lecturas por batch 

	//si es interleaved el batch tendria pares de lectura 

	ReadBatches rb{ interleaved }; //creamos un ReadBatches con el valor de interleaved para indicar si es single o paired 
	rb.ids.reserve(batch_size); //reservar memoria 
	rb.seqs.reserve(batch_size); //reservar memoria 
	if ( interleaved ){
		rb.seqs2.reserve(batch_size); //reservar memoria para la segunda secuencia de la pareja 
	}

	while ( cursor < logical_end ){
		if (cursor >= end_ptr) break; // archivo cortado no hay pa leer un id
					      //lectura 1 forward 
					      // ID 
		char* id_start = cursor +1; // saltar el @
		char* id_end = find_next_newline( id_start, end_ptr ); //buscar el final de la linea del ID 
		if ( !id_end ) break; //si no hay mas lineas pa fuera
				      //!
				      //el ID hasta el primer espacio 
		char* space_pos = (char*) memchr( id_start, ' ', id_end - id_start ); //buscar el primer espacio en la linea del ID


		size_t id_len = (space_pos) ? (space_pos - id_start) : (id_end - id_start); //calcular la longitud del ID hasta el primer espacio o hasta el final de la linea si no hay espacio

		rb.ids.emplace_back( id_start, id_len ); //construir el string del ID directamente desde el puntero y la longitud, esto evita copias innecesarias y mejora el rendimiento

		// secuecnia 

		char* seq_start = id_end +1; //inicio de la secuencia, es la siguiente linea despues del ID 
		char* seq_end = find_next_newline( seq_start, end_ptr ); //buscar el final de la linea de la secuencia 
		if ( !seq_end ) break; //si no hay mas lineas pa fuera 

		std::string_view seq_view( seq_start, seq_end - seq_start ); //crear un string_view de la secuencia directamente desde el puntero y la longitud, esto evita copias innecesarias y mejora el rendimiento

		auto dna_view = seq_view | seqan3::views::char_to< seqan3::dna4 >; //convertir el string_view a una vista de dna4 usando char_to, esto evita copias innecesarias y mejora el rendimiento

		rb.seqs.emplace_back( dna_view.begin(), dna_view.end() ); //construir el vector de dna4 directamente desde la vista, esto evita copias innecesarias y mejora el rendimiento
									  // skil + calidad 
		char* plus_line_end = find_next_newline( seq_end +1, end_ptr ); //buscar el final de la linea del +

		if ( !plus_line_end ) break; //si no hay mas lineas pa fuera 
		char* qual_line_end = find_next_newline( plus_line_end +1, end_ptr ); //buscar el final de la linea de la calidad 
		if ( !qual_line_end ) break; //si no hay mas lineas pa fuera

		cursor = qual_line_end +1; //avanzar el cursor al inicio de la siguiente lectura 

		//si es interleaved toca leer la pareja
		if ( interleaved ){
			if (cursor >= end_ptr) break; //archivo cortado/ impar raro 
						      //ignorar el id, es el mismo pero con un 1 o un 2 
			char* id2_end = find_next_newline( cursor, end_ptr ); //buscar el final de la linea del ID de la pareja 
			if ( !id2_end ) break; //si no hay mas lineas pa fuera 

			//seq2 
			char* seq2_start = id2_end +1; //inicio de la secuencia de la pareja, es la siguiente linea despues del ID 
			char* seq2_end = find_next_newline( seq2_start, end_ptr ); //buscar el final de la linea de la secuencia de la pareja 

			if ( !seq2_end ) break; //si no hay mas lineas pa fuera 

			std::string_view seq2_view( seq2_start, seq2_end - seq2_start ); //crear un string_view de la secuencia de la pareja directamente desde el puntero y la longitud, esto evita copias innecesarias y mejora el rendimiento
			auto dna2_view = seq2_view | seqan3::views::char_to< seqan3::dna4 >; //convertir el string_view a una vista de dna4 usando char_to, esto evita copias innecesarias y mejora el rendimiento
			rb.seqs2.emplace_back( dna2_view.begin(), dna2_view.end() ); //construir el vector de dna4 de la pareja directamente desde la vista, esto evita copias innecesarias y mejora el rendimiento

			char* plus2_line_end = find_next_newline( seq2_end +1, end_ptr ); //buscar el final de la linea del + de la pareja
			if ( !plus2_line_end ) break; //si no hay mas lineas pa fuera
			char* qual2_line_end = find_next_newline( plus2_line_end +1, end_ptr ); //buscar el final de la linea de la calidad de la pareja
			if ( !qual2_line_end ) break; //si no hay mas lineas pa fuera
			cursor = qual2_line_end +1; //avanzar el cursor al inicio de la siguiente lectura, que seria la siguiente pareja
		}
		// enviar batch 
		if ( rb.ids.size() >= batch_size){
			stats.input_reads += rb.ids.size(); //actualizamos el contador de lecturas procesadas con el lote actualizamos
			queue1.push( std::move(rb) ) ; //enviar el batch a la cola, moviendo el ReadBatches para evitar copias innecesarias y mejorar el rendimiento
			rb = ReadBatches{ interleaved }; //crear un nuevo ReadBatches para el siguiente lote, con el mismo valor de interleaved, esto evita copias innecesarias y mejora el rendimiento
			rb.ids.reserve(batch_size); //reservar memoria para el nuevo lote 
			rb.seqs.reserve(batch_size); //reservar memoria para el nuevo lote 
			if ( interleaved ){
				rb.seqs2.reserve(batch_size); //reservar memoria para la segunda secuencia de la pareja en el nuevo lote
			}
		}

	}
	// ultimo push 
	if ( !rb.ids.empty() ){
		stats.input_reads += rb.ids.size(); //actualizamos el contador de lecturas procesadas con el lote actualizamos
		queue1.push( std::move(rb) ) ; //enviar el batch a la cola, moviendo el ReadBatches para evitar copias innecesarias y mejorar el rendimiento
	}
	queue1.notify_push_over(); //notificar que ya no se van a enviar mas batches, esto es importante para que el hilo de procesamiento sepa cuando terminar y no se quede esperando indefinidamente en la cola 
	munmap( map_ptr, map_length ); //desmapear la memoria mapeada, esto es importante para liberar los recursos del sistema y evitar fugas de memoria 
	close(fd); //cerrar el descriptor de fichero, esto es importante para liberar los recursos del sistema y evitar fugas de descriptores de fichero



	} catch (const std::exception& e){
		// ATRAPAR EXCEPCIONES PARA QUE NO SE HUNDAN EN SILENCIO
		std::cerr << "\n[RANK " << rank << " FATAL ERROR EN HILO]: " << e.what() << "\n" << std::endl;
		queue1.notify_push_over(); // Evitar deadlock
	} catch(...){
		std::cerr << "\n[RANK " << rank << " FATAL ERROR EN HILO]: Excepción desconocida (OOM/Segfault).\n" << std::endl;
		queue1.notify_push_over(); // Evitar deadlock
	}

	}
	*/

		/*
		   void parse_reads_parallel( SafeQueue< ReadBatches >& queue1 , Stats& stats, Config const& config, uint64_t start, uint64_t end, int rank){
		   try{
		   long long lecturas_procesadas = 0;
		//de momento implementado para simple reads para paired reads habria que sincronizar ambos archivos
		if ( config.single_reads.empty() ){
		return; //si no hay single reads no hay nada q hacer
		}

		std::string filename = config.single_reads[0]; //nombre del archivo 

		//abrir el archivo 
		int fd = open(filename.c_str(), O_RDONLY); //abro el archivo en modo lectura, me da un descriptor de fichero,  .c_str() es necesario para convertir el string a un puntero a char que es lo q espera open 

		std::ifstream file; //stream de lectura para leer el archivo 

		if (fd == -1){
		std::cerr << "[RANK " << rank << "] Error opening file, open() fail: " << filename << std::endl; //mensaje de error si no de pude abirir el archivo 
		return;
		}

		//obtener tamano real 
		struct stat st; //struct para almacenar la info del archivo 

		if (fstat ( fd, &st ) == -1){
		std::cerr << "[RANK " << rank << "] Error getting file size with fstat(): " << filename << std::endl; //mensaje de error si no de pude obtener el tamano del archivo 
		close(fd); //cerrar el descriptor de fichero antes de salir
		return;
		}

		uint64_t file_size = st.st_size; //tamanao del archivo en bytes 

		//calcular la alineacion de pagina, a mmap hay q pasarle un multiplo del tamano de pagina pa q funcione, tamano de pagina suele ser 4kB si varia , la solucion es sacarla de la confi del sistema 

		long page_size = sysconf(_SC_PAGESIZE); //obtenemos el tamano de pagina del sistema con sysconf, esto es necesario para calcular la posicion de inicio alineada a pagina

		uint64_t map_offset = (start / page_size) * page_size; // calcular la posicion de incio de mapeao alineada, redondear hacia abajo claro 

		uint64_t local_skip = start - map_offset; //calcular cuanto hay q saltar desde el inicio pa leer la parte que nos toca 

		// ahora poner un margen de seguridad para no leer fuera de memoria por si el final que nos toca leer nos quedamos a medias de una lectura 
		// con 64 kB es mas que suficiente 
		size_t margin = 64 * 1024; // 64kB 
		size_t map_length = (end - map_offset) + margin; //calcular el tamano a mapear del archivo a memoria ram 

		if ( map_offset + map_length > file_size ){
		map_length = file_size - map_offset; //ajustar el tamano a mapear para no leer fuera del archivo 
		}

		//mapear el archivo a memoria RAM 
		char* map_ptr = static_cast<char*>( mmap( nullptr, map_length, PROT_READ, MAP_PRIVATE, fd, map_offset ) ); //mmap devuelve un puntero a la memoria mapeada, el primer parametro es nullptr para que el sistema elija la direccion de memoria, el segundo es el tamano a mapear, el tercero son las protecciones (lectura), el cuarto son las banderas (MAP_PRIVATE para que los cambios no se reflejen en el archivo), el quinto es el descriptor de fichero y el sexto es el offset del archivo a mapear

		if ( map_ptr == MAP_FAILED ){
		std::cerr << "[RANK " << rank << "] Error mapping file to memory with mmap(): " << filename << std::endl; //mensaje de error si no de pude mapear el archivo a memoria ram 
		close(fd); //cerrar el descriptor de fichero antes de salir
		return;
		}

		//ahora para optimizar mas , como vamos a leer todo secuencialmente sin ir a saltos (la memoria mapeada ) le vamos a decir al kernel que se preprare, que va a mapear agresivamente, que meta muchos datos a memoria no vaya poco a poco que no nos vamos a saltar nada si no que vamos a leer todo de corrido y secuencialmente sin saltos 

		if ( madvise( map_ptr, map_length, MADV_SEQUENTIAL | MADV_WILLNEED ) == -1 ){
		std::cerr << "[RANK " << rank << "] Warning: madvise() failed: " << filename << std::endl; //mensaje de advertencia si no se pudo usar madvise, no es un error crítico pero puede afectar al rendimiento 
		} // mas info en https://man7.org/linux/man-pages/man2/madvise.2.html


		// punteros utiles 

		char* cursor = map_ptr + local_skip; //puntero al inicio de la parte que nos toca leer, esto es lo que usaremos para leer el archivo 
		char* end_ptr = map_ptr + map_length; //puntero al final de la parte que nos toca leer, esto es lo que usaremos para no leer fuera de nuestra parte asignad
		char* logical_end = map_ptr + (end - map_offset); //puntero al final logico de nuestra parte, esto es para saber cuando hemos llegado al final de nuestra parte asignada aunque el mapa sea un poco mas grande por el margen de seguridad

		if( start > 0 ){
			// si no empezamos al principio, buscar la primera sig linea 

			//toca buscar el siguiente encabezado (@) para empezar, simplemente toca avabzar
			while ( cursor < end_ptr && *cursor != '\n' )
				cursor++; //avanzar hasta el final de la linea actual
			if ( cursor < end_ptr){
				cursor++; //avanzar al inicio de la siguiente linea
			}
			//buscar el patron de inici '@'
			while ( cursor < logical_end){
				if (*cursor == '@')
					break; //inicio de linea encontrado 
				while (cursor < end_ptr && *cursor != '\n')
					cursor++; //avanzar hasta el final de la linea actual
				if (cursor < end_ptr){
					cursor++; //avanzar al inicio de la siguiente linea
				}
			}
		}	

	//reservar memoria para evitar que el vector se redimensione constantemente y disminuya el rendimiento 
	size_t batch_size = config.n_reads; // numero de lecturas por batch 

	//toca leer hasta final de bloque o mas para no cortar ninguna lectura 
	std::vector< std::string > ids;
	std::vector< std::vector< seqan3::dna4 > > seqs;

	ids.reserve(batch_size); 
	seqs.reserve(batch_size);

	//parsear el bloque 
	while ( cursor < logical_end){
		if ( cursor >= end_ptr)
			break; //llegamos al final 
			       //leer el ID 
		if (*cursor != '@'){
			while (cursor < end_ptr && *cursor != '\n')
				cursor++; //avanzar hasta el final de la linea actual
			cursor++; continue; 
		}

		//saltar el '@'
		cursor++;

		char* id_start = cursor; //inicio del id_start 
		while ( cursor < end_ptr && *cursor != '\n')
			cursor++; //avanzar hasta el final de la linea actual 
				  //copiar el id al vector 
		ids.emplace_back( id_start, cursor - id_start ); //construir el string directamente desde el puntero y la longitud, esto evita copias innecesarias y mejora el rendimiento
		cursor++; //avanzar, saltar el '\n'
			  //leer la secuencialmente
			  //
		char* seq_start = cursor; //inicio de la secuencia 
		while ( cursor < end_ptr && *cursor != '\n')
			cursor++; //avanzar hasta el final de la linea actual
		size_t seq_len = cursor - seq_start; //longitud de la secuencialmente

		//usare string_view para no copiar la cadena si no pasarle el puntero directamente a seqan3 
		std::string_view seq_view( seq_start, seq_len ); //crear un string_view directamente desde el puntero y la longitud, esto evita copias innecesarias y mejora el rendimiento
		auto dna_view = seq_view | seqan3::views::char_to< seqan3::dna4 >; //convertir el string_view a una vista de dna4 usando char_to, esto evita copias innecesarias y mejora el rendimiento
		seqs.emplace_back( dna_view.begin(), dna_view.end() ); //construir el vector de dna4 directamente desde la vista, esto evita copias innecesarias y mejora el rendimiento
		cursor++; //avanzar, saltar el '\n'

		// el +, la linea del separador 
		while ( cursor < end_ptr && *cursor != '\n')
			cursor++; //avanzar hasta el final de la linea actual

		//la calidad 
		while ( cursor < end_ptr && *cursor != '\n')
			cursor++; //avanzar hasta el final de la linea actual 
		cursor++; 

		//enviar el batch 

		if ( ids.size() >= batch_size ){
			stats.input_reads += ids.size(); //actualizamos el contador de lecturas procesadas con el lote actual
			queue1.push( ReadBatches( false, std::move(ids), std::move(seqs) ) ) ; //enviamos el lote a la cola, moviendo los vectores para evitar copias y mejorar el rendimiento
			ids.clear(); //limpiar el vector de ids para el siguiente lote 
			seqs.clear(); //limpiar el vector de secuencias para el siguiente lote
		}
	}


	// ahora toca enviar el ultimo lote que no esta completo 
	if ( !ids.empty() ){
		stats.input_reads += ids.size(); //actualizamos el contador de lecturas procesadas con el ultimo lote
		queue1.push( ReadBatches( false, std::move(ids), std::move(seqs) ) ) ; //enviamos el ultimo lote a la cola
	}
	//notificar a al consumidor (cola) q no hay na mas 
	queue1.notify_push_over(); // no se van a meter mas lotes

	//limpiar, demapear la memoria 
	if ( munmap( map_ptr, map_length ) == -1 ){
		std::cerr << "[RANK " << rank << "] Warning: munmap() failed: " << filename << std::endl; //mensaje de advertencia si no se pudo demapear la memoria, no es un error crítico pero puede afectar al rendimiento 
	}
	close(fd); //cerrar el descriptor de fichero

	} catch (const std::exception& e){
		// ATRAPAR EXCEPCIONES PARA QUE NO SE HUNDAN EN SILENCIO
		std::cerr << "\n[RANK " << rank << " FATAL ERROR EN HILO]: " << e.what() << "\n" << std::endl;
		queue1.notify_push_over(); // Evitar deadlock
	} catch(...){
		std::cerr << "\n[RANK " << rank << " FATAL ERROR EN HILO]: Excepción desconocida (OOM/Segfault).\n" << std::endl;
		queue1.notify_push_over(); // Evitar deadlock
	}
	}

	*/

		std::map< std::string, HierarchyConfig > parse_hierarchy( Config& config )
		{

			std::map< std::string, HierarchyConfig > parsed_hierarchy;

			std::vector< std::string > sorted_hierarchy = config.hierarchy_labels;
			std::sort( sorted_hierarchy.begin(), sorted_hierarchy.end() );
			// get unique hierarcy labels
			const size_t unique_hierarchy =
				std::unique( sorted_hierarchy.begin(), sorted_hierarchy.end() ) - sorted_hierarchy.begin();

			size_t hierarchy_count = 0;
			for ( size_t h = 0; h < config.hierarchy_labels.size(); ++h )
			{

				auto filter_cfg = FilterConfig{ config.ibf[h], config.rel_cutoff[h] };
				if ( config.tax.size() > 0 )
					filter_cfg.tax_file = config.tax[h];

				if ( parsed_hierarchy.find( config.hierarchy_labels[h] ) == parsed_hierarchy.end() )
				{ // not found
				  // validate by hiearchy
					std::vector< FilterConfig > fc;
					fc.push_back( filter_cfg );
					std::string output_file_lca = "";
					std::string output_file_all = "";
					if ( !config.output_prefix.empty() && unique_hierarchy > 1 && !config.output_single )
					{
						output_file_lca = config.output_prefix + "." + config.hierarchy_labels[h] + ".one";
						output_file_all = config.output_prefix + "." + config.hierarchy_labels[h] + ".all";
					}
					else if ( !config.output_prefix.empty() )
					{
						output_file_lca = config.output_prefix + ".one";
						output_file_all = config.output_prefix + ".all";
					}

					parsed_hierarchy[config.hierarchy_labels[h]] = HierarchyConfig{ fc,
						0,
						0,
						config.rel_filter[hierarchy_count],
						config.fpr_query[hierarchy_count],
						output_file_lca,
						output_file_all };
					++hierarchy_count;
				}
				else
				{ // found
					parsed_hierarchy[config.hierarchy_labels[h]].filters.push_back( filter_cfg );
				}
			}

			return parsed_hierarchy;
		}

	void print_hierarchy( Config const& config, auto const& parsed_hierarchy )
	{

		constexpr auto newl{ "\n" };
		for ( auto const& hierarchy_config : parsed_hierarchy )
		{
			std::cerr << hierarchy_config.first << newl;
			std::cerr << "--rel-filter " << hierarchy_config.second.rel_filter << newl;
			std::cerr << "--fpr-query " << hierarchy_config.second.fpr_query << newl;
			for ( auto const& filter_config : hierarchy_config.second.filters )
			{
				std::cerr << "    " << filter_config.ibf_file;
				if ( !filter_config.tax_file.empty() )
					std::cerr << ", " << filter_config.tax_file;
				if ( filter_config.rel_cutoff > -1 )
					std::cerr << " --rel-cutoff " << filter_config.rel_cutoff;
				std::cerr << newl;
			}
			if ( !config.output_prefix.empty() )
			{
				std::cerr << "    Output files: ";
				std::cerr << config.output_prefix + ".rep";
				if ( config.output_lca )
					std::cerr << ", " << hierarchy_config.second.output_file_lca;
				if ( config.output_all )
					std::cerr << ", " << hierarchy_config.second.output_file_all;
				std::cerr << newl;
			}
		}
		std::cerr << "----------------------------------------------------------------------" << newl;
	}

	inline TRep sum_reports( std::vector< TRep > const& reports )
	{
		TRep report_sum;
		for ( auto const& report : reports )
		{
			for ( auto const& [target, r] : report )
			{
				report_sum[target].matches += r.matches;
				report_sum[target].lca_reads += r.lca_reads;
				report_sum[target].unique_reads += r.unique_reads;
			}
		}
		return report_sum;
	}

	inline size_t threshold_rel( size_t n_hashes, double p )
	{
		return std::ceil( n_hashes * p );
	}

	// https://stackoverflow.com/questions/44718971/calculate-binomial-coffeficient-very-reliably
	inline double binom( double n, double k ) noexcept
	{
		return std::exp( std::lgamma( n + 1 ) - std::lgamma( n - k + 1 ) - std::lgamma( k + 1 ) );
	}


	void select_matches( Filter< TIBF >&        filter,
			TMatches&              matches,
			std::vector< size_t >& hashes,
			auto&                  agent,
			size_t                 threshold_cutoff,
			size_t&                max_count_read,
			size_t&                min_count_read,
			size_t                 n_hashes )
	{
		// Count every occurrence on IBF
		seqan3::counting_vector< detail::TIntCount > counts = agent.bulk_count( hashes );

		for ( auto const& [target, bins] : filter.map )
		{
			// Sum counts among bins (split target (user bins) into several technical bins)
			size_t summed_count = 0;
			for ( auto const& binno : bins )
			{
				summed_count += counts[binno];
			}
			// summed_count can be higher than n_hashes for matches in several technical bins
			if ( summed_count > n_hashes )
				summed_count = n_hashes;
			if ( summed_count >= threshold_cutoff )
			{
				// ensure that count was not already found for target with higher count
				// can happen in case of ambiguos targets in multiple filters
				if ( summed_count > std::get< 0 >( matches[target] ) )
				{
					matches[target] = std::make_tuple( summed_count, filter.filter_config.target_fpr[target] );
					if ( summed_count > max_count_read )
						max_count_read = summed_count;
					if ( summed_count < min_count_read )
						min_count_read = summed_count;
				}
			}
		}
	}

	void select_matches( Filter< THIBF >&       filter,
			TMatches&              matches,
			std::vector< size_t >& hashes,
			auto&                  agent,
			size_t                 threshold_cutoff,
			size_t&                max_count_read,
			size_t&                min_count_read,
			size_t                 n_hashes )
	{
		// Count only matches above threhsold
		seqan3::counting_vector< detail::TIntCount > counts = agent.bulk_count( hashes, threshold_cutoff );

		// Only one bin per target
		for ( auto const& [target, bins] : filter.map )
		{
			if ( counts[bins[0]] > 0 )
			{
				// Sum counts among bins (split target (user bins) into several technical bins)
				size_t summed_count = counts[bins[0]];
				// summed_count can be higher than n_hashes for matches in several technical bins
				if ( summed_count > n_hashes )
					summed_count = n_hashes;
				// ensure that count was not already found for target with higher count
				// can happen in case of ambiguous targets in multiple filters
				if ( summed_count > std::get< 0 >( matches[target] ) )
				{
					matches[target] = std::make_tuple( summed_count, filter.filter_config.target_fpr[target] );
					if ( summed_count > max_count_read )
						max_count_read = summed_count;
					if ( summed_count < min_count_read )
						min_count_read = summed_count;
				}
			}
		}
	}

	size_t filter_matches(
			ReadOut& read_out, TMatches& matches, TRep& rep, size_t n_hashes, double threshold_filter, double min_fpr_query )
	{

		for ( auto const& [target, count_fpr] : matches )
		{
			if ( std::get< 0 >( count_fpr ) >= threshold_filter )
			{
				// Filter by fpr-query
				if ( min_fpr_query < 1.0 )
				{
					double q = 1;
					for ( size_t i = 0; i <= std::get< 0 >( count_fpr ); i++ )
					{
						q -= binom( n_hashes, i ) * pow( std::get< 1 >( count_fpr ), i )
							* pow( 1 - std::get< 1 >( count_fpr ), n_hashes - i );
					}
					if ( q > min_fpr_query )
					{
						continue;
					}
				}

				rep[target].matches++;
				read_out.matches.push_back( ReadMatch{ target, std::get< 0 >( count_fpr ) } );
			}
		}

		return read_out.matches.size();
	}

	void lca_matches( ReadOut& read_out, ReadOut& read_out_lca, size_t max_count_read, LCA& lca, TRep& rep )
	{

		std::vector< std::string > targets;
		for ( auto const& r : read_out.matches )
		{
			targets.push_back( r.target );
		}

		std::string target_lca = lca.getLCA( targets );
		rep[target_lca].lca_reads++;
		read_out_lca.matches.push_back( ReadMatch{ target_lca, max_count_read } );
	}


	template < typename TFilter >
		void classify( std::vector< Filter< TFilter > >& filters,
				LCA&                              lca,
				TRep&                             rep,
				Total&                            total,
				SafeQueue< ReadOut >&             classified_all_queue,
				SafeQueue< ReadOut >&             classified_lca_queue,
				SafeQueue< ReadOut >&             unclassified_queue,
				Config const&                     config,
				SafeQueue< ReadBatches >*         pointer_current,
				SafeQueue< ReadBatches >*         pointer_helper,
				HierarchyConfig const&            hierarchy_config,
				bool                              hierarchy_first,
				bool                              hierarchy_last )
		{

			// oner hash adaptor per thread
			const auto minimiser_hash =
				seqan3::views::minimiser_hash( seqan3::shape{ seqan3::ungapped{ hierarchy_config.kmer_size } },
						seqan3::window_size{ hierarchy_config.window_size },
						seqan3::seed{ raptor::adjust_seed( hierarchy_config.kmer_size ) } );

			// one agent per thread per filter
			using TAgent = std::conditional_t< std::same_as< TFilter, THIBF >,
			      THIBF::counting_agent_type< detail::TIntCount >,
			      TIBF::counting_agent_type< detail::TIntCount > >;
			std::vector< TAgent > agents;
			for ( Filter< TFilter >& filter : filters )
			{
				agents.push_back( filter.ibf.template counting_agent< detail::TIntCount >() );
			}

			while ( true )
			{
				// Wait here until reads are available or push is over and queue is empty
				ReadBatches rb = pointer_current->pop();

				// If batch is empty exit thread
				if ( rb.ids.empty() )
					break;

				// store unclassified reads for next iteration
				ReadBatches left_over_reads{ rb.paired };

				const size_t hashes_limit = std::numeric_limits< detail::TIntCount >::max();

				for ( size_t readID = 0; readID < rb.ids.size(); ++readID )
				{
					// read lenghts
					const size_t read1_len = rb.seqs[readID].size();
					const size_t read2_len = rb.paired ? rb.seqs2[readID].size() : 0;

					// Store matches for this read
					TMatches matches;

					// Best scoring kmer count
					size_t max_count_read = 0;
					size_t min_count_read = 0;
					size_t n_hashes       = 0;
					// if length is smaller than window, skip read
					if ( read1_len >= hierarchy_config.window_size )
					{
						// Count hashes
						std::vector< size_t > hashes = rb.seqs[readID] | minimiser_hash | seqan3::ranges::to< std::vector >();
						// Count hashes from both pairs if second is given
						if ( read2_len >= hierarchy_config.window_size )
						{
							// Add hashes of second pair
							const auto h2 = rb.seqs2[readID] | minimiser_hash | std::views::common;
							hashes.insert( hashes.end(), h2.begin(), h2.end() );
						}

						n_hashes = hashes.size();
						// set min as max. possible hashes
						min_count_read = n_hashes;
						// if n_hashes are bigger than int limit, skip read
						if ( n_hashes <= hashes_limit )
						{
							// Sum sequence to totals
							if ( hierarchy_first )
							{
								total.reads_processed++;
								total.length_processed += read1_len + read2_len;
							}

							// For each filter in the hierarchy
							for ( size_t i = 0; i < filters.size(); ++i )
							{
								// Calculate threshold for cutoff (keep matches above)
								size_t threshold_cutoff = threshold_rel( n_hashes, filters[i].filter_config.rel_cutoff );

								// reset low threshold_cutoff to just one kmer (0 would match everywhere)
								if ( threshold_cutoff == 0 )
									threshold_cutoff = 1;

								// count and select matches
								select_matches( filters[i],
										matches,
										hashes,
										agents[i],
										threshold_cutoff,
										max_count_read,
										min_count_read,
										n_hashes );
							}
						}
					}

					// store read and matches to be printed
					ReadOut read_out( rb.ids[readID] );

					// if read got valid matches (above cutoff)
					if ( max_count_read > 0 )
					{

						// Calculate threshold for filtering (keep matches above)
						const size_t threshold_filter =
							max_count_read - threshold_rel( max_count_read - min_count_read, hierarchy_config.rel_filter );

						// Filter matches
						const size_t count_filtered_matches =
							filter_matches( read_out, matches, rep, n_hashes, threshold_filter, hierarchy_config.fpr_query );

						if ( count_filtered_matches > 0 )
						{

							total.reads_classified++;

							if ( !config.skip_lca )
							{
								ReadOut read_out_lca( rb.ids[readID] );
								if ( count_filtered_matches == 1 )
								{
									// just one match, copy read read_out and set as unique
									read_out_lca = read_out;
									rep[read_out.matches[0].target].unique_reads++;
								}
								else
								{
									lca_matches( read_out, read_out_lca, max_count_read, lca, rep );
								}

								if ( config.output_lca )
									classified_lca_queue.push( read_out_lca );
							}
							else
							{
								// Not running lca and has unique match
								if ( count_filtered_matches == 1 )
								{
									rep[read_out.matches[0].target].unique_reads++;
								}
								else
								{
									// without tax, no lca, count multi-matches to a root node
									// to keep consistency among reports (no. of classified reads)
									rep[config.tax_root_node].lca_reads++;
								}
							}

							if ( config.output_all )
								classified_all_queue.push( read_out );

							// read classified, continue to the next
							continue;
						}
					}

					// not classified
					if ( !hierarchy_last ) // if there is more levels, store read
					{
						left_over_reads.ids.push_back( std::move( rb.ids[readID] ) );
						left_over_reads.seqs.push_back( std::move( rb.seqs[readID] ) );

						if ( rb.paired )
						{
							// seqan::appendValue( left_over_reads.seqs2, rb.seqs2[readID] );
							left_over_reads.seqs2.push_back( std::move( rb.seqs2[readID] ) );
						}
					}
					else if ( config.output_unclassified ) // no more levels and no classification, add to
									       // unclassified printing queue
					{
						unclassified_queue.push( read_out );
					}
				}

				// if there are more levels to classify and something was left, keep reads in memory
				if ( !hierarchy_last && left_over_reads.ids.size() > 0 )
					pointer_helper->push( std::move( left_over_reads ) );
			}
		}

	void write_report( TRep& rep, TTax& tax, std::ofstream& out_rep, std::string hierarchy_label )
	{
		for ( auto const& [target, report] : rep )
		{
			if ( report.matches || report.lca_reads || report.unique_reads )
			{
				out_rep << hierarchy_label << '\t' << target << '\t' << report.matches << '\t' << report.unique_reads
					<< '\t' << report.lca_reads;

				if ( !tax.empty() )
				{
					out_rep << '\t' << tax.at( target ).rank << '\t' << tax.at( target ).name;
				}
				out_rep << '\n';
			}
		}
	}

	static inline void replace_all( std::string& str, const std::string& from, const std::string& to )
	{
		size_t start_pos = 0;
		while ( ( start_pos = str.find( from, start_pos ) ) != std::string::npos )
		{
			str.replace( start_pos, from.length(), to );
			start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
		}
	}

	size_t load_filter( THIBF&             filter,
			IBFConfig&         ibf_config,
			TBinMap&           bin_map,
			std::string const& input_filter_file,
			TTargetFpr&        target_fpr )
	{
		std::ifstream              is( input_filter_file, std::ios::binary );
		cereal::BinaryInputArchive archive( is );

		uint32_t                                  parsed_version;
		uint64_t                                  window_size;
		seqan3::shape                             shape{};
		uint8_t                                   parts;
		bool                                      compressed;
		std::vector< std::vector< std::string > > bin_path{};
		double                                    fpr{};
		bool                                      is_hibf{};

		archive( parsed_version );
		archive( window_size );
		archive( shape );
		archive( parts );
		archive( compressed );
		archive( bin_path );
		archive( fpr );
		archive( is_hibf );
		archive( filter );

		// load ibf_config from raptor params
		ibf_config.window_size = window_size;
		ibf_config.kmer_size   = shape.count();
		ibf_config.max_fp      = fpr;

		// Create map from paths
		size_t binno{};
		for ( auto const& file_list : bin_path )
		{
			for ( auto const& filename : file_list )
			{
				// based on the filename, get target.minimiser
				// (e.g. 562.minimiser or GCF_013391805.1.minimiser), otherwise use filename as target
				auto   f     = std::filesystem::path( filename ).filename().string();
				size_t found = f.find( ".minimiser" );
				if ( found != std::string::npos )
					f = f.substr( 0, found );

				// workaround when file has a . (e.g. GCF_013391805.1)
				// "." replaced by "|||" in ganon build wrapper
				// fixed on ganon v2.0.0 (+ raptor 3.0.1) but kept for compatibility (>= ganon v1.8.0)
				replace_all( f, "|||", "." );

				// workaround when target has a space (e.g. s__Pectobacterium carotovorum)
				// " " replaced by "---" in ganon build wrapper
				replace_all( f, "---", " " );

				bin_map.push_back( std::make_tuple( binno, f ) );
				// same fpr for all
				target_fpr[f] = fpr;
			}
			++binno;
		}

		return filter.user_bins.num_user_bins();
	}

	inline double false_positive( uint64_t bin_size_bits, uint8_t hash_functions, uint64_t n_hashes )
	{
		/*
		 * calculates the theoretical false positive of a bin (bf) based on parameters
		 */
		return std::pow( 1 - std::exp( -hash_functions / ( bin_size_bits / static_cast< double >( n_hashes ) ) ),
				hash_functions );
	}

	size_t load_filter( TIBF&              filter,
			IBFConfig&         ibf_config,
			TBinMap&           bin_map,
			std::string const& input_filter_file,
			TTargetFpr&        target_fpr )
	{
		std::ifstream              is( input_filter_file, std::ios::binary );
		cereal::BinaryInputArchive archive( is );

		std::tuple< int, int, int >                        parsed_version;
		std::vector< std::tuple< std::string, uint64_t > > hashes_count_std;

		archive( parsed_version );
		archive( ibf_config );
		archive( hashes_count_std );
		archive( bin_map );
		archive( filter );


		// generate fpr for each bin
		for ( auto const& [target, count] : hashes_count_std )
		{
			// Use average number of hashes for each bin to calculate fp
			uint64_t n_bins_target = std::ceil( count / static_cast< double >( ibf_config.max_hashes_bin ) );
			// this can be off by a very small number (rounding ceil on multiple bins)
			uint64_t n_hashes_bin = std::ceil( count / static_cast< double >( n_bins_target ) );

			// false positive for the current target, considering split bins
			target_fpr[target] =
				1.0
				- std::pow( 1.0 - false_positive( ibf_config.bin_size_bits, ibf_config.hash_functions, n_hashes_bin ),
						n_bins_target );
			;
		}


		return filter.bin_count();
	}

	TTax load_tax( std::string tax_file )
	{
		TTax          tax;
		std::string   line;
		std::ifstream infile;
		infile.open( tax_file );
		while ( std::getline( infile, line, '\n' ) )
		{
			std::istringstream         stream_line( line );
			std::vector< std::string > fields;
			std::string                field;
			while ( std::getline( stream_line, field, '\t' ) )
				fields.push_back( field );
			tax[fields[0]] = Node{ fields[1], fields[2], fields[3] };
		}
		infile.close();
		return tax;
	}

	template < typename TFilter >
		bool load_files( std::vector< Filter< TFilter > >& filters, std::vector< FilterConfig >& fconf )
		{
			size_t filter_cnt = 0;
			for ( auto& filter_config : fconf )
			{
				TTax       tax;
				IBFConfig  ibf_config;
				TBinMap    bin_map;
				TFilter    filter;
				TTargetFpr target_fpr;
				auto       bin_count = load_filter( filter, ibf_config, bin_map, filter_config.ibf_file, target_fpr );

				// Parse vector with bin_map to the old map
				TMap map;
				for ( auto const& [binno, target] : bin_map )
				{
					map[target].push_back( binno );
				}

				filter_config.ibf_config = ibf_config;
				filter_config.target_fpr = target_fpr;

				if ( filter_config.tax_file != "" )
					tax = load_tax( filter_config.tax_file );

				filters.push_back(
						Filter< TFilter >{ std::move( filter ), std::move( map ), std::move( tax ), bin_count, filter_config } );
				filter_cnt++;
			}

			return true;
		}

	void print_time( const StopClock& timeGanon, const StopClock& timeLoadFilters, const StopClock& timeClassPrint )
	{
		using ::operator<<;
		std::cerr << "ganon-classify        start time: " << StopClock_datetime( timeGanon.begin() ) << std::endl;
		std::cerr << "loading filters      elapsed (s): " << timeLoadFilters.elapsed() << " seconds" << std::endl;
		std::cerr << "classifying+printing elapsed (s): " << timeClassPrint.elapsed() << " seconds" << std::endl;
		std::cerr << "ganon-classify       elapsed (s): " << timeGanon.elapsed() << " seconds" << std::endl;
		std::cerr << "ganon-classify          end time: " << StopClock_datetime( timeGanon.end() ) << std::endl;
		std::cerr << std::endl;
	}

	void print_stats( Stats& stats, const StopClock& timeClassPrint, auto const& parsed_hierarchy )
	{
		const double elapsed_classification = timeClassPrint.elapsed();
		const double total_reads_processed  = stats.total.reads_processed > 0
			? static_cast< double >( stats.total.reads_processed )
			: 1; // to not report nan on divisions
		std::cerr << "ganon-classify processed " << stats.total.reads_processed << " sequences ("
			<< stats.total.length_processed / 1000000.0 << " Mbp) in " << elapsed_classification << " seconds ("
			<< ( stats.total.length_processed / 1000000.0 ) / ( elapsed_classification / 60.0 ) << " Mbp/m)"
			<< std::endl;
		std::cerr << " - " << stats.total.reads_classified << " reads classified ("
			<< ( stats.total.reads_classified / total_reads_processed ) * 100 << "%)" << std::endl;
		std::cerr << "   - " << stats.total.unique_matches << " with unique matches ("
			<< ( stats.total.unique_matches / total_reads_processed ) * 100 << "%)" << std::endl;
		std::cerr << "   - " << stats.total.reads_classified - stats.total.unique_matches << " with multiple matches ("
			<< ( ( stats.total.reads_classified - stats.total.unique_matches ) / total_reads_processed ) * 100 << "%)"
			<< std::endl;

		double avg_matches = stats.total.reads_classified
			? ( stats.total.matches / static_cast< double >( stats.total.reads_classified ) )
			: 0;
		std::cerr << " - " << stats.total.matches << " matches (avg. " << avg_matches << " match/read classified)"
			<< std::endl;
		const size_t total_reads_unclassified = stats.total.reads_processed - stats.total.reads_classified;
		std::cerr << " - " << total_reads_unclassified << " reads unclassified ("
			<< ( total_reads_unclassified / total_reads_processed ) * 100 << "%)" << std::endl;

		if ( stats.total.reads_processed < stats.input_reads )
		{
			std::cerr << " - " << stats.input_reads - stats.total.reads_processed
				<< " reads skipped (too long or too short (< window size))" << std::endl;
		}

		if ( parsed_hierarchy.size() > 1 )
		{
			std::cerr << std::endl;
			std::cerr << "By database hierarchy:" << std::endl;
			for ( auto const& h : parsed_hierarchy )
			{
				std::string hierarchy_label = h.first;
				avg_matches                 = stats.hierarchy_total[hierarchy_label].reads_classified
					? ( stats.hierarchy_total[hierarchy_label].matches
							/ static_cast< double >( stats.hierarchy_total[hierarchy_label].reads_classified ) )
					: 0;
				std::cerr << " - " << hierarchy_label << ": " << stats.hierarchy_total[hierarchy_label].reads_classified
					<< " classified ("
					<< ( stats.hierarchy_total[hierarchy_label].reads_classified / total_reads_processed ) * 100
					<< "%) " << stats.hierarchy_total[hierarchy_label].unique_matches << " unique ("
					<< ( stats.hierarchy_total[hierarchy_label].unique_matches / total_reads_processed ) * 100
					<< "%) "
					<< stats.hierarchy_total[hierarchy_label].reads_classified
					- stats.hierarchy_total[hierarchy_label].unique_matches
					<< " multiple ("
					<< ( ( stats.hierarchy_total[hierarchy_label].reads_classified
								- stats.hierarchy_total[hierarchy_label].unique_matches )
							/ total_reads_processed )
					* 100
					<< "%) " << stats.hierarchy_total[hierarchy_label].matches << " matches (avg. " << avg_matches
					<< ")" << std::endl;
			}
		}
	}

	void parse_reads( SafeQueue< ReadBatches >& queue1, Stats& stats, Config const& config )
	{
		for ( auto const& reads_file : config.single_reads )
		{
			try
			{
				seqan3::sequence_file_input< raptor::dna4_traits, seqan3::fields< seqan3::field::id, seqan3::field::seq > >
					fin1{ reads_file };
				for ( auto&& rec : fin1 | seqan3::views::chunk( config.n_reads ) )
				{
					ReadBatches rb{ false };
					for ( auto& [id, seq] : rec )
					{
						rb.ids.push_back( std::move( id ) );
						rb.seqs.push_back( std::move( seq ) );
					}
					stats.input_reads += rb.ids.size();
					queue1.push( std::move( rb ) );
				}
			}
			catch ( seqan3::parse_error const& e )
			{
				std::cerr << "Error parsing file [" << reads_file << "]. " << e.what() << std::endl;
				continue;
			}
		}
		if ( config.paired_reads.size() > 0 )
		{
			for ( size_t pair_cnt = 0; pair_cnt < config.paired_reads.size(); pair_cnt += 2 )
			{
				try
				{
					seqan3::sequence_file_input< raptor::dna4_traits,
						seqan3::fields< seqan3::field::id, seqan3::field::seq > >
							fin1{ config.paired_reads[pair_cnt] };
					seqan3::sequence_file_input< raptor::dna4_traits,
						seqan3::fields< seqan3::field::id, seqan3::field::seq > >
							fin2{ config.paired_reads[pair_cnt + 1] };
					for ( auto&& rec : fin1 | seqan3::views::chunk( config.n_reads ) )
					{
						ReadBatches rb{ true };
						for ( auto& [id, seq] : rec )
						{
							rb.ids.push_back( std::move( id ) );
							rb.seqs.push_back( std::move( seq ) );
						}
						// loop in the second file and get same amount of reads
						for ( auto& [id, seq] : fin2 | std::views::take( config.n_reads ) )
						{
							rb.seqs2.push_back( std::move( seq ) );
						}
						stats.input_reads += rb.ids.size();
						queue1.push( std::move( rb ) );
					}
				}
				catch ( seqan3::parse_error const& ext )
				{
					std::cerr << "Error parsing files [" << config.paired_reads[pair_cnt] << "/"
						<< config.paired_reads[pair_cnt + 1] << "]. " << ext.what() << std::endl;
					continue;
				}
			}
		}
		queue1.notify_push_over();
	}

	void write_classified( SafeQueue< ReadOut >& classified_queue, std::ofstream& out )
	{
		while ( true )
		{
			ReadOut ro = classified_queue.pop();
			if ( ro.readID != "" )
			{
				for ( size_t i = 0; i < ro.matches.size(); ++i )
				{
					out << ro.readID << '\t' << ro.matches[i].target << '\t' << ro.matches[i].kmer_count << '\n';
				}
			}
			else
			{
				break;
			}
		}
	}

	void write_unclassified( SafeQueue< ReadOut >& unclassified_queue, std::string out_unclassified_file )
	{
		std::ofstream out_unclassified( out_unclassified_file );
		while ( true )
		{
			ReadOut rou = unclassified_queue.pop();
			if ( rou.readID != "" )
			{
				out_unclassified << rou.readID << '\n';
			}
			else
			{
				out_unclassified.close();
				break;
			}
		}
	}

	template < typename TFilter >
		TTax merge_tax( std::vector< Filter< TFilter > > const& filters )
		{
			if ( filters.size() == 1 )
			{
				return filters[0].tax;
			}
			else
			{
				TTax merged_tax = filters[0].tax;
				for ( size_t i = 1; i < filters.size(); ++i )
				{
					// merge taxonomies keeping the first one as a default
					merged_tax.insert( filters[i].tax.begin(), filters[i].tax.end() );
				}
				return merged_tax;
			}
		}

	template < typename TFilter >
		void validate_targets_tax( std::vector< Filter< TFilter > > const& filters,
				TTax&                                   tax,
				bool                                    quiet,
				const std::string                       tax_root_node )
		{
			for ( auto const& filter : filters )
			{
				for ( auto const& [target, bins] : filter.map )
				{
					if ( tax.count( target ) == 0 )
					{
						tax[target] = Node{ tax_root_node, "no rank", target };
						if ( !quiet )
							std::cerr << "WARNING: target [" << target << "] without tax entry, setting parent as root node ["
								<< tax_root_node << "]" << std::endl;
					}
				}
			}
		}

	void pre_process_lca( LCA& lca, TTax& tax, std::string tax_root_node )
	{
		for ( auto const& [target, node] : tax )
		{
			lca.addEdge( node.parent, target );
		}
		lca.doEulerWalk( tax_root_node );
	}

	} // namespace detail

	template < typename TFilter >
		bool ganon_classify( Config config , int rank, int numeroProcesos)
		{
			// Start time count
			StopClock timeGanon;
			timeGanon.start();

			auto parsed_hierarchy = detail::parse_hierarchy( config );

			if ( config.verbose && rank == 0)
				detail::print_hierarchy( config, parsed_hierarchy );

			// Initialize variables
			StopClock timeLoadFilters;
			StopClock timeClassPrint;

			detail::Stats stats;
			std::ofstream out_rep; // Set default output stream (file or stdout)
			std::ofstream out_all; // output all file
			std::ofstream out_lca; // output lca file

			if (rank == 0){
				// If there's no output prefix, redirect to STDOUT
				if ( config.output_prefix.empty() )
				{
					out_rep.copyfmt( std::cout ); // STDOUT
					out_rep.clear( std::cout.rdstate() );
					out_rep.basic_ios< char >::rdbuf( std::cout.rdbuf() );
				}
				else
				{
					out_rep.open( config.output_prefix + ".rep" );
				}
			}

			// calcular rangos de lectura 
			uint64_t start_read1 = 0;
			uint64_t end_read1 = 0;

			uint64_t start_read2 = 0;
			uint64_t end_read2 = 0;

			bool interleaved = false;

			bool two_files = false; 

			// comprobamos si son paired reads y 2 archivos, si es paired reads pero un archivo (interleaved) o si es single reads. 
			if ( config.paired_reads.size() == 2 ){
				two_files = true;
				try{
					//calcular chunks pa el archivo 1 y 2 , los calculo a parte aunque deberian ser iguales ya que tienen q tener mismas lecturas y mismos numero de bytes si no fallara ya que no seria posible de una manera correcta, rapida (eficinte) y simple sincronizar cada archivo en cada rank 
					uint64_t file_size1 = std::filesystem::file_size(config.paired_reads[0]);
					uint64_t file_size2 = std::filesystem::file_size(config.paired_reads[1]);
					if (file_size1 != file_size2){
						std::cerr << "Error: paired read files have different sizes, cannot be processed in parallel" << std::endl;
						return false;
					}
					uint64_t chunk_size = file_size1 / numeroProcesos; //cuanto le va a tocar a cada proceso 
											  //caculo el resto 
											  //uint64_t resto = file_size % numeroProcesos;
											  // no hace falta calcularlo, el ultimo proceso se queda con el resto y listo, poca diferencia va a haber en la practica y evito andar con un bucle para añadirle unos bytes a los primeros x procesos 
					start_read1 = rank * chunk_size;
					end_read1 = (rank == numeroProcesos -1) ? file_size1 : (rank+1) * chunk_size; // el ultimo proceso se queda con el resto 
					start_read2 = rank * chunk_size;
					end_read2 = (rank == numeroProcesos -1) ? file_size2 : (rank+1) * chunk_size; // el ultimo proceso se queda con el resto 

				} catch ( const std::exception& e ){
					std::cerr << "Error getting file size: " << e.what() << std::endl;
					return false;
				}
			}else if ( config.paired_reads.size() == 1 ) {
				//interleaved 
				interleaved = true;
				try{
					uint64_t file_size = std::filesystem::file_size(config.paired_reads[0]);
					uint64_t chunk_size = file_size / numeroProcesos; //cuanto le va a tocar a cada proceso 
											  //caculo el resto 
											  //uint64_t resto = file_size % numeroProcesos;
											  // no hace falta calcularlo, el ultimo proceso se queda con el resto y listo, poca diferencia va a haber en la practica y evito andar con un bucle para añadirle unos bytes a los primeros x procesos 
					start_read1 = rank * chunk_size;
					end_read1 = (rank == numeroProcesos -1) ? file_size : (rank+1) * chunk_size; // el ultimo proceso se queda con el resto 
				} catch ( const std::exception& e ){
					std::cerr << "Error getting file size: " << e.what() << std::endl;
					return false;	
				}
			}
			//si es single reads 
			if (!two_files && !config.single_reads.empty()){
				try{
					uint64_t file_size = std::filesystem::file_size(config.single_reads[0]);
					uint64_t chunk_size = file_size / numeroProcesos; //cuanto le va a tocar a cada proceso 
											  //caculo el resto 
											  //uint64_t resto = file_size % numeroProcesos;
											  // no hace falta calcularlo, el ultimo proceso se queda con el resto y listo, poca diferencia va a haber en la practica y evito andar con un bucle para añadirle unos bytes a los primeros x procesos 
					start_read1 = rank * chunk_size; 
					end_read1 = (rank == numeroProcesos -1) ? file_size : (rank+1) * chunk_size; // el ultimo proceso se queda con el resto
				} catch ( const std::exception& e ){
					std::cerr << "Error getting file size: " << e.what() << std::endl;
					return false;
				}
			}

			/*
			if ( !config.single_reads.empty() ){
				std::string filename = config.single_reads[0];
				try{
					uint64_t file_size = std::filesystem::file_size(filename);
					uint64_t chunk_size = file_size / numeroProcesos; //cuanto le va a tocar a cada proceso 
											  //caculo el resto 
											  //uint64_t resto = file_size % numeroProcesos;
											  // no hace falta calcularlo, el ultimo proceso se queda con el resto y listo, poca diferencia va a haber en la practica y evito andar con un bucle para añadirle unos bytes a los primeros x procesos 
					start_read1 = rank * chunk_size; 
					end_read1 = (rank == numeroProcesos -1) ? file_size : (rank+1) * chunk_size; // el ultimo proceso se queda con el resto
				} catch ( const std::exception& e ){
					std::cerr << "Error getting file size: " << e.what() << std::endl;
					return false;
				}
			}*/

			// Queues for internal read handling
			// queue1 get reads from file
			// queue2 will get unclassified reads if hierachy == 2
			// if hierachy == 3 queue1 is used for unclassified and so on
			// config.n_batches*config.n_reads = max. amount of reads in memory
			SafeQueue< detail::ReadBatches >  queue1( config.n_batches );
			SafeQueue< detail::ReadBatches >  queue2;
			SafeQueue< detail::ReadBatches >* pointer_current = &queue1; // pointer to the queues
			SafeQueue< detail::ReadBatches >* pointer_helper  = &queue2; // pointer to the queues
			SafeQueue< detail::ReadBatches >* pointer_extra;             // pointer to the queues

			// Define one threads for decompress bgzf files
			seqan3::contrib::bgzf_thread_count = 1u;

			/*
			// mirar si es paired o no  
			if ( !config.paired_reads.empty() ){
				//es un solo archivo inteleaved, lo meto en single reads pero interleaved a true 
				config.single_reads.push_back(config.paired_reads[0]);
				config.paired_reads.clear();
				interleaved = true;
			}*/

			// debug para saber q ha llegado
			std::cerr << "[RANK " << rank << "] Empezando a leer archivo)..." << std::endl;			

			// Thread for reading input files
			std::future< void > read_task = std::async(
					std::launch::async, detail::parse_reads_parallel, std::ref( queue1 ), std::ref( stats ), std::ref( config ), start_read1, end_read1, start_read2, end_read2, rank, interleaved, two_files ); //lanzar hilo de lecturas en paralelo

			// Thread for printing unclassified reads
			SafeQueue< detail::ReadOut > unclassified_queue;
			std::future< void >          write_unclassified_task;
			if ( config.output_unclassified && !config.output_prefix.empty() )
			{
				write_unclassified_task = std::async( std::launch::async,
						detail::write_unclassified,
						std::ref( unclassified_queue ),
						config.output_prefix + ".unc" );
			}


			// Classify reads iteractively for each hierarchy level
			size_t       hierarchy_id   = 0;
			const size_t hierarchy_size = parsed_hierarchy.size();
			for ( auto& [hierarchy_label, hierarchy_config] : parsed_hierarchy )
			{
				++hierarchy_id;
				bool                                     hierarchy_first = ( hierarchy_id == 1 );
				bool                                     hierarchy_last  = ( hierarchy_id == hierarchy_size );
				std::vector< detail::Filter< TFilter > > filters;
				detail::TTax                             tax;
				LCA                                      lca;

				timeLoadFilters.start();
				const bool loaded = detail::load_files( filters, parsed_hierarchy[hierarchy_label].filters );
				if ( !loaded )
				{
					std::cerr << "ERROR: loading ibf or tax files" << std::endl;
					return false;
				}
				timeLoadFilters.stop();

				hierarchy_config.kmer_size   = filters[0].filter_config.ibf_config.kmer_size;
				hierarchy_config.window_size = filters[0].filter_config.ibf_config.window_size;
				if ( filters.size() > 1 )
				{
					// Check if all filters share the same k,w
					for ( auto const& filter : filters )
					{
						if ( filter.filter_config.ibf_config.kmer_size != hierarchy_config.kmer_size
								|| filter.filter_config.ibf_config.window_size != hierarchy_config.window_size )
						{
							std::cerr << "ERROR: databases on the same hierarchy should share same k-mer and window sizes"
								<< std::endl;
							return false;
						}
					}
				}


				// Parse tax if provided
				if ( filters[0].filter_config.tax_file != "" )
				{
					// merge repeated elements from multiple filters
					tax = detail::merge_tax( filters );
					// if target not found in tax, add node target with parent root
					detail::validate_targets_tax( filters, tax, config.quiet, config.tax_root_node );
				}

				if ( !config.skip_lca )
				{
					if ( tax.count( config.tax_root_node ) == 0 )
					{
						std::cerr << "Root node [" << config.tax_root_node << "] not found (--tax-root-node)" << std::endl;
						return false;
					}
					// pre-processing of nodes to generate LCA
					detail::pre_process_lca( lca, tax, config.tax_root_node );
				}


				// Thread for printing classified reads (.lca, .all)
				std::vector< std::future< void > > write_tasks;

				// hierarchy_id = 1
				//  pointer_current=queue1, data comes from file in a limited size queue
				//  pointer_helper=queue2, empty
				// hierarchy_id > 1
				//  pointer_current=queue2, with all data already in from last iteration
				//  pointer_helper=queue1, empty
				// Exchange queues instance pointers for each hierachy (if not first)
				if ( !hierarchy_first )
				{
					pointer_extra   = pointer_current;
					pointer_current = pointer_helper;
					pointer_helper  = pointer_extra;

					// Remove size limit from reading since it's always already loaded
					if ( hierarchy_id == 2 )
						queue1.set_max_size( std::numeric_limits< size_t >::max() );
				}

				SafeQueue< detail::ReadOut > classified_all_queue;
				SafeQueue< detail::ReadOut > classified_lca_queue;

				if ( !config.output_prefix.empty() )
				{
					// nombres de ficihero de salida con el rank del procesos
					std::string hierarchy_output_file_lca = hierarchy_config.output_file_lca;
					std::string hierarchy_output_file_all = hierarchy_config.output_file_all;
					if ( rank > 0 ){
						if ( hierarchy_config.output_file_lca.empty() == false){
							hierarchy_output_file_lca += ".rank" + std::to_string(rank);
						}
						if ( hierarchy_config.output_file_all.empty() == false){
							hierarchy_output_file_all += ".rank" + std::to_string(rank);
						}
					}
					if ( config.output_lca && !config.skip_lca )
					{
						if ( hierarchy_first || !config.output_single )
							out_lca.open( hierarchy_output_file_lca );
						else // append if not first and output_single
							out_lca.open( hierarchy_output_file_lca, std::ofstream::app );

						// Start writing thread for lca matches
						write_tasks.emplace_back( std::async( std::launch::async,
									detail::write_classified,
									std::ref( classified_lca_queue ),
									std::ref( out_lca ) ) );
					}
					if ( config.output_all )
					{
						if ( hierarchy_first || !config.output_single )
							out_all.open( hierarchy_output_file_all );
						else // append if not first and output_single
							out_all.open( hierarchy_output_file_all, std::ofstream::app );

						// Start writing thread for all matches
						write_tasks.emplace_back( std::async( std::launch::async,
									detail::write_classified,
									std::ref( classified_all_queue ),
									std::ref( out_all ) ) );
					}
				}

				// One report and total counters for each thread
				std::vector< detail::TRep >  reports( config.threads );
				std::vector< detail::Total > totals( config.threads );

				std::vector< std::future< void > > tasks;
				// Threads for classification
				timeClassPrint.start();

				std::cerr << "[RANK " << rank << "] Empezando a lanzar a los hilos clasificadores (tasks)..." << std::endl;

				for ( size_t taskNo = 0; taskNo < config.threads; ++taskNo )
				{

					tasks.emplace_back( std::async( std::launch::async,
								detail::classify< TFilter >,
								std::ref( filters ),
								std::ref( lca ),
								std::ref( reports[taskNo] ),
								std::ref( totals[taskNo] ),
								std::ref( classified_all_queue ),
								std::ref( classified_lca_queue ),
								std::ref( unclassified_queue ),
								std::ref( config ),
								pointer_current,
								pointer_helper,
								hierarchy_config,
								hierarchy_first,
								hierarchy_last ));
				}

				// Wait here until classification is over
				//

				std::cerr << "[RANK " << rank << "] Empezando a esperar a los workers (tasks)..." << std::endl;

				for ( auto&& task : tasks )
				{
					task.get();
				}

				std::cerr << "[RANK " << rank << "] Fin de  esperar a los workers (tasks)..." << std::endl;

				// After classification, no more reads are going to be pushed to the output
				classified_all_queue.notify_push_over();
				classified_lca_queue.notify_push_over();



				// Sum reports of each threads into one
				detail::TRep rep = sum_reports( reports ); //cada thread tiene su propio reporte, juntarolo todo

				// Sum totals for each thread and report into stats
				stats.add_totals( hierarchy_label, totals );
				//stats.add_reports( hierarchy_label, rep ); los reportes pa el final 

				//mandar al rank 0 los resultados 
				size_t local_classified = stats.total.reads_classified;
				size_t local_matches = stats.total.matches;
				size_t local_unique = stats.total.unique_matches;
				size_t local_processed = stats.total.reads_processed;
				size_t local_input = stats.input_reads;

				size_t global_classified = 0;
				size_t global_matches = 0;
				size_t global_unique = 0;
				size_t global_processed = 0;
				size_t global_input = 0;

				std::cerr << "[RANK " << rank << "] Empezando a esperar a hacer los reduce..." << std::endl;

				//hacer los MPI reduce para sumar los resultados
				MPI_Reduce(&local_classified, &global_classified, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(&local_matches, &global_matches, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(&local_unique, &global_unique, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(&local_processed, &global_processed, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(&local_input, &global_input, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

				//una vez tenemos el reduce 
				if ( rank == 0 ){
					stats.total.reads_classified = global_classified;
					stats.total.matches = 0; //ya estan en el mapa 
					stats.total.unique_matches = 0;
					stats.total.reads_processed = global_processed;
					stats.input_reads = global_input;
					// limpiar los contadores por jerarquia 
					for ( auto& [etiqueta, j] : stats.hierarchy_total){
						j.matches = 0;
						j.unique_matches = 0;
					}
				}


				//unificar el mapa rep usando cereal para serializar y deserializar el mapa y enviarlo al rank 0. Podria intentar usar Mpi pack y unpack y un objetivo creado para el mapa pero el mapa es una coleccion de punteros distribuidos en memoria y por tanto no estan en memoria con un desplazamiento fijo y se vuelve todo muy complejo y complicado
				std::cerr << "[RANK " << rank << "] Empezando a serializar el mapa con cereal..." << std::endl;

				// seralizar y enviar los mapas a los ranks 
				// para hacerlo mas rapido que secuencial lo voy a hacer estilo arbol binario, asi seran pasos logaritmicos y no lineales. 
				int paso = 1; //indicar el paso en el que estamos 

				while ( paso < numeroProcesos ){
					// recibe datos y los suma a los suyos 
					if (rank % ( 2 * paso) == 0){
						// si el rank es multiplo de 2*paso recibe los datos del rank rank+paso y los suma a los suyos 
						int source = rank + paso;
						if ( source < numeroProcesos ){

							//recibir el tamano 
							int len = 0;
							MPI_Recv(&len, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							//recibir el mapa serializado (los datos )
							std::cerr << "[DEBUG RANK " << rank << "] - Recibiendo la cabecera (tamano) del rank " << source << " . Bytes a leer: " << len << std::endl;
							std::vector<char> buffer(len);
							MPI_Recv(buffer.data(), len, MPI_CHAR, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							//deserializar el mapa recibido
							std::string serialized_str(buffer.data(), len);
							std::stringstream ss_in(serialized_str);
							//std::vector<std::pair<std::string, detail::Rep> > temp_vec;
							std::vector<std::string> keys;
							std::vector<detail::Rep > values;
							try{
								cereal::BinaryInputArchive iarchive(ss_in);
								//iarchive(temp_vec);
								iarchive(keys);
								iarchive(values);
							}catch ( std::exception& e ){
								std::cerr << "[ERROR RANK " << rank << " ]. Error al deserializar CERAL: " << e.what() << std::endl;
							}

							//std::cerr << "[DEBUG RANK " << rank << " ]. Desempaquetados " << temp_vec.size() << " pares del rank " << source << std::endl;

							size_t cont_nuevos = 0;
							size_t matches_recibidos_total = 0;

							/*
							// fusionar en el mapa local sumando contadores 
							for ( auto const& [target, rep_received] : temp_vec ){
							if ( rep.find(target) == rep.end() ) cont_nuevos++;
							if ( target.empty() ) continue; //si esta vacio siguiente
							detail::Rep& rep_local = rep[target];
							rep_local.matches += rep_received.matches;
							rep_local.lca_reads += rep_received.lca_reads;
							rep_local.unique_reads += rep_received.unique_reads;
							}*/

							for ( size_t k = 0; k < keys.size(); ++k ){
								auto& target = keys[k];
								auto& r_recv = values[k];

								matches_recibidos_total += r_recv.matches;
								std::cerr << "[DEBUG RANK " << rank << "] Sumando al mapa local, matches recibidos: " << r_recv.matches << " , lca reads: " << r_recv.lca_reads << " , unique reads: " << r_recv.unique_reads << std::endl;
								//sumar al mapa local 
								rep[target].matches += r_recv.matches;
								rep[target].lca_reads += r_recv.lca_reads;
								rep[target].unique_reads += r_recv.unique_reads;

							}

							//std::cerr << "[DEBUG RANK " << rank << " ]. Fusion completada, anadidos " << cont_nuevos << " virus nuevos. Total actual es: " << rep.size() << std::endl;

							std::cerr << "[DEBUG RANK " << rank << "] Recibidos " << keys.size() << " virus del Rank " << source << ". Suma de matches recibidos: " << matches_recibidos_total << std::endl;
						}

					}else{
						// si el rank no es multiplo, es el que envia los datos 
						int dest = rank - paso;

						size_t num_elementos = 0;
						size_t total_matches_enviados = 0;

						//serializar 
						std::stringstream ss;

						try{
							cereal::BinaryOutputArchive oarchive(ss);
							//std::vector<std::pair<std::string, detail::Rep> > temp_vec;
							std::vector<std::string> keys;
							std::vector<detail::Rep> values;
							keys.reserve(rep.size());
							values.reserve(rep.size());
							//temp_vec.reserve(rep.size());


							for ( auto const& [key,val] : rep ) {
								keys.push_back(key);
								values.push_back(val);
								total_matches_enviados += val.matches;

							}

							oarchive(keys);
							oarchive(values);


							/*
							   for ( auto const& kv : rep ){
							   temp_vec.emplace_back(kv.first, kv.second);
							   }
							   */

							//num_elementos = temp_vec.size();
							//oarchive(temp_vec);
						}catch ( std::exception& e ){
							std::cerr << "[ ERROR RANK " << rank << "]" << rank << " Fallo al serializar CEREAL: " << e.what() << std::endl;
						}

						std::string serialized_str = ss.str();
						int len = static_cast<int>(serialized_str.size());

						//std::cerr << "[ERROR RANK " << rank << "] -- > Enviando " << num_elementos << " elementos ( " << len << " bytes) al rank " << dest << std::endl;

						std::cerr << "[DEBUG RANK " << rank << "] Enviando " << rep.size() << " virus al Rank " << dest << ". Suma de matches a enviar: " << total_matches_enviados << std::endl;

						//enviar el tamano y los datos con un Isend y un Waitall (no bloqueantes) 
						MPI_Request requests[2]; // array de request para el Isend, el primero pa el tamano y el segundo pa los datos 
						MPI_Isend(&len, 1, MPI_INT, dest, 0, MPI_COMM_WORLD, &requests[0]);
						MPI_Isend(serialized_str.c_str(), len, MPI_CHAR, dest, 0, MPI_COMM_WORLD, &requests[1]);
						MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
						// una vez enviado el mapa, salir del bucle 
						break;
					}
					paso *= 2; //pasar al siguiente paso (multiplicar por 2)
				}


				if (rank == 0){
					std::cerr << "[RANK " << rank << "] Rank 0 anadir al reporte y escribir el reporte ..." << std::endl;
					stats.add_reports( hierarchy_label, rep );

					// write reports
					detail::write_report( rep, tax, out_rep, hierarchy_label );
				}

				std::cerr << "[RANK " << rank << "] Empezando a esperar a los hilos de escritura (tasks)..." << std::endl;

				// Wait here until all files are written
				for ( auto&& task : write_tasks )
				{
					task.get();
				}

				std::cerr << "[RANK " << rank << "] FIN de esperar a los workers (tasks)..." << std::endl;
				timeClassPrint.stop();

				std::cerr << "[RANK " << rank << "] Cerrar archivos escritura)..." << std::endl;
				// Close file for writing (if not STDOUT)
				if ( !config.output_prefix.empty() )
				{
					if ( config.output_lca && out_lca.is_open())
						out_lca.close();
					if ( config.output_all && out_all.is_open())
						out_all.close();
				}

				std::cerr << "[RANK " << rank << "] Fin cerrar archivos de escritura..." << std::endl;

				// ahora todos los procesos esperan 
				if ( hierarchy_first )
				{
					read_task.get();                    // get reading tasks at the end of the first hierarchy
					pointer_helper->notify_push_over(); // notify push is over, only on first time (will be always set over for
									    // next iterations)
				}

			}

			// Wait here until all unclassified reads are written
			if ( config.output_unclassified )
			{
				unclassified_queue.notify_push_over();
				write_unclassified_task.get();
			}

			if ( rank == 0 ) {
				out_rep << "#total_classified\t" << stats.total.reads_classified << '\n';
				// account for unclassified and skipped sequences
				out_rep << "#total_unclassified\t" << stats.input_reads - stats.total.reads_classified << '\n';
				if ( !config.output_prefix.empty() )
				{
					out_rep.close();
				}
			}

			timeGanon.stop();

			if ( rank == 0 ) {
				if ( !config.quiet )
				{
					std::cerr << std::endl;
					if ( config.verbose )
					{
						detail::print_time( timeGanon, timeLoadFilters, timeClassPrint );
					}
					detail::print_stats( stats, timeClassPrint, parsed_hierarchy );
				}
			}
			return true;
		}

	bool run( Config config, int numeroProcesos, int rank )
	{

		// Validate configuration input
		if ( !config.validate() )
			return false;

		// Print config
		if ( config.verbose )
			std::cerr << config;

		if ( config.hibf )
			return ganon_classify< detail::THIBF >( config, rank, numeroProcesos );
		else
			return ganon_classify< detail::TIBF >( config , rank, numeroProcesos);
	}

} // namespace GanonClassify
