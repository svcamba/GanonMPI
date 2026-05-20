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

#include <cinttypes>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

//añadir MPI 
#include <mpi.h>


// --- PARCHE PARA SERIALIZAR PAIRS (std::pair) ---
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

		//funcion auxiliar para empaquetar el ReadBatches

		void pack_readbatch( const ReadBatches& rb, char* buffer, int bufferSize, int& position, MPI_Comm comm)
		{
			int paired = rb.paired ? 1 : 0; //paired es un bool y lo vamos a enviar como un int
			int number_reads = rb.ids.size();

			MPI_Pack( &paired, 1, MPI_INT, buffer, bufferSize, &position, comm ); //empaquetar el valor de paired 
			MPI_Pack( &number_reads, 1, MPI_INT, buffer, bufferSize, &position, comm ); //empaquetar el numero de lecturas
												    //empaquetar los ids de las lecturas que son strings 
			for (const std::string& id : rb.ids)
			{
				int idlen = id.size();
				MPI_Pack( &idlen, 1, MPI_INT, buffer, bufferSize, &position, comm ); //empaquetar la longitud del id
				MPI_Pack( id.c_str(), idlen, MPI_CHAR, buffer, bufferSize, &position, comm ); //empaquetar el id, c_str() devuelve un puntero un const char * con el contenido del 
													      //string pero terminado en '\0' y MPI_Pack lo empaqueta como un array de caracteres
													      //std:string no es reconocido por MPI directamente y para empaquetarlo necesito un 
													      //const char * con el contenido del string, el idlen no enviara el \0 para ello 
													      //tendria que ser idlen +1 ya que size no cuenta el \0

			}
			//Espaquetar las secuencias de los nucleotidos 
			for (const auto& seq : rb.seqs)
			{
				int seqlen = seq.size();
				MPI_Pack( &seqlen, 1, MPI_INT, buffer, bufferSize, &position, comm ); //empaquetar la longitud de la secuencia
				MPI_Pack( seq.data(), seqlen, MPI_CHAR, buffer, bufferSize, &position, comm ); //empaquetar la secuencia de nucleotidos
			}

			//empaquetar las secuencias de los nucleotidos del segundo par si existe (modo paired)
			if(paired){
				for (const auto& seq2 : rb.seqs2){
					int seqlen2 = seq2.size();
					MPI_Pack( &seqlen2, 1, MPI_INT, buffer, bufferSize, &position, comm ); //empaquetar la longitud de la secuencia del segundo para
					MPI_Pack( seq2.data(), seqlen2, MPI_CHAR, buffer, bufferSize, &position, comm ); //empaquetar la secuencia del segundo par
				}
			}
		}

		//funcion auxiliar para desempaquetar el ReadBatches 

		void unpack_readbatch(ReadBatches& rb, const char* buffer, int bufferSize, int& position, MPI_Comm comm)
		{
			int paired = 0;
			int number_reads = 0;

			MPI_Unpack(buffer, bufferSize, &position, &paired, 1, MPI_INT, comm ); //desempaquetar el valor de paired 
			MPI_Unpack(buffer, bufferSize, &position, &number_reads, 1, MPI_INT, comm ); //desempaquetar el numero de lecturas
												     //
			rb.paired = (paired == 1); //convierto el int de nuevo a bool true o false (1 o 0)
			rb.ids.resize(number_reads); //redimensionar el vector de ids al numero de lecturas 
			rb.seqs.resize(number_reads); //redimensionar el vector de secuencias al numero de lecturas 
			if(rb.paired)
				rb.seqs2.resize(number_reads); //redimensionar el vector de secuencias del segundo par al numero de lecturas 

			//toca desempaquetar los ids 
			for (int i = 0; i < number_reads; ++i){
				int idlen = 0;
				MPI_Unpack(buffer, bufferSize, &position, &idlen, 1, MPI_INT, comm ); //desempaquetar la longitud del id 
				std::vector<char> temp_id(idlen); //vector temporal para almacenar el id
				MPI_Unpack(buffer, bufferSize, &position, temp_id.data(), idlen, MPI_CHAR, comm ); //desempaquetar el id 
				rb.ids[i] = std::string(temp_id.data(), idlen); //convertir el vector de char a string y asignarlo al vector de ids 
			}

			//toca desempaquetar las secuencias de nucleotidos (seqs)
			for (int i = 0; i < number_reads; ++i){
				int seqlen = 0;
				MPI_Unpack(buffer, bufferSize, &position, &seqlen, 1, MPI_INT, comm ); //desempaquetar la longitud de la secuencia 
				std::vector<seqan3::dna4> temp_seq(seqlen); //vector temporal para almacenar la secuencia 
				MPI_Unpack(buffer, bufferSize, &position, temp_seq.data(), seqlen, MPI_CHAR, comm ); //desempaquetar la secuencia 
				rb.seqs[i] = std::move(temp_seq); //mover la secuencia al vector de secuencias
			}

			//toca desempaquetar las secuencias del segundo par (seqs2) si existe (modo paired) 
			if(rb.paired){
				for (int i = 0; i < number_reads; ++i){
					int seqlen2 = 0;
					MPI_Unpack(buffer, bufferSize, &position, &seqlen2, 1, MPI_INT, comm ); //desempaquetar la longitud de la secuencia del segundo par
					std::vector<seqan3::dna4> temp_seq2(seqlen2); //vector temporal para almacenar la secuencia del segundo par
					MPI_Unpack(buffer, bufferSize, &position, temp_seq2.data(), seqlen2, MPI_CHAR, comm ); //desempaquetar la secuencia del segundo par
					rb.seqs2[i] = std::move(temp_seq2); //mover la secuencia del segundo par al vector de secuencias del segundo par
				}
			}
		}



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

			// funcion serialize para decirle a Cereal como tiene que guardar el struct
			template <class Archive>
				void serialize( Archive& ar ){
					ar (matches, lca_reads, unique_reads);
				}
		};

		typedef robin_hood::unordered_map< std::string, detail::Rep >  TRep;
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
		bool ganon_classify( Config config, int rank, int numeroProcesos )
		{
			// Start time count
			StopClock timeGanon;
			timeGanon.start();

			auto parsed_hierarchy = detail::parse_hierarchy( config );

			if ( config.verbose && rank == 0 ) //que solo la imprima el rank 0, no todos
				detail::print_hierarchy( config, parsed_hierarchy );

			// Initialize variables
			StopClock timeLoadFilters;
			StopClock timeClassPrint;

			//solo el proceso 0 va a escribir el reporte final solo el debe abrirlo


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


			//Va a ser el proceso 0 el que lea las lecturas y las reparta a todos los demas procesos inlcuyendose a si mismo para clasificiarlas

			// Queues for internal read handling
			// queue1 get reads from file
			// queue2 will get unclassified reads if hierachy == 2
			// if hierachy == 3 queue1 is used for unclassified and so on
			// config.n_batches*config.n_reads = max. amount of reads in memory

			// cola local para el rank 0 
			SafeQueue< detail::ReadBatches > queue_local( config.n_batches ); //cola local para el proceso 0 donde almacenara las lecturas que le tocan a el para clasificar

			SafeQueue< detail::ReadBatches >  queue1( config.n_batches );
			SafeQueue< detail::ReadBatches >  queue2;
			//SafeQueue< detail::ReadBatches >* pointer_current = &queue1; // pointer to the queues
			//el puntero a la cola dependera del rank 
			SafeQueue< detail::ReadBatches >* pointer_current = ( rank == 0) ? &queue_local : &queue1; // pointer to the queues, el proceso 0 usara la cola local y los demas procesos usaran la cola 1 para recibir las lecturas del proceso 0
			SafeQueue< detail::ReadBatches >* pointer_helper  = &queue2; // pointer to the queues
			SafeQueue< detail::ReadBatches >* pointer_extra;             // pointer to the queues

			// Define one threads for decompress bgzf files
			seqan3::contrib::bgzf_thread_count = 1u;

			/*
			// Thread for reading input files
			std::future< void > read_task = std::async(
			std::launch::async, detail::parse_reads, std::ref( queue1 ), std::ref( stats ), std::ref( config ) );

*/
			std::future<void> read_task; //tarea asincrona para leer las lecturas, solo el proceso 0 la va a ejecutar
			std::future<void> distribute_task; //tarea asincrona para distribuir las lecturas, solo la ejecuta el proceso 0, para que no haya deadlock con el hilo principal si hay muchas lecturas 

			//el proceso 0 va a leer los archivos de entrada
			if ( rank == 0 )
			{
				//detail::parse_reads( queue1, stats, config );
				//lanzar la tarea asincrona para leer las lecturas en un hilo separado y paralelo para evitar que el proceso 0 se quede esperando a que se lean todas las lecturas y evitar un deadlock por si no hay suficiente memoria nunca se vaciaria la cola (queue1) y los demas procesos se quedarian esperando a que el proceso 0 les envie las lecturas pero este no las enviaria al no haber leido todas las lecturas (insuficiente memoria) y al no enviarlas nunca vaciaria la cola (queue1)
				read_task = std::async(
						std::launch::async, detail::parse_reads, std::ref( queue1 ), std::ref( stats ), std::ref( config ) );
			}

			//para repartir las lecturas a los demas procesos usarios Send y Rec pero estas funciones no valen para pasar estructuras de datos complejas que puedan tener punteros
			//a direcciones de memoria donde la memoria no esta compartida entre los diferentes procesos y solo pasarias el puntero y no la info que contiene. 
			//La solucion es usar las funciones de OpenMPI; MPI_Pack y MPI_Unpack que empaquetan y desempaquetan los datos en un buffer de memoria compartida entre todos los procesos,
			//usare dos funciones auxialiares que son pack_readbatch y unpack_readbatch que empaquetan y desempaquetan los datos de ReadBatches.
			//

			//el proceso 0 es el que empaqueta y envia las lecturas a los demas procesos 

			if ( rank == 0 )
			{
				//la explicacion del tamano para el bufferSize viene explicado de que, seqan3:dna4 para expresar las bases del ADN (A,C, G, T) usa 2 bits por base, 
				//(https://docs.seqan.de/seqan3/3.0.1/classseqan3_1_1dna4.html)
				//entonces, en 1 byte ( 8 bits ) me caben 4 bases. Entonces el vector de seqan3:dna4 , cada base son 2 bits, el tamano de lote se configura en n_reads que por defecto son 400
				//el tamano por defecto de 400 del n_reads viene en la linea 45 de (desde ganon_classify) en /include/ganon_classify/Config.hpp [ linea 45 ]
				// No sabemos cuantas bases va a tener cada lectura, pero si sabemos que el tamano de lote es de 400 lecturas,
				//Los IDS son std::string , No sabemos cuanto va a ser la cabecera de promedio por ID, 400 x X  bytes, 
				//si es el proceso 0 , empaquetamos las lecturas y las enviamos
				//Si se establece un tamano de buffer arbitrario puede pasar que en una lectura determinada, el buffer no llegue, aunque esto no ocurriera, estariamos desperdiciando memoria.
				// La solucion es usar un buffer de tamano variable, que se ajuste al numero de lecturas y a la longitud de las lecturas.

				//PENSAR EN LO DEL TAMANO DEL BUFFER; DE MOMENTO QUEDA CON UN BUFFERSIZE DE 1 MB (deberia llegar)

			distribute_task = std::async( std::launch::async, [&queue1, &queue_local, numeroProcesos, rank](){
				int proceso = 0; //variable para controlar a que proceso se envian las lecturas
				//en vez de queue1.empty a veces si el lector aun no metio nada en la cola, salta la condiion y acaba inmediatamente el bucle, para ello lo he cambiado a while true y dentro que rompa con una condicion si ya no hay ninguna lectura mas, uso queue.pop que tiene en cuenta si el lector ha terminado o no
				while ( true )
				{
					//sacamos un lote de lecturas de la cola 
					detail::ReadBatches rb = queue1.pop();
					//si el lote de lecturas esta vacio (seqs, seqs2 e ids) es que el lector ha terminado
					if ( rb.ids.empty() && rb.seqs.empty() && rb.seqs2.empty() )
					{
						break;
					}

					// si el proceso destino es el proceso 0 no hay que enviar lecturas que se las queda el rank 0 oara clasificiarlas
					if ( proceso == 0 ){
						//el proceso 0 se queda con la lectura para clasificarla, la mete en su cola local y no la envian por tanto
						queue_local.push( std::move( rb ) ); 
					}
					else {
						//toca empaquetar 
						int pack_size = 0; //variable para controlar el tamano del buffer empaquetado 
						int temp_size; //variable auxiliar temporal para MPI PackSize 

						//si es paired (es un int)
						MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &temp_size ); //tamano necesario para empaquetar el int de paired
						pack_size += temp_size;
						//numero de reads (es un int )
						MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &temp_size ); //tamano necesario para empaquetar el int de numero de reads
						pack_size += temp_size;
						//ids (vector de strings)
						for (const auto& id : rb.ids){
							int id_size = id.size();
							MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &temp_size ); //tamano necesario para empaquetar el int de tamano de id
							pack_size += temp_size;
							MPI_Pack_size( id_size, MPI_CHAR, MPI_COMM_WORLD, &temp_size ); //tamano necesario para empaquetar el string del id
							pack_size += temp_size;
						}
						//seqs (Read1) (vector de seqan3:dna4)
						for (const auto& seq : rb.seqs){
							int seq_size = seq.size();
							MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &temp_size ); //tamano necesario para empaquetar el int de tamano de seq
							pack_size += temp_size;
							MPI_Pack_size( seq_size, MPI_CHAR, MPI_COMM_WORLD, &temp_size ); //tamano necesario para empaquetar el vector de seqan3:dna4
							pack_size += temp_size;
						}
						//si es paired, seqs2 (Read2) (vector de seqan3:dna4)
						if ( rb.paired ){
							for (const auto& seq2 : rb.seqs2){
								int seq2_size = seq2.size();
								MPI_Pack_size( 1, MPI_INT, MPI_COMM_WORLD, &temp_size ); //tamano necesario para empaquetar el int de tamano de seq2
								pack_size += temp_size;
								MPI_Pack_size( seq2_size, MPI_CHAR, MPI_COMM_WORLD, &temp_size ); //tamano necesario para empaquetar el vector de seqan3:dna4_traits
								pack_size += temp_size;
							}
						}

						//creo un buffer dinmaico y empaqueto 

						std::vector<char> dynamic_buffer(pack_size); // buffer dinamico con el tamano necesario
						int pos = 0;
						pack_readbatch( rb, dynamic_buffer.data(), pack_size, pos, MPI_COMM_WORLD ); // empaquetamos el lote de lecturas en el buffer 

						//enviamos el buffer al proceso que toque en esta iteraccion del bucle solo si hay mas de un proceso claro 
						if ( numeroProcesos > 1 ){
							MPI_Send( &pos, 1, MPI_INT, proceso, 0, MPI_COMM_WORLD ); // enviamos pos para que receptor sepa cuantos bytes tiene que leer 
							MPI_Send( dynamic_buffer.data(), pos, MPI_PACKED, proceso, 0, MPI_COMM_WORLD ); // enviamos el buffer empaquetado al proceso que toque (pack_size bytes)
						}														  //incrementamos el proceso al que enviamos las lecturas con round robin 
					}

					proceso++; //incrementamos el proceso al que enviamos las lecturas
					if ( proceso >= numeroProcesos ) proceso = 0; //si el proceso es mayor o igual al numero de procesos, volvemos al proceso 0 (round robin)
				}
				//queue2 es una cola pa almacenar las lecturas no clasificadas si hay mas de una jerarquia, se gestiona en cada memoria local de cada proceso
				/*
				//si hay paired reads en queue2
				int proceso = 1; //empezamos a enviar al proceso 1
				while ( !queue2.empty() )
				{
				//sacamos un lote de lecturas de la cola 
				detail::ReadBatches rb = queue2.pop();
				//si el lote de lecturas esta vacio (seqs, seqs2 e ids) hacemos continue para sacar otro lote de lecturas en la siguiente iteracion del bucle
				if ( rb.ids.empty() && rb.seqs.empty() && rb.seqs2.empty() )
				{
				continue;
				}

				//toca empaquetar 
				int pos = 0;
				pack_readbatch( rb, buffer, BUFFERSIZE, pos, MPI_COMM_WORLD ); // empaquetamos el lote de lecturas en el buffer 
											       //enviamos el buffer al proceso que toque en esta iteraccion del bucle 
											       MPI_Send( &pos, 1, MPI_INT, proceso, 1, MPI_COMM_WORLD ); // enviamos pos para que receptor sepa cuantos bytes tiene que leer , le pongo 1 (etiqueta distintia al queue1)
											       MPI_Send( buffer, pos, MPI_PACKED, proceso, 1, MPI_COMM_WORLD ); // enviamos el buffer empaquetado al proceso que toque (pos bytes)
																				//incrementamos el proceso al que enviamos las lecturas con round robin 
																				proceso = (proceso % (numeroProcesos - 1)) + 1; // el proceso 0 no envia lecturas a si mismo, empieza a enviar al proceso 1 y luego al 2, etc. ( del 1 al numeroProcesos -1 y una vez llega pues otra vez al 1 y asi)
																				}*/
				// Enviar senal a todos los demas procesos para indicar que no hay mas lecturas
				//
				//senal fin lectura queue1 
				int fin = -1; //-1 indicara que no hay mas lecturas 
				for ( int i = 1; i < numeroProcesos; ++i )
				{
					MPI_Send( &fin, 1, MPI_INT, i, 0, MPI_COMM_WORLD ); // enviamos fin de lectura queue1
				}
				/*
				//senal fin lectura queue2
				for ( int i = 1; i < numeroProcesos; ++i )
				{
				MPI_Send( &fin, 1, MPI_INT, i, 1, MPI_COMM_WORLD ); // enviamos fin de lectura queue2 (la etiqueta pa queue1 era 0 y para queue2 es 1)
				}*/
				//esperar a que el hilo de lectura termine 
				//read_task.get();
				//notificar a la cola local del proceso 0 que no hay mas lecturas para clasificarlas
				queue_local.notify_push_over(); 
			});} //fin del hilo paralelo de distribucion de lecturas 


			//declarar receive task 
			std::future< void > receive_task;

			// toca recibir las lecturas de los demas procesos para proceder a clasificarlas 
			if ( rank != 0 )
			{
				//al recibir no vamos a recibir y procesar si no que vamos a recibir y guardar en una cola ya que si no bloquearia a MPI_REcv
				//recibir lecturas de queue1 
				//lanzar la recepcion de lecturas en un hilo paralelo, asi el proceso puede empezar a clasificar lecturas mientras sigue recibiendo mas lecturas 
				receive_task = std::async( std::launch::async, [&queue1]() {

						while ( true ){
						int pos;
						MPI_Recv( &pos, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); // recibimos el numero de bytes que tenemos que leer 
						if ( pos == -1 ) // si pos es -1, significa que no hay mas lecturas
						{
						break; // salimos del bucle
						}
						//buffer dinamico para recibir
						std::vector<char> buffer(pos); // creamos un buffer de tamano pos que sera donde recibimos los datos 
						MPI_Recv( buffer.data(), pos, MPI_PACKED, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); // recibimos el buffer empaquetado 
						int offset = 0; // desempaquetamos el buffer en un ReadBatches 
						detail::ReadBatches rb;
						unpack_readbatch( rb, buffer.data(), pos, offset, MPI_COMM_WORLD ); // desempaquetamos el buffer en un ReadBatches
														    //metemos el ReadBatches en la cola de lectura si no esta vacio 
						if ( !rb.ids.empty() || !rb.seqs.empty() || !rb.seqs2.empty() ) // si el ReadBatches no esta vacio
						{
						queue1.push( std::move( rb ) ); // metemos el ReadBatches en la cola de lectura 
						}
						}
						/*
						//recibir lecturas de queue2 ( si es que hay paired reads ) 
						while ( true )
						{
						int pos;
						MPI_Recv( &pos, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); // recibimos el numero de bytes que tenemos que leer 
						if ( pos == -1 ) // si pos es -1, significa que no hay mas lecturas
						{
						break; // salimos del bucle
						}
						char buffer[pos]; // creamos un buffer de tamano pos que sera donde recibimos los datos 
						MPI_Recv( buffer, pos, MPI_PACKED, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE ); // recibimos el buffer empaquetado 
						int offset = 0; // desempaquetamos el buffer en un ReadBatches 
						detail::ReadBatches rb;
						unpack_readbatch( rb, buffer, pos, offset, MPI_COMM_WORLD ); // desempaquetamos el buffer en un ReadBatches
													     //metemos el ReadBatches en la cola de lectura si no esta vacio 
													     if ( !rb.ids.empty() || !rb.seqs.empty() || !rb.seqs2.empty() ) // si el ReadBatches no esta vacio
													     {
													     pointer_current->push( std::move( rb ) ); // metemos el ReadBatches en la cola de lectura 
													     }
													     }*/ 
						//notificar que no hay mas lecturas 
						queue1.notify_push_over();
				}); //fin del hilo paralelo de recepcion de lecturas
			}



			// Thread for printing unclassified reads
			SafeQueue< detail::ReadOut > unclassified_queue;
			std::future< void >          write_unclassified_task;
			//cada rank que escriba en su propio archivo de unclassified
			std::string out_unclassified_file = config.output_prefix + ".unc.rank";
			//si el rank no es el 0 (que escribe en el .unc) le añado el numero de rank al nombre del archivo de unclassified 
			if (rank > 0)
				out_unclassified_file += ".rank" + std::to_string(rank);
			if ( config.output_unclassified && !config.output_prefix.empty() )
			{
				write_unclassified_task = std::async( std::launch::async,
						detail::write_unclassified,
						std::ref( unclassified_queue ),
						out_unclassified_file );
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
					// nombre de archivo de salida distinto para cada rank  
					std::string hierarchy_output_file_lca = hierarchy_config.output_file_lca;
					std::string hierarchy_output_file_all = hierarchy_config.output_file_all;	
					if (rank > 0){
						if( hierarchy_config.output_file_lca.empty() == false ){
							hierarchy_output_file_lca += ".rank" + std::to_string(rank);
						}
						if( hierarchy_config.output_file_all.empty() == false ){
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
								hierarchy_last ) );
				}

				// Wait here until classification is over
				for ( auto&& task : tasks )
				{
					task.get();
				}

				// After classification, no more reads are going to be pushed to the output
				classified_all_queue.notify_push_over();
				classified_lca_queue.notify_push_over();

				// Sum reports of each threads into one
				detail::TRep rep = sum_reports( reports );

				//el mapa rep contiene las bacterias que han hecho match

				// Sum totals for each thread and report into stats
				stats.add_totals( hierarchy_label, totals );
				stats.add_reports( hierarchy_label, rep );

				//esperar por el hilo de distribucion, lectura y recibo 
				if ( rank == 0 && hierarchy_first ){
					if ( read_task.valid() ) read_task.get();
					if ( distribute_task.valid() ) distribute_task.get();
				}
				if ( rank != 0 && hierarchy_first ){
					if ( receive_task.valid() ) receive_task.get();
					//esperar a que el hilo de recepcion de lecturas termine para asegurarnos que ya no quedan mas lecturas por recibir antes de proceder a clasificar las lecturas recibidas y evitar un posible deadlock
				}

				//una vez estan todos los reports y stats, las mandamos al rank 0 para que tenga un report y stats global.
				size_t local_classified = stats.total.reads_classified;
				size_t local_unique = stats.total.unique_matches;
				size_t local_matches = stats.total.matches;
				size_t local_processed = stats.total.reads_processed;

				size_t global_classified;
				size_t global_unique;
				size_t global_matches;
				size_t global_processed;

				//ahora toca reducir todos los totales de los stats y rep al global en el rank 0 
				MPI_Reduce( &local_classified, &global_classified, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
				MPI_Reduce( &local_unique, &global_unique, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
				MPI_Reduce( &local_matches, &global_matches, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
				MPI_Reduce( &local_processed, &global_processed, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

				//una vez hecho el reduce
				if ( rank == 0 ){
					stats.total.reads_classified = global_classified;
					stats.total.reads_processed = global_processed;
					//el total reads_classified y el reads_processed el mapa no lo sabe, el unique matches y total matches si, pongo a 0 los del rank 0
					//stats.total.unique_matches = global_unique;
					//stats.total.matches = global_matches;
					stats.total.unique_matches = 0;
					stats.total.matches = 0;
					//limpiar los contadores por jerarquia
					for (auto& [etiqueta,j] : stats.hierarchy_total){
						j.matches = 0;
						j.unique_matches = 0;
					}
				}


				//unificar el mapa rep usando cereal para serializar y deserializar el mapa y enviarlo al rank 0. Podria intentar usar Mpi pack y unpack y un objetivo creado para el mapa pero el mapa es una coleccion de punteros distribuidos en memoria y por tanto no estan en memoria con un desplazamiento fijo y se vuelve todo muy complejo y complicado.


				if (rank != 0){
					//si no soy el rank 0
					std::stringstream ss;
					{
						cereal::BinaryOutputArchive oarchive( ss );

						//como cereal no sabe trabajar con robin_hood map directamente, lo convierto temporalmente en un vector
						//convierto el mapa a un vector con un par clave valor para que cereal lo entienda
						std::vector<std::pair<std::string, detail::Rep>> temp_vec;
						temp_vec.reserve(rep.size());
						for (auto const& [clave, valor] : rep) {
							temp_vec.push_back({clave,valor});
						}

						//oarchive( rep);
						oarchive(temp_vec); //serializa el vector
					}
					std::string serialized_str = ss.str();
					int len = serialized_str.size();
					//envio el tamano del string serializado
					MPI_Send( &len, 1, MPI_INT, 0, 0, MPI_COMM_WORLD );
					//envio el string serializado
					MPI_Send( serialized_str.c_str(), len, MPI_CHAR, 0, 0, MPI_COMM_WORLD );
				}else{
					//si soy el rank 0
					for ( int i = 1; i < numeroProcesos; ++i ){
						//recibo el tamano del string serializado
						int len;
						MPI_Recv( &len, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
						//recibo el string serializado
						std::vector<char> buffer(len);
						MPI_Recv( buffer.data(), len, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
						std::string serialized_str( buffer.data(), len );
						//deserializo el string recibido y lo uno al rep global
						std::stringstream ss_in( serialized_str );

						//vector temporal para recibir los datos
						std::vector<std::pair<std::string, detail::Rep>> temp_vec;

						{
							cereal::BinaryInputArchive iarchive(ss_in);
							//deserializar el vector
							iarchive(temp_vec);
						}

						//sumo temp_rep a rep (fusionar los dos mapas)
						for (auto const& [target, r] : temp_vec){
							if(target.empty()){
								continue;
							}
							rep[target].matches += r.matches;
							rep[target].lca_reads += r.lca_reads;
							rep[target].unique_reads += r.unique_reads;
						}
					}
				}

				//ahora rep en el rank 0 tiene el reporte global de todos los demas procesos, añadimos los totales globales 
				//solo el rank 0 escribe el reporte final 
				if( rank == 0){
					stats.add_reports( hierarchy_label, rep ); //añadimos el reporte global a los stats globales
										   // write reports, escribe el reporte final solo el rank 0 que ya tiene todos los reports de todos los demas procesos
					detail::write_report( rep, tax, out_rep, hierarchy_label );
				}

				//todos los procesos esperan a que se terminen de escribir sus hilos de escritura de archivos .lca y .all
				// Wait here until all files are written
				for ( auto&& task : write_tasks )
				{
					task.get();
				}
				timeClassPrint.stop();


				// Close file for writing (if not STDOUT), cada rank cierra sus propios archivos de salida .lca y .all
				if ( !config.output_prefix.empty() )
				{
					if ( config.output_lca  && out_lca.is_open() )
						out_lca.close();
					if ( config.output_all && out_all.is_open() )
						out_all.close();
				}

				if(rank == 0){
					//sincronizar todos los procesos antes de continuar a la siguiente jerarquia
					if ( hierarchy_first )
					{
						//read_task.get();                    // get reading tasks at the end of the first hierarchy
						pointer_helper->notify_push_over(); // notify push is over, only on first time (will be always set over for
										    // next iterations)
					}}
			}//fin del bucle de jerarquias

			// Wait here until all unclassified reads are written
			if ( config.output_unclassified )
			{
				unclassified_queue.notify_push_over();
				write_unclassified_task.get();
			}

			//solo el rank 0 que imprima las estadisticas finales totales 
			if( rank == 0){
				out_rep << "#total_classified\t" << stats.total.reads_classified << '\n';
				// account for unclassified and skipped sequences
				out_rep << "#total_unclassified\t" << stats.input_reads - stats.total.reads_classified << '\n';
				if ( !config.output_prefix.empty() )
				{
					out_rep.close();
				}
			}

			timeGanon.stop();

			//solo el rank 0 imprime las estadisticas 
			if( rank == 0 ){
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
			return ganon_classify< detail::THIBF >( config , rank, numeroProcesos );
		else
			return ganon_classify< detail::TIBF >( config, rank, numeroProcesos );
	}

} // namespace GanonClassify
