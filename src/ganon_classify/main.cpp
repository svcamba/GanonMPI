#include <CommandLineParser.hpp>
#include <GanonClassify.hpp>
#include <mpi.h> //libreria MPI

#include <cstdlib>
#include <utility>

int main( int argc, char** argv )
{
	int procesos; //numero de procesos MPI 
	int idProceso; //id del proceso MPI 

	//inicliazar MPI 
	MPI_Init( &argc, &argv ); 
	MPI_Comm_size( MPI_COMM_WORLD, &procesos ); //obtener el numero de procesos
	MPI_Comm_rank( MPI_COMM_WORLD, &idProceso ); //obtener el id del proceso

	int result = 0;

	if ( auto config = GanonClassify::CommandLineParser::parse( argc, argv ); config.has_value() )
	{
		if ( GanonClassify::run( std::move( config.value() ), procesos, idProceso )) 
			result = EXIT_SUCCESS;
		else
			result = EXIT_FAILURE;
	}
	else
	{
		result = (argc == 1 ? EXIT_FAILURE : EXIT_SUCCESS);
	}
	//esperar a todos los procesos 
	MPI_Barrier( MPI_COMM_WORLD );
	//finalizar MPI 
	MPI_Finalize();
	return result;
}
