CC = mpic++

# Debugging CFLAGS
CFLAGS = -O3 -march=native -DNDEBUG -funroll-loops -Wall -std=c++23
# Non-debugging CFLAGS
#CFLAGS = -Wall 
INCLUDES = -I/usr/local/include -I. -I./include -I../../libs -I../../libs/cxxopts -I../../libs/robin_hood
LDFLAGS = -O3 -march=native -Wall -lz -lbz2

OBJS = \
  CommandLineParser.o \
  GanonClassify.o \
  main.o 
  
EXE = ganon_classify  

all: ${OBJS}
	${CC} ${OBJS} ${LDFLAGS} -o ${EXE}

CommandLineParser.o: CommandLineParser.cpp include/CommandLineParser.hpp
	${CC} -c ${CFLAGS} ${INCLUDES} CommandLineParser.cpp -o CommandLineParser.o

GanonClassify.o: GanonClassify.cpp include/GanonClassify.hpp include/hierarchical_interleaved_bloom_filter.hpp
	${CC} -c ${CFLAGS} ${INCLUDES} GanonClassify.cpp -o GanonClassify.o
	
main.o: main.cpp 
	${CC} -c ${CFLAGS} ${INCLUDES} main.cpp -o main.o

clean:
	rm -f ${OBJS} ${EXE}  

