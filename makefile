CC=g++
CFLAGS=-Wall -O3 -fopenmp -std=c++0x
LIB=-lboost_system -lboost_program_options -lboost_filesystem -lboost_date_time -lgomp

UTIL=src/io.hpp src/logging.hpp src/stringutil.hpp

all: bin/reduce_fib

bin/reduce_fib: src/reduce_fib.cpp src/reduce_fib.hpp $(UTIL) src/rules.hpp
	$(CC) $(CFLAGS) src/reduce_fib.cpp $(LIB) -o bin/reduce_fib

clean:
	rm -rf bin/reduce_fib
