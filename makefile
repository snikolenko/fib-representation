CC=g++
CFLAGS=-Wall -O3 -fopenmp -std=c++0x
LIB=-lboost_system -lboost_program_options -lboost_filesystem -lboost_date_time -lgomp

UTIL=src/io.hpp src/logging.hpp src/stringutil.hpp

all: bin/fib_onebyone

bin/reduce_fib: src/reduce_fib.cpp src/reduce_fib.hpp $(UTIL) src/rules.hpp
	$(CC) $(CFLAGS) src/reduce_fib.cpp $(LIB) -o bin/reduce_fib

bin/fib_onebyone: src/fib_onebyone.cpp src/reduce_fib.hpp $(UTIL) src/rules.hpp
	$(CC) $(CFLAGS) src/fib_onebyone.cpp $(LIB) -o bin/fib_onebyone

clean:
	rm -rf bin/reduce_fib
