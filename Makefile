CXXFLAGS=-std=c++20 -Wall -Wextra -pedantic
CXXFLAGS_DEBUG=$(CXXFLAGS) -g3 -O0
CXXFLAGS_SANITIZE=$(CXXFLAGS) -O0 -g3 \
				  -fsanitize=address,undefined -fno-omit-frame-pointer \
				  -fno-optimize-sibling-calls
CXXFLAGS_RELEASE=$(CXXFLAGS) -Ofast -march=native -DNDEBUG #-g3
INC=
LIB=
SOURCE=main.cpp
HEADER=params.hpp

main: $(SOURCE) $(HEADER)
	#clang++ $(CXXFLAGS_SANITIZE) -o $@ $(SOURCE) $(INC) $(LIB)
	#clang++ $(CXXFLAGS_DEBUG) -o $@ $(SOURCE) $(INC) $(LIB)
	g++-10 $(CXXFLAGS_RELEASE) -o $@ $(SOURCE) $(INC) $(LIB) #-lprofiler

test: main
	AQTFHE3_VERBOSE=1 ./main

.PHONY: test
