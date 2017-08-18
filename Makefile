CXXFLAGS=-ggdb -march=native -O3 -fopenmp
CPPFLAGS=-Wall -Wextra -std=c++11

.PHONY: all clean

all: fswm

fswm: Bucket.o Seed.o Fswm.o pattern.o Sequence.o Word.o
	$(CXX) $(CXXFLAGS) -o $@ $^

clean:
	rm -f fswm *.o
