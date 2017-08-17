
all: 
		g++ -ggdb -Wall -Wextra -march=native -O3 -std=c++11 Bucket.cpp Seed.cpp Fswm.cpp pattern.cpp Sequence.cpp Word.cpp -o fswm -fopenmp
 
