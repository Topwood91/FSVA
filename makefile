all:	index fsva
index:	index.cpp
	g++ index.cpp -o index -O3 -Wall -static
fsva:	fsva.cpp ssw.c ssw_cpp.cpp radix_sort.cpp
	g++ fsva.cpp ssw.c ssw_cpp.cpp -o fsva -O3 -Wall -fopenmp -lz -march=native -static

