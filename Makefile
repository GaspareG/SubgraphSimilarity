CXXFLAGS += --std=c++11 -Wall -pedantic -g -DMAKE_VALGRIND_HAPPY -fopenmp

graph-gen: graph_generator

induced-naive: k-induced-path-naive

induced-color: k-induced-path-color-coding

path-naive: k-path-naive

path-color: k-path-color-coding

path-color-parallel: k-path-color-coding-parallel

path-divide: k-path-divide-color

slashburn: slash-burn

all: graph-gen induced-naive induced-color path-naive path-color path-color-parallel path-divide slashburn

clean:
	rm graph_generator
	rm k-induced-path-naive
	rm k-induced-path-color-coding
	rm k-path-naive
	rm k-path-color-coding
	rm k-path-color-coding-parallel
	rm k-path-divide-color
	rm slash-burn

test: all
	k=2
	r=2
	while [ "$r" -ne 124 ]; do
		timeout 5m ./k-path-naive $k < input/graph.1.txt
		retval=$?
		k = $k + 1
	done
