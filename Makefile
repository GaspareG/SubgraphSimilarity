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

clean-bin:
	-rm graph_generator
	-rm k-induced-path-naive
	-rm k-induced-path-color-coding
	-rm k-path-naive
	-rm k-path-color-coding
	-rm k-path-color-coding-parallel
	-rm k-path-divide-color
	-rm slash-burn

clean-dataset:
	-rm input/*
	-rm input/snap/*
	-rm input/gen/*
	-rmdir input/gen
	-rmdir input/snap
	-rmdir input

clean: clean-dataset clean-bin

dataset-snap: graph-gen
	mkdir -p input | true
	mkdir -p input/snap | true
	wget -P input/snap https://snap.stanford.edu/data/web-BerkStan.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-Google.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-NotreDame.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-Stanford.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/facebook_combined.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/twitter_combined.txt.gz
	gunzip input/snap/*.gz

dataset-gen: graph-gen
	mkdir -p input | true
	mkdir -p input/gen | true
	./graph_generator     1000      5000 input/gen/graph-1k-5k.nme.bin
	./graph_generator     1000     10000 input/gen/graph-1k-10k.nme.bin
	./graph_generator    10000     50000 input/gen/graph-10k-50k.nme.bin
	./graph_generator    10000    100000 input/gen/graph-10k-100k.nme.bin
	./graph_generator   100000    500000 input/gen/graph-100k-500k.nme.bin
	./graph_generator   100000   1000000 input/gen/graph-100k-1M.nme.bin
	./graph_generator  1000000   5000000 input/gen/graph-1M-5M.nme.bin
	./graph_generator  1000000  10000000 input/gen/graph-1M-10M.nme.bin
	./graph_generator 10000000  50000000 input/gen/graph-10M-50M.nme.bin
	./graph_generator 10000000 100000000 input/gen/graph-10M-100M.nme.bin

dataset: dataset-snap dataset-gen

test-gen: path-color-parallel
	./k-path-color-coding-parallel -k 6 -g input/gen/graph-1k-5k.nme.bin     -f nme --verbose
	./k-path-color-coding-parallel -k 6 -g input/gen/graph-10k-50k.nme.bin   -f nme --verbose
	./k-path-color-coding-parallel -k 6 -g input/gen/graph-100k-500k.nme.bin -f nme --verbose
	./k-path-color-coding-parallel -k 6 -g input/gen/graph-1M-5M.nme.bin     -f nme --verbose
	./k-path-color-coding-parallel -k 6 -g input/gen/graph-10M-50M.nme.bin   -f nme --verbose
