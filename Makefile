CXXFLAGS += --std=c++11 -Wall -pedantic -O2 -DMAKE_VALGRIND_HAPPY -fopenmp
objects = graph_generator k-path-color-coding k-path-color-coding-parallel k-induced-path-color-coding k-path-divide-color slash-burn k-induced-path-naive k-path-naive k-labeled-dpc k-labeled-dpl k-labeled-dpl-tau k-labeled-dplw k-path-jaccard-naive k-path-color-coding-jaccard final final_benchmark1 final_benchmark2 final-freq final_node

$(objects): %: %.cpp
#	$(CC) $(CXXFLAGS) -o $@ $<

all: $(objects)

clean-bin:
	rm $(objects)

clean-dataset:
	-rm input/*
	-rm input/snap/*
	-rm input/gen/*
	-rm input/label/*
	-rm input/jaccard/*
	-rmdir input/gen
	-rmdir input/label
	-rmdir input/snap
	-rmdir input/jaccard
	-rmdir input

clean-output:
	-rm output/*
	-rmdir output

clean: clean-dataset clean-bin clean-output

dataset-snap: graph_generator
	mkdir -p input | true
	mkdir -p input/snap | true
	wget -P input/snap https://snap.stanford.edu/data/web-BerkStan.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-Google.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-NotreDame.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-Stanford.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/facebook_combined.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/twitter_combined.txt.gz
	gunzip input/snap/*.gz

dataset-gen: graph_generator
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
	#./graph_generator 10000000  50000000 input/gen/graph-10M-50M.nme.bin
	#./graph_generator 10000000 100000000 input/gen/graph-10M-100M.nme.bin

dataset-label: graph_generator
	mkdir -p input | true
	mkdir -p input/label | true
	./graph_generator     1000      5000 input/label/graph-label-1k-5k.nme.bin     8
	./graph_generator     1000     10000 input/label/graph-label-1k-10k.nme.bin    8
	./graph_generator    10000     50000 input/label/graph-label-10k-50k.nme.bin   8
	./graph_generator    10000    100000 input/label/graph-label-10k-100k.nme.bin  8
	./graph_generator   100000    500000 input/label/graph-label-100k-500k.nme.bin 8
	./graph_generator   100000   1000000 input/label/graph-label-100k-1M.nme.bin   8
	./graph_generator  1000000   5000000 input/label/graph-label-1M-5M.nme.bin     8
	./graph_generator  1000000  10000000 input/label/graph-label-1M-10M.nme.bin    8
	#./graph_generator 10000000  50000000 input/gen/graph-10M-50M.nme.bin
	#./graph_generator 10000000 100000000 input/gen/graph-10M-100M.nme.bin

dataset-jaccard: graph_generator
	mkdir -p input | true
	mkdir -p input/jaccard | true
	./graph_generator     1000      5000 input/jaccard/graph-label-1k-5k-50.nme.bin      4     50     50
	./graph_generator     1000      5000 input/jaccard/graph-label-1k-5k-100.nme.bin     4    100    100
	./graph_generator    10000     50000 input/jaccard/graph-label-10k-50k-500.nme.bin   4    500    500
	./graph_generator    10000     50000 input/jaccard/graph-label-10k-50k-1k.nme.bin    4   1000   1000
	./graph_generator   100000    500000 input/jaccard/graph-label-100k-500k-5k.nme.bin  4   5000   5000
	./graph_generator   100000    500000 input/jaccard/graph-label-100k-500k-10k.nme.bin 4  10000  10000
	./graph_generator  1000000   5000000 input/jaccard/graph-label-1M-5M-50k.nme.bin     4  50000  50000
	./graph_generator  1000000   5000000 input/jaccard/graph-label-1M-5M-100k.nme.bin    4 100000 100000

dataset: dataset-snap dataset-gen dataset-label dataset-jaccard

test-gen: k-path-color-coding-parallel
	./k-path-color-coding-parallel -k 4 -g input/gen/graph-1k-5k.nme.bin     -f nme --verbose
	./k-path-color-coding-parallel -k 4 -g input/gen/graph-10k-50k.nme.bin   -f nme --verbose
	./k-path-color-coding-parallel -k 4 -g input/gen/graph-100k-500k.nme.bin -f nme --verbose
	./k-path-color-coding-parallel -k 4 -g input/gen/graph-1M-5M.nme.bin     -f nme --verbose
	#./k-path-color-coding-parallel -k 6 -g input/gen/graph-10M-50M.nme.bin   -f nme --verbose
	#./k-path-color-coding-parallel -k 6 -g input/gen/graph-10M-100M.nme.bin   -f nme --verbose

test-snap: k-path-color-coding-parallel
	./k-path-color-coding-parallel -k 4 -g input/snap/facebook_combined.txt -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/twitter_combined.txt  -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/web-Google.txt        -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/web-NotreDame.txt     -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/web-Stanford.txt      -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/web-BerkStan.txt      -f snap --verbose
