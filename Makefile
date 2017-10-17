CXXFLAGS += --std=c++11 -Wall -pedantic -g -DMAKE_VALGRIND_HAPPY -fopenmp
objects = graph_generator k-induced-path-color-coding \
					k-induced-path-naive k-labeled-path-cc-conte \
					k-path-color-coding-parallel k-path-color-coding \
					k-path-divide-color k-path-naive slash-burn \
					k-labeled-path-cc-marino

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
	-rmdir input/gen
	-rmdir input/label
	-rmdir input/snap
	-rmdir input

clean: clean-dataset clean-bin

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

dataset: dataset-snap dataset-gen dataset-label

test-gen: k-path-color-coding-parallel
	./k-path-color-coding-parallel -k 4 -g input/gen/graph-1k-5k.nme.bin     -f nme --verbose
	./k-path-color-coding-parallel -k 4 -g input/gen/graph-10k-50k.nme.bin   -f nme --verbose
	./k-path-color-coding-parallel -k 4 -g input/gen/graph-100k-500k.nme.bin -f nme --verbose
	./k-path-color-coding-parallel -k 4 -g input/gen/graph-1M-5M.nme.bin     -f nme --verbose
	#./k-path-color-coding-parallel -k 6 -g input/gen/graph-10M-50M.nme.bin   -f nme --verbose

test-snap: k-path-color-coding-parallel
	./k-path-color-coding-parallel -k 4 -g input/snap/facebook_combined.txt -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/twitter_combined.txt  -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/web-Google.txt        -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/web-NotreDame.txt     -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/web-Stanford.txt      -f snap --verbose
	./k-path-color-coding-parallel -k 4 -g input/snap/web-BerkStan.txt      -f snap --verbose

test-labeled-conte: k-labeled-path-cc-conte
	./k-labeled-path-cc-conte -k 4 -g input/label/graph-label-1k-5k.nme.bin     --verbose
	./k-labeled-path-cc-conte -k 4 -g input/label/graph-label-10k-50k.nme.bin   --verbose
	./k-labeled-path-cc-conte -k 4 -g input/label/graph-label-100k-500k.nme.bin --verbose
	./k-labeled-path-cc-conte -k 4 -g input/label/graph-label-1M-5M.nme.bin     --verbose

test-labeled-marino: k-labeled-path-cc-marino
	./k-labeled-path-cc-marino -k 4 -g input/label/graph-label-1k-5k.nme.bin     --verbose
	./k-labeled-path-cc-marino -k 4 -g input/label/graph-label-10k-50k.nme.bin   --verbose
	./k-labeled-path-cc-marino -k 4 -g input/label/graph-label-100k-500k.nme.bin --verbose
	./k-labeled-path-cc-marino -k 4 -g input/label/graph-label-1M-5M.nme.bin     --verbose
