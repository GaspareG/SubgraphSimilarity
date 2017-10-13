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
	rm -R input/*.txt | true
	rm -R input/snap/*.txt | true

download:
	mkdir -p input
	mkdir -p input/snap
	wget -P input/snap https://snap.stanford.edu/data/web-BerkStan.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-Google.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-NotreDame.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/web-Stanford.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/facebook_combined.txt.gz
	wget -P input/snap https://snap.stanford.edu/data/twitter_combined.txt.gz
	gunzip input/snap/*.gz
	sed -i '/^[[:blank:]]*#/d;s/#.*//' input/snap/*.txt
	sed -i '1s/^/685231 7600595\n/' input/snap/web-BerkStan.txt
	sed -i '1s/^/875714 5105039\n/' input/snap/web-Google.txt
	sed -i '1s/^/325729 1497134\n/' input/snap/web-NotreDame.txt
	sed -i '1s/^/281903 2312497\n/' input/snap/web-Stanford.txt
	sed -i '1s/^/4040 88234\n/' input/snap/facebook_combined.txt
	sed -i '1s/^/81307 1768149\n/' input/snap/twitter_combined.txt
