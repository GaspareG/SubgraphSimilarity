graph-gen: graph_generator

induced-naive: k-induced-path-naive

induced-color: k-induced-path-color-coding

path-naive: k-path-naive

path-color: k-path-color-coding

path-divide: k-path-divide-color

all: graph-gen induced-naive induced-color path-naive path-color path-divide

clean:
	rm graph_generator
	rm k-induced-path-naive
	rm k-induced-path-color-coding
	rm k-path-naive
	rm k-path-color-coding
	rm k-path-divide-color
