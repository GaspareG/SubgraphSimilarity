CXXFLAGS += --std=c++17 -g -Wall -pedantic -O0 -DMAKE_VALGRIND_HAPPY -fopenmp
objects = similarity

$(objects): %: %.cpp
#	$(CC) $(CXXFLAGS) -o $@ $<

all: $(objects)

clean-bin:
	rm $(objects)
