CXXFLAGS += --std=c++11 -Wall -pedantic -O2 -DMAKE_VALGRIND_HAPPY -fopenmp
objects = sim-facebook facebook

$(objects): %: %.cpp
#	$(CC) $(CXXFLAGS) -o $@ $<

all: $(objects)

clean-bin:
	rm $(objects)
