CXXFLAGS += --std=c++17 -Wall -pedantic -O2 -DMAKE_VALGRIND_HAPPY -fopenmp
objects = similarity similarity-author similarity-APC

$(objects): %: %.cpp
#	$(CC) $(CXXFLAGS) -o $@ $<

all: $(objects)

clean-bin:
	rm $(objects)
