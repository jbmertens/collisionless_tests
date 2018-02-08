
include Makefile.inc

CXXFLAGS += -std=c++11

LDLIBS += -lfftw3f_omp
LDLIBS += -lfftw3f
LDLIBS += -lz

# Comment this out if you are using the Intel compiler (`icpc`).
LDLIBS += -lm

EXE = collisionless_tests
TEST_EXE = run_tests

all: main.cc *.h
	$(CXX) main.cc -o $(EXE) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

tests: tests.cc *.h
	$(CXX) tests.cc -o $(TEST_EXE) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(EXE) $(TEST_EXE)
