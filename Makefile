ifdef USE_INT
MACRO = -DUSE_INT
endif

#compiler setup
CXX = g++
MPICXX = mpic++
CXXFLAGS = -std=c++14 -O3 -pthread $(MACRO)

COMMON= core/utils.h core/cxxopts.h core/get_time.h core/graph.h core/quick_sort.h
SERIAL= apsp_serial
THREADS= apsp_threads
MPI= apsp_mpi
ALL= $(SERIAL) $(THREADS)

all : $(ALL)

$(SERIAL): %: %.cpp $(COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(THREADS): %: %.cpp $(COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(MPI): %: %.cpp $(COMMON)
	$(MPICXX) $(CXXFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)
