#modied from htslib makefile
FLAGS2=-O3 -std=c++11 

CFLAGS += $(FLAGS2)
CXXFLAGS += $(FLAGS2)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)

all: ngsLCA


# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
endif


-include $(OBJ:.o=.d)

ifdef LDFLAGS
FLAGS += $(LDFLAGS)
else
FLAGS += $(FLAGS2)
endif

ifdef HTSSRC
%.o: %.c
	$(CC) -c  $(CFLAGS) -I$(HTS_INCDIR) $*.c
	$(CC) -MM $(CFLAGS)  -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS)  -I$(HTS_INCDIR) $*.cpp >$*.d

ngsLCA: $(OBJ)
	$(CXX) $(FLAGS)  -o ngsLCA *.o $(HTS_LIBDIR) -lz -llzma -lbz2 -lpthread -lcurl 
else
%.o: %.c
	$(CC) -c  $(CFLAGS)  $*.c
	$(CC) -MM $(CFLAGS)  $*.c >$*.d

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS)  $*.cpp
	$(CXX) -MM $(CXXFLAGS)  $*.cpp >$*.d

ngsLCA: $(OBJ)
	$(CXX) $(FLAGS)  -o ngsLCA *.o -lz -llzma -lbz2 -lpthread -lcurl -lhts
endif

clean:	
	rm  -f ngsLCA *.o *.d

