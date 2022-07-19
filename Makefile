#modied from htslib makefile
FLAGS2=-O3 -std=c++11 

CFLAGS += $(FLAGS2)
CXXFLAGS += $(FLAGS2)

CSRC = $(wildcard *.c) 
CXXSRC = $(wildcard *.cpp)
OBJ = $(CSRC:.c=.o) $(CXXSRC:.cpp=.o)



LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl

CRYPTO_TRY=$(shell echo 'int main(){}'|g++ -x c++ - -lcrypto 2>/dev/null; echo $$?)
ifeq "$(CRYPTO_TRY)" "0"
$(info Crypto library is available to link; adding -lcrypto to LIBS)
LIBS += -lcrypto
else
$(info Crypto library is not available to link; will not use -lcrypto)
endif



all: ngsLCA

PACKAGE_VERSION  = 0.9

ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define NGSLCA_VERSION "$(PACKAGE_VERSION)"' > $@




# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
$(info HTSSRC defined)
HTS_INCDIR=$(realpath $(HTSSRC))
HTS_LIBDIR=$(realpath $(HTSSRC))/libhts.a
else
$(info HTSSRC not defined, assuming systemwide installation -lhts)
LIBS += -lhts
endif

-include $(OBJ:.o=.d)

ifdef LDFLAGS
FLAGS += $(LDFLAGS)
else
FLAGS += $(FLAGS2)
endif

ifdef HTSSRC
%.o: %.c
	$(CC) -c $(CFLAGS) -I$(HTS_INCDIR) $*.c
	$(CC) -MM $(CFLAGS) -I$(HTS_INCDIR) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -I$(HTS_INCDIR) $*.cpp
	$(CXX) -MM $(CXXFLAGS) -I$(HTS_INCDIR) $*.cpp >$*.d

ngsLCA: version.h $(OBJ)
	$(CXX) $(FLAGS) -o ngsLCA *.o $(HTS_LIBDIR) $(LIBS)
else
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c
	$(CC) -MM $(CFLAGS) $*.c >$*.d

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d

ngsLCA: version.h $(OBJ)
	$(CXX) $(FLAGS) -o ngsLCA *.o $(LIBS)
endif

clean:	
	rm -f ngsLCA *.o *.d version.h

force:

