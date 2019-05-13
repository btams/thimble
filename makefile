# name of the library
LIBRARY = thimble

# compiler
CCC = g++

# C++ compiler flags
CCFLAGS = -Wall -Wwrite-strings -ansi -pedantic -O2

# Local stuff for compilation
INCLUDEFLAGS = -I./include/
LIBRARYFLAGS= -L./
LINKFLAGS = -l$(LIBRARY) -lm

# Global directories for installation
INCDIR = /usr/local/include/
LIBDIR = /usr/local/lib/

# library source files.
SRCS = $(wildcard ./src/*.cpp)

# binary source files
BSRC = $(wildcard ./*.cpp)

# library source files
OBJS = $(SRCS:.cpp=.o)

# binary source files
BINS = $(BSRC:.cpp=)

all:	$(OBJS) $(BOBJS)
	ar rcs lib$(LIBRARY).a $(OBJS)
	make binaries
	
.cpp.o:	
	$(CCC) $(INCLUDEFLAGS) $(CCFLAGS) -c $< -o $@

binaries: $(BINS)

.cpp:
	$(CCC) $(CCFLAGS) $(INCLUDEFLAGS) $(LIBRARYFLAGS) $< -o $@ $(LINKFLAGS)

doxy:
	doxygen ./doc/doxy.cfg

install:
	mkdir -p -m 755 $(INCDIR)
	rm -rf $(INCDIR)/$(LIBRARY)
	mkdir -m 755 $(INCDIR)/$(LIBRARY)
	cp -r ./include/$(LIBRARY) $(INCDIR)
	- chmod -R a+r $(INCDIR)/$(LIBRARY)
	mkdir -p -m 755 $(LIBDIR)
	cp -p lib$(LIBRARY).a $(LIBDIR)/lib$(LIBRARY).a
	- chmod a+r $(LIBDIR)/lib$(LIBRARY).a


uninstall:
	rm -f $(LIBDIR)/lib$(LIBRARY).a
	rm -rf $(INCDIR)/$(LIBRARY)

clean:
	rm -f $(OBJS)
	rm -f lib$(LIBRARY).a
	rm -f $(BINS)
	rm -f ./include/*~
	rm -f ./include/thimble/*~
	rm -f ./include/thimble/ecc/*~
	rm -f ./include/thimble/finger/*~
	rm -f ./include/thimble/image/*~
	rm -f ./include/thimble/math/*~
	rm -f ./include/thimble/math/linalg/*~
	rm -f ./include/thimble/math/numbertheory/*~
	rm -f ./include/thimble/math/numerical/*~
	rm -f ./include/thimble/misc/*~
	rm -f ./include/thimble/security/*~
	rm -f ./src/*~
	rm -f ./doc/*~
	rm -f *~

deldoxy:
	rm -f -r ./doc/html
	rm -f -r ./doc/latex	

clobber:
	make clean
	make deldoxy

	

