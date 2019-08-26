
#CC = gcc-5
CC = gcc
# CFLAGS = -fopenmp -fno-strict-aliasing -arch x86_64 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes
CFLAGS = -fopenmp -arch x86_64
LFLAGS = -fopenmp -arch x86_64 -bundle -undefined dynamic_lookup
LIBS = -L/Users/daniel/anaconda/lib -pthread
# INCS = -I/Users/daniel/anaconda/lib/python2.7/site-packages/numpy/core/include -I. -I/Users/daniel/anaconda/include/python2.7
INCS = -I/Users/daniel/anaconda/lib/python3.5/site-packages/numpy/core/include -I. -I/Users/daniel/anaconda/include/python3.5

SRC = $(wildcard *.c)
SO_LIBS = $(SRC:.c=.so)

.PHONY: clean

all: $(SO_LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCS) \
		 -c $< -o build/temp.macosx-10.5-x86_64-2.7/$@

%.so: %.o
	$(CC) $(LFLAGS) $(LIBS) \
		 build/temp.macosx-10.5-x86_64-2.7/$< -o $(shell pwd)/$@
clean:
	rm *.so



# gcc -pthread -fno-strict-aliasing -g -O2 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -fPIC -I/home/darothen/miniconda/lib/python2.7/site-packages/numpy/core/include -I. -I/home/darothen/miniconda/include/python2.7 -c hall.c -o build/temp.linux-x86_64-2.7/hall.o


# gcc -pthread -shared -L/home/darothen/miniconda/lib -Wl,-rpath=/home/darothen/miniconda/lib,--no-as-needed build/temp.linux-x86_64-2.7/hall.o -L/home/darothen/miniconda/lib -lpython2.7 -o /home/darothen/workspace/sd-coal/hall.so
