#Defines the compiler
CC = g++
#Displays all compilation warnings
CFLAGS += -Wall
#Optimize the code for speed
CFLAGS += -O3 -ffast-math
#Enable debuging with gdb (makes the executable bigger)
#CFLAGS += -ggdb
#Enable gprof based profiling
CFLAGS += -pg
#Allow to use intrinsic functions
CFLAGS += -march=native
#On SSE2 compatible processors this will compile a faster program, in float 
#precision. Comment for a slower program in double precision.
#DEFINES += -DVECTORIZE

all : simulator

simulator : coltree.o atoms.o potential.o constants.o integrator.o common.o \
	main.o
	$(CC) $(CFLAGS) $(DEFINES) $^ -o $@ && mv $@ ../bin/

coltree.o potential.o atoms.o constants.o integrator.o common.o main.o : \
	%.o : %.cpp
	$(CC) $(CFLAGS) $(DEFINES) -c $<

clean :
	rm -rf *.o
