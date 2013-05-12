CPP=g++
CFLAGS=-c -Wall -Wno-deprecated -march=amdfam10 -O2
LDLIBS=-L./triangulate-1.4 -ltriangulate -L/usr/lib -lz -lvte -lvtkCommon -lvtkGraphics \
	-lvtkIO -lvtkFiltering -lvtkRendering -lvtkImaging
INCLUDES=-I/usr/include/vtk-5.10
OPENMP=-fopenmp

OBJECTS=main.o functions.o globals.o interp.o border.o \
		debug.o output.o

main : $(OBJECTS)
	$(CPP) -o main $(OBJECTS) $(LDLIBS) $(OPENMP)

main.o : main.cpp globals.cpp main.h
	$(CPP) $(CFLAGS) $(INCLUDES) $(OPENMP) main.cpp

functions.o : functions.cpp main.h
	$(CPP) $(CFLAGS) functions.cpp

globals.o : globals.cpp main.h
	$(CPP) $(CFLAGS) globals.cpp

interp.o : interp.c interp.h
	$(CPP) $(CFLAGS) interp.c

debug.o : debug.c debug.h
	$(CPP) $(CFLAGS) debug.c

border.o : border.c border.h
	$(CPP) $(CFLAGS) border.c

output.o : output.c output.h
	$(CPP) $(CFLAGS) $(INCLUDES) output.c

clean:
	rm -rf *o main
