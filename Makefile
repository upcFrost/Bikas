CPP=g++
CFLAGS=-c -Wall -Wno-deprecated -march=native -mtune=native -O3
LDLIBS=-L. -L./lib -L/usr/lib -lvtkCommon -lvtkGraphics \
	-lvtkIO -lvtkFiltering -lvtkRendering -lvtkImaging

INCLUDES=-I. -I./include 

ifeq ($(OS),Windows_NT)
	CFLAGS+=-static
	INCLUDES+=-I"C:\Program Files (x86)\VTK\include\vtk-5.10" 
	LDLIBS+=-L"C:\Program Files (x86)\VTK\lib\vtk-5.10"
	TARGET=main.exe
else
	INCLUDES+=-I/usr/include/vtk-5.10
	LDLIBS+=-L./usr/lib 
	TARGET=main 
endif

OPENMP=-fopenmp

OBJECTS=main.o init.o math_functions.o globals.o interp.o border.o \
		debug.o output.o projectile.o triangulate.o \
		Line2D.o
		
all : $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CPP) -o bin/$(TARGET) $(OBJECTS) $(OPENMP) $(LDLIBS)


main.o : main.cpp globals.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) $(OPENMP) main.cpp

math_functions.o : math_functions.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) math_functions.cpp

globals.o : globals.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) globals.cpp

interp.o : interp.c
	$(CPP) $(CFLAGS) $(INCLUDES) interp.c

debug.o : debug.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) debug.cpp

border.o : border.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) border.cpp

output.o : output.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) $(OPENMP) output.cpp

init.o : init.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) init.cpp
	
projectile.o : projectile.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) projectile.cpp
	
triangulate.o : triangulate.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) triangulate.cpp
	
Line2D.o : Line2D.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) Line2D.cpp

clean:
	rm -rf *o $(TARGET)
