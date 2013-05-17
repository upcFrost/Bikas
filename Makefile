CPP=g++
CFLAGS=-c -g -Wall -Wno-deprecated -Wno-unused-but-set-variable -march=amdfam10
LDLIBS=-L. -L./lib -L/usr/lib -lvtkCommon -lvtkGraphics \
	-lvtkIO -lvtkFiltering -lvtkRendering -lvtkImaging

INCLUDES=-I. -I./include 

ifeq ($(OS),Windows_NT)
	INCLUDES+=-I"C:\Program Files (x86)\VTK\include\vtk-5.10" 
	LDLIBS+=-L"C:\Program Files (x86)\VTK\lib\vtk-5.10"
	TARGET=main.exe
else
	INCLUDES+=-I/usr/include/vtk-5.10
	LDLIBS+=-L./usr/lib 
	TARGET=main 
endif

OPENMP=-fopenmp

OBJECTS=main.o functions.o globals.o interp.o border.o \
		debug.o output.o
		
all : $(TARGET)

$(TARGET) : $(OBJECTS)
	$(CPP) -o bin/$(TARGET) $(OBJECTS) $(LDLIBS) $(OPENMP)


main.o : main.cpp globals.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) $(OPENMP) main.cpp

functions.o : functions.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) functions.cpp

globals.o : globals.cpp
	$(CPP) $(CFLAGS) $(INCLUDES) globals.cpp

interp.o : interp.c
	$(CPP) $(CFLAGS) $(INCLUDES) interp.c

debug.o : debug.c
	$(CPP) $(CFLAGS) $(INCLUDES) debug.c

border.o : border.c
	$(CPP) $(CFLAGS) $(INCLUDES) border.c

output.o : output.c
	$(CPP) $(CFLAGS) $(INCLUDES) output.c

clean:
	rm -rf *o $(TARGET)
