all: main.exe

main.exe: SurfaceClass.o BindingSiteClass.o main.o
	g++ -std=c++11 -gdwarf-3 -O3 SurfaceClass.o BindingSiteClass.o main.o -o main.exe

SurfaceClass.o: SurfaceClass.cpp BindingSiteClass.cpp SurfaceClass.h BindingSiteClass.h
	g++ -std=c++11 -gdwarf-3 -c -O3 -Wall SurfaceClass.cpp -o SurfaceClass.o

BindingSiteClass.o: BindingSiteClass.cpp BindingSiteClass.h
	g++ -std=c++11 -gdwarf-3 -c -O3 -Wall BindingSiteClass.cpp -o BindingSiteClass.o

main.o: main.cpp
	g++ -std=c++11 -gdwarf-3 -c -O3 -Wall main.cpp -o main.o

clean: 
	rm -rf main.exe SurfaceClass.o BindingSiteClass.o main.o
