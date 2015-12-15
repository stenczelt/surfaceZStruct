all: main.exe

main.exe: SurfaceClass.o BindingSiteClass.o main.o
	g++ -std=c++11 -gdwarf-3 SurfaceClass.o BindingSiteClass.o main.o -o main.exe
	#g++ -std=c++11 -g -O3 SurfaceClass.o BindingSiteClass.o main.o -o main.exe

SurfaceClass.o: SurfaceClass.cpp
	#g++ -std=c++11 -g -c -O3 -Wall SurfaceClass.cpp -o SurfaceClass.o
	g++ -std=c++11 -gdwarf-3 -c -Wall SurfaceClass.cpp -o SurfaceClass.o

BindingSiteClass.o: BindingSiteClass.cpp
	#g++ -std=c++11 -g -c -O3 -Wall BindingSiteClass.cpp -o BindingSiteClass.o
	g++ -std=c++11 -gdwarf-3 -c -Wall BindingSiteClass.cpp -o BindingSiteClass.o

main.o: main.cpp
	#g++ -std=c++11 -g -c -O3 -Wall main.cpp -o main.o
	g++ -std=c++11 -gdwarf-3 -c -Wall main.cpp -o main.o

clean: 
	rm -rf main.exe SurfaceClass.o BindingSiteClass.o main.o
