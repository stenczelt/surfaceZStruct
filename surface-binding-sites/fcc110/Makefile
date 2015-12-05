all: main.exe

main.exe: fcc110.o main.o
	g++ -std=c++11 -g -O3 fcc110.o main.o -o main.exe

fcc110.o: fcc110.cpp
	g++ -std=c++11 -g -c -O3 -Wall fcc110.cpp -o fcc110.o

main.o: main.cpp
	g++ -std=c++11 -g -c -O3 -Wall main.cpp -o main.o

clean: 
	rm -rf main.exe fcc110.o main.o 
