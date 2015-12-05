all: main.exe

main.exe: fcc100.o main.o
	g++ -std=c++11 fcc100.o main.o -o main.exe

fcc100.o: fcc100.cpp
	g++ -std=c++11 -c -O3 -Wall fcc100.cpp -o fcc100.o

main.o: main.cpp
	g++ -std=c++11 -c -O3 -Wall main.cpp -o main.o

clean: 
	rm -rf main.exe fcc100.o main.o 
