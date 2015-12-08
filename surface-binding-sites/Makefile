all: main.exe

main.exe: fcc100.o fcc110.o fcc111.o main.o
	g++ -std=c++11 -g -O3 fcc100.o fcc110.o fcc111.o main.o -o main.exe

fcc100.o: fcc100.cpp
	g++ -std=c++11 -g -c -O3 -Wall fcc100.cpp -o fcc100.o

fcc110.o: fcc110.cpp
	g++ -std=c++11 -g -c -O3 -Wall fcc110.cpp -o fcc110.o

fcc111.o: fcc111.cpp
	g++ -std=c++11 -g -c -O3 -Wall fcc111.cpp -o fcc111.o

main.o: main.cpp
	g++ -std=c++11 -g -c -O3 -Wall main.cpp -o main.o

clean: 
	rm -rf main.exe fcc100.o fcc110.o fcc111.o main.o 
