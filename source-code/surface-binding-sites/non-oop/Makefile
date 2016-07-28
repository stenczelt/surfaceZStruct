all: main.exe

main.exe: bcc100.o bcc110.o bcc111.o hcp0001.o diamond100.o fcc100.o fcc110.o fcc111.o main.o
	g++ -std=c++11 -g -O3 bcc100.o bcc110.o bcc111.o hcp0001.o diamond100.o fcc100.o fcc110.o fcc111.o main.o -o main.exe

bcc100.o: bcc100.cpp
	g++ -std=c++11 -g -c -O3 -Wall bcc100.cpp -o bcc100.o

bcc110.o: bcc110.cpp
	g++ -std=c++11 -g -c -O3 -Wall bcc110.cpp -o bcc110.o

bcc111.o: bcc111.cpp
	g++ -std=c++11 -g -c -O3 -Wall bcc111.cpp -o bcc111.o

hcp0001.o: hcp0001.cpp
	g++ -std=c++11 -g -c -O3 -Wall hcp0001.cpp -o hcp0001.o

diamond100.o: diamond100.cpp
	g++ -std=c++11 -g -c -O3 -Wall diamond100.cpp -o diamond100.o

fcc100.o: fcc100.cpp
	g++ -std=c++11 -g -c -O3 -Wall fcc100.cpp -o fcc100.o

fcc110.o: fcc110.cpp
	g++ -std=c++11 -g -c -O3 -Wall fcc110.cpp -o fcc110.o

fcc111.o: fcc111.cpp
	g++ -std=c++11 -g -c -O3 -Wall fcc111.cpp -o fcc111.o

main.o: main.cpp
	g++ -std=c++11 -g -c -O3 -Wall main.cpp -o main.o

clean: 
	rm -rf main.exe bcc100.o bcc110.o bcc111.o hcp0001.o diamond100.o fcc100.o fcc110.o fcc111.o main.o 
