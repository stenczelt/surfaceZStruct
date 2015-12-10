all: main.exe

#main.exe: Termination100.o Termination111.o bcc110.o bcc111.o diamond100.o fcc110.o main.o
#	g++ -std=c++11 -g -O3 Termination100.o Termination111.o bcc110.o bcc111.o diamond100.o fcc110.o main.o -o main.exe

main.exe: Termination100.o Termination111.o bcc110.o bcc111.o fcc110.o main.o
	g++ -std=c++11 -g -O3 Termination100.o Termination111.o bcc110.o bcc111.o fcc110.o main.o -o main.exe

Termination100.o: Termination100.cpp
	g++ -std=c++11 -g -c -O3 -Wall Termination100.cpp -o Termination100.o

Termination111.o: Termination111.cpp
	g++ -std=c++11 -g -c -O3 -Wall Termination111.cpp -o Termination111.o

bcc110.o: bcc110.cpp
	g++ -std=c++11 -g -c -O3 -Wall bcc110.cpp -o bcc110.o

bcc111.o: bcc111.cpp
	g++ -std=c++11 -g -c -O3 -Wall bcc111.cpp -o bcc111.o

#diamond100.o: diamond100.cpp
#	g++ -std=c++11 -g -c -O3 -Wall diamond100.cpp -o diamond100.o

fcc110.o: fcc110.cpp
	g++ -std=c++11 -g -c -O3 -Wall fcc110.cpp -o fcc110.o

main.o: main.cpp
	g++ -std=c++11 -g -c -O3 -Wall main.cpp -o main.o

clean: 
	rm -rf main.exe Termination100.o Termination111.o bcc110.o bcc111.o fcc110.o main.o
#	rm -rf main.exe Termination100.o Termination111.o bcc110.o bcc111.o diamond100.o fcc110.o main.o
