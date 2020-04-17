all: exacteulerrp

exacteulerrp: exacteulerrp.o
	g++ exacteulerrp.o -o exacteulerrp

exacteulerrp.o: exacteulerrp.cc
	g++ -std=c++14 -O3 -Wall -c exacteulerrp.cc

clean:
	rm -rf exacteulerrp exacteulerrp.o
