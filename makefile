CC= g++
CGLAG=  -ggdb -g -m64 -O0

all: main

main: main.o util.o my_curve.o my_vector.o lsh.o g_i.o h_i.o GridHash.o init_k++.o
	$(CC) $(CFLAG) -o main main.o util.o my_curve.o my_vector.o lsh.o g_i.o h_i.o GridHash.o init_k++.o

main.o: main.cpp
	$(CC) -c main.cpp

my_curve.o: my_curve.cpp my_curve.hpp my_vector.o
	$(CC) -c my_curve.cpp

my_vector.o: my_vector.cpp my_vector.hpp
	$(CC) -c my_vector.cpp

util.o: util.cpp util.hpp my_curve.o
	$(CC) -c util.cpp

lsh.o: lsh.cpp lsh.hpp my_vector.o my_curve.o util.o
	$(CC) -c lsh.cpp

GridHash.o: GridHash.cpp GridHash.hpp my_curve.o
	$(CC) -c GridHash.cpp

g_i.o: g_i.cpp g_i.hpp my_vector.o
	$(CC) -c g_i.cpp

h_i.o: h_i.cpp h_i.hpp my_vector.o
	$(CC) -c h_i.cpp

init_k++.o:	init_k++.hpp init_k++.cpp
		$(CC) $(CFLAG)	-c	init_k++.cpp


.PHONY: clean
clean:
	rm -f main main.o my_curve.o util.o my_vector.o lsh.o GridHash.o g_i.o h_i.o init_k++.o
