CC=g++
CFLAGS=-std=c++11 -g3
exec=lsh cube lsh.o cube.o
.PHONY:clean

all: lsh cube

lsh:	lsh.o
	@echo " Compile lsh ...";
	$(CC)	$(CFLAGS)	 lsh.o 	-o lsh	;
	rm lsh.o;

cube:	cube.o
	@echo " Compile cube...";
	$(CC)	$(CFLAGS)	cube.o	-o cube;
	rm cube.o;

lsh.o:
	$(CC)	$(CFLAGS)	-c lsh.cpp
cube.o:
	$(CC)	$(CFLAGS) 	-c cube.cpp

clean:
	rm $(exec)





