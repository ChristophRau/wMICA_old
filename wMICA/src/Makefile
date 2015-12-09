all: ICMg.iterations.so

ICMg.iterations.so: ICMg.iterations.c
	MAKEFLAGS='CFLAGS=-O2\ -Wall' R CMD SHLIB ICMg.iterations.c
	
clean:
	
