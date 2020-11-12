# Makefile for LMBM-B
FF = gfortran -O2 -fPIC
CC = gcc

all:	shared_lib

#tlmbmb.o: tlmbmb.f
#	$(FF) -c tlmbmb.f

lmbmb.o: lmbmb.f
	$(FF) -c lmbmb.f

lmbmbs.o: lmbmbs.f
	$(FF) -c lmbmbs.f

lmbmbmex.o: lmbmbmex.f
	$(FF) -c lmbmbmex.f

matca2.o: matca2.f
	$(FF) -c matca2.f

objfunc.o: objfunc.c
	$(CC) -c objfunc.c -o objfunc.o

clean:
	rm *.o *.so

shared_lib: lmbmb.o lmbmbs.o lmbmbmex.o matca2.o objfunc.o
	gfortran -O3 -fPIC -shared objfunc.o lmbmb.o lmbmbs.o lmbmbmex.o matca2.o -o liblmbmb.so
