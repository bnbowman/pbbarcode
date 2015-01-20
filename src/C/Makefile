.PHONY: clean all
SHELL = /bin/bash -e

all: mkdir build/sw.so build/sw_new.so

mkdir:
	mkdir -p ./build

build/sw.so: sw.c
	gcc -O4 -DGETPROB -shared -fPIC sw.c -o build/sw.so;

build/sw_new.so: sw_new.c
	gcc -O4 -DGETPROB -shared -fPIC sw_new.c -o build/sw_new.so;

clean:
	rm -rf build

