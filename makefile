# makefile
CC = g++
CFLAGS = -O3 -march=native -std=c++11

Simulate : main.cpp CASCL.cpp CASCL.h COO.cpp COO.h
	$(CC) $(CFLAGS) main.cpp CASCL.cpp COO.cpp -o simulate
