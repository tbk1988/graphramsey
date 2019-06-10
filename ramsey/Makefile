CC=gcc 
#CFLAGS= -O3 -DNDEBUG -I./
CFLAGS= -Wall -g -I./ 

all: edge_interval

edge_interval: edge_interval.c
	gcc $(CFLAGS) edge_interval.c -lgmp -o edge_interval 

clean:
	rm -rf edge_interval
