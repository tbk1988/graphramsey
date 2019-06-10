# 

NAUTYDIR=/home/okr/packages/nauty26r7
VNAUTYDIR=/home/okr/packages/very_nauty-1.1
CC=gcc 
#CFLAGS= -O3 -DNDEBUG -I$(VNAUTYDIR) -I$(NAUTYDIR) -I/usr/include/igraph -I./
CFLAGS= -Wall -g -I$(VNAUTYDIR) -I$(NAUTYDIR) -I/usr/include/igraph -I./
LDFLAGS=-L$(VNAUTYDIR) -L$(NAUTYDIR) -L/usr/local/lib 

all: edge_interval #numberofC4 independence alphacounts minindependence shearertest edge_interval #iv1test #local

local: localstructuretable.c
	${CC} ${CFLAGS} -o local localstructuretable.c \
	$(NAUTYDIR)/naututil.o \
	$(NAUTYDIR)/nauty.o \
	$(NAUTYDIR)/nautil.o \
	$(NAUTYDIR)/naugraph.o \
	$(NAUTYDIR)/schreier.o $(NAUTYDIR)/naurng.o \
	$(NAUTYDIR)/gtools.o \
	$(NAUTYDIR)/nausparse.o ${LDFLAGS} -lvn_graph -lm

numberofC4: numberofC4.c graph6_utils.o gutils.o
	gcc $(CFLAGS) numberofC4.c graph6_utils.o gutils.o -ligraph \
	-o numberofC4

independence: independence.c graph6_utils.o gutils.o
	gcc $(CFLAGS) independence.c graph6_utils.o gutils.o -ligraph -o independence 

shearertest: shearertest.c graph6_utils.o gutils.o
	gcc $(CFLAGS) shearertest.c graph6_utils.o gutils.o -ligraph -lgmp -o shearertest 

iv1test: iv1test.c graph6_utils.o gutils.o
	gcc $(CFLAGS) iv1test.c graph6_utils.o gutils.o -ligraph -lgmp -o iv1test 

edge_interval: edge_interval.c
	gcc $(CFLAGS) edge_interval.c -lgmp -o edge_interval

minindependence: minindependence.c graph6_utils.o 
	gcc $(CFLAGS) minindependence.c graph6_utils.o -ligraph -o minindependence 

alphacounts: alphacounts.c graph6_utils.o 
	gcc $(CFLAGS) alphacounts.c graph6_utils.o -ligraph -o alphacounts

graph6_utils.o: graph6_utils.c
	${CC} -c ${CFLAGS} graph6_utils.c -lligraph -o graph6_utils.o

gutils.o: gutils.c
	${CC} -c ${CFLAGS} gutils.c -ligraph -o gutils.o

clean:
	rm -rf local numberofC4 graph6_utils.o independence gutils.o \
	alphacounts shearertest minindependence iv1test
