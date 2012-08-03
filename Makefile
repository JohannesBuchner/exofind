
CC := gcc
CFLAGS := -O3 -fPIC -std=c99 -Wall -Wextra -pedantic -lm
CFLAGS += -DINTERPOL=1 
LFLAGS += -L$(MULTINEST) -lnest3 -llapack -lpthread

exo.so: exo.o interpol.o rungekutta.o
	${CC} ${CFLAGS} ${LFLAGS} $^ -o $@ -shared

%.o: %.c
	${CC} ${CFLAGS} $^ -o $@ -c

testexo: testexo.c exo.o interpol.o rungekutta.o
	${CC} ${CFLAGS} ${LFLAGS} $^ -pg -g -o $@ -lgsl

