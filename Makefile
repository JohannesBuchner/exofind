
CC := gcc
CFLAGS := -O3 -fPIC -std=c99 -Wall -Wextra -pedantic -lm
CFLAGS += -DINTERPOL=1 
LFLAGS += -L$(MULTINEST) -lnest3 -latllapack -lpthread

all: exo.so exo_noprior.so testexo

exo.so: exo.o interpol.o rungekutta.o
	${CC} -DUSE_PRIORS -DNOISE ${CFLAGS} ${LFLAGS} $^ -o $@ -shared
exo_noprior.so: exo.o interpol.o rungekutta.o
	${CC} -DNOISE ${CFLAGS} ${LFLAGS} $^ -o $@ -shared

%.o: %.c
	${CC} ${CFLAGS} $^ -o $@ -c

testexo: testexo.c exo.o interpol.o rungekutta.o
	${CC} ${CFLAGS} ${LFLAGS} $^ -pg -g -o $@ -lgsl


