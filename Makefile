CC=gcc
CFLAGS=-c -Wall -Werror -O3 -ffast_math -march=native
COPTS=
LDFLAGS=-lz -lm
SOURCES=seqPrep.c utils.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=SeqPrep

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) ${COPTS} $(LDFLAGS) $(OBJECTS) -o $@

install: all
	-cp $(EXECUTABLE) $(HOME)/bin

.c.o:
	$(CC) ${COPTS} $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS) $(EXECUTABLE)

check-syntax:
	$(CC) ${CFLAGS} -o .nul -S ${CHK_SOURCES}