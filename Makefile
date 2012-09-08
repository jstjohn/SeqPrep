CC=gcc
CFLAGS=-c -Wall -Werror -O0 -g
#recommended options: -ffast-math -ftree-vectorize -march=core2 -mssse3 -O3
COPTS=
LDFLAGS=-lz -lm
SOURCES=SeqPrep.c utils.c stdaln.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=SeqPrep

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) ${COPTS} $(OBJECTS) $(LDFLAGS) -o $@

install: all
	-cp $(EXECUTABLE) $(HOME)/bin

.c.o:
	$(CC) ${COPTS} $(CFLAGS) $< -o $@

clean:
	-rm $(OBJECTS) $(EXECUTABLE)

check-syntax:
	$(CC) ${CFLAGS} -o .nul -S ${CHK_SOURCES}
