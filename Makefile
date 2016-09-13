CC=gcc
SRCDIR= src
CFLAGS=-Wall -O3 -std=c99 -g  -lm -lfftw3 -lmpi
CFLAGS2=-c
LDFLAGS=
SOURCES=main.c	cal_miu_HS_ex.c  iteration.c  namelist.c  omega.c  parallel.c  parameter.c  setup.c  vext.c

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=cdft

all:
	$(CC) $(SRCDIR)/*.c -o $(EXECUTABLE) $(CFLAGS) 
	rm -rf *.o *~


verbose: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
	rm -rf *.o *~

.c.o:
	$(CC) $(CFLAGS2) $(CFLAGS) $< -o $@


clean:
	rm -rf *.o *~
	rm -f $(EXECUTABLE)
