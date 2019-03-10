CC      =  gcc
CFLAGS  =  -g -Wall #-O2 
LDFLAGS = -lz

PROG = segmerge

.SUFFIXS:.c .o

all:$(PROG)

segmerge: opt.o paf.o sdict.o segmerge.o
	$(CC) $(CFLAGS) $^ -o $@  $(LDFLAGS) 

.c .o:
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o $(PROG)

	
