CC=gcc
CFLAGS= -Wall

Diff1DExplicit: Diff1DExplicit.c write_to_VTK.c
	$(CC) Diff1DExplicit.c write_to_VTK.c -o Diff1DExplicit.out $(CFLAGS)

clean:
	-rm -f Diff1DExplicit.o
	-rm -f Diff1DExplicit
