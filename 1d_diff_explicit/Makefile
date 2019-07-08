CC=gcc
CFLAGS= -Wall
TARGET=Diff1DExplicit

$(TARGET): $(TARGET).c write_to_VTK.c
	$(CC) $(TARGET).c write_to_VTK.c -o $(TARGET).out $(CFLAGS)

.PHONY: clean
	clean:
	-rm -f $(TARGET).o
	-rm -f $(TARGET)
