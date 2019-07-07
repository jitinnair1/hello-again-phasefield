CC=gcc
CFLAGS= -Wall
TARGET=spinodal

$(TARGET): $(TARGET).c write_to_VTK.c
	$(CC) $(TARGET).c write_to_VTK.c -o $(TARGET).out $(CFLAGS)

clean:
	-rm -f $(TARGET).o
	-rm -f $(TARGET)
