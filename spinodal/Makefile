CC=gcc
CFLAGS= -Wall
TARGET=spinodal
DEPENDENCIES= free_energy.c \
							laplacian.c \
							rand_ZeroToOne.c \
							solve.c \
							solve2.c \
							write_to_VTK.c

$(TARGET): $(TARGET).c $(DEPENDENCIES)
	$(CC) $(TARGET).c $(DEPENDENCIES) -o $(TARGET).out $(CFLAGS)

clean:
	-rm -f $(TARGET).o
	-rm -f $(TARGET)
