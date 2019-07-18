CC=gcc
DEBUG=-ggdb
CFLAGS= -Wall -std=c99
TARGET=spinodal
DEPENDENCIES= array_deallocate.c \
							free_energy.c \
							laplacian.c \
							prep_microstructure.c \
							rand_ZeroToOne.c \
							set_array_zero.c \
							solve.c \
							solve2.c \
							write_to_VTK.c

$(TARGET): $(TARGET).c $(DEPENDENCIES)
	$(CC) $(TARGET).c $(DEPENDENCIES) -o $(TARGET).out $(CFLAGS)

debug: $(TARGET).c $(DEPENDENCIES)
	$(CC) $(DEBUG) $(TARGET).c $(DEPENDENCIES) -o $(TARGET).out $(CFLAGS)

clean:
	-rm -rf $(TARGET).out.*
	-rm -f $(TARGET).out
	-rm -f $(TARGET).o
	-rm -f $(TARGET)

purge:
	-rm -rf $(TARGET).out.*
	-rm -f $(TARGET).out
	-rm -f $(TARGET).o
	-rm -f $(TARGET)
	-rm -rf output/
