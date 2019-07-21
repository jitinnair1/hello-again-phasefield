CC=gcc
DEBUG=-ggdb
CFLAGS= -Wall -std=c99
LMFLAGS=-lm -lfftw3
TARGET=spinodal_spect
LIBRARIES=-iquote /usr/local/include
LINKER=-L/usr/local/lib
DEPENDENCIES= array_allocate.c \
							prep_microstructure.c \
							rand_ZeroToOne.c


$(TARGET): $(TARGET).c $(DEPENDENCIES)
	$(CC) $(LMFLAGS) $(TARGET).c $(DEPENDENCIES) -o $(TARGET).out $(CFLAGS) $(LIBRARIES) $(LINKER)

debug: $(TARGET).c $(DEPENDENCIES)
	$(CC) $(LMFLAGS) $(DEBUG) $(TARGET).c $(DEPENDENCIES) -o $(TARGET).out $(CFLAGS) $(LIBRARIES) $(LINKER)

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
