CC=gcc-9
DEBUG=-ggdb
CFLAGS= -Wall -std=c99
LMFLAGS=-lm -lfftw3 -lprofiler
TARGET=spinodal_spect
LIBRARIES=-I /usr/local/include/
LINKER=-L /usr/local/lib/

DEPENDENCIES= array_allocate.c \
							array_deallocate.c \
							prep_microstructure.c \
							rand_ZeroToOne.c \
							write_to_VTK.c

$(TARGET): $(TARGET).c $(DEPENDENCIES)
	$(CC) $(LMFLAGS) $(TARGET).c $(DEPENDENCIES) -o $(TARGET).out $(CFLAGS) $(LIBRARIES) $(LINKER)

debug: $(TARGET).c $(DEPENDENCIES)
	$(CC) $(LMFLAGS) $(DEBUG) $(TARGET).c $(DEPENDENCIES) -o $(TARGET).out $(CFLAGS) $(LIBRARIES) $(LINKER)

time: $(TARGET).c $(DEPENDENCIES)
		$(CC) $(LMFLAGS) $(TARGET).c $(DEPENDENCIES) -o $(TARGET).out $(CFLAGS) $(LIBRARIES) $(LINKER)
		CPUPROFILE=$(TARGET).prof ./$(TARGET).out
		pprof --functions --text $(TARGET).out $(TARGET).prof > $(TARGET).txt

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
	-rm -f $(TARGET).prof
	-rm -f $(TARGET).pdf
