CC = g++

CFLAGS = -g
INCLUDES = -I/usr/include
LIBS = -L/usr/lib -lm -lgsl -lgslcblas
OBJECTFILES = def.o crystal.o

# Rule to create .o files from .c files
.c.o:
	$(RM) $@
	$(CC) -c $(CFLAGS) $(INCLUDES) $*.c

#Targets:

all:: crystal

#Dependencies
def.o: LR_vector.h def.cpp
crystal.o: crystal.cpp

#Executable
crystal: $(OBJECTFILES)
		 $(RM) $@
		 $(CC) $(CFLAGS) $(INCLUDES) $(OBJECTFILES) $(LIBS) -o $@
		 rm -rf *.o

clean:
	rm -rf *.o crystal