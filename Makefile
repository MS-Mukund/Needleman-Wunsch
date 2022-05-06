# indicate how the object files are to be created
CC         := icc 
# CFLAGS     := -g -O3 -axCORE-AVX2  -qopenmp
CFLAGS     := -g -O3 -axCORE-AVX2 

NRUNS = 10

OBJECT_FILES := helper.o without_parallelisation.o
wparl: $(OBJECT_FILES)
	$(CC) $(CFLAGS) $(OBJECT_FILES) -o wparl

OBJECT_FILES := helper.o naive.o
naive: $(OBJECT_FILES)
	$(CC) $(CFLAGS) $(OBJECT_FILES) -o naive

naive.o: naive.c
	$(CC) $(CFLAGS) -c naive.c

without_parallelisation.o : without_parallelisation.c
	$(CC) $(CFLAGS) -c without_parallelisation.c

helper.o : helper.c
	$(CC) $(CFLAGS) -c helper.c

clean:
	rm *.o 