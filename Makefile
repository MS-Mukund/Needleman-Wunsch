# indicate how the object files are to be created

#icc
CC         := icc
CFLAGS     := -g -O3 -axCORE-AVX2 -qopenmp
# CFLAGS     := -O3 -axCORE-AVX2

#gcc O0
# CC         := gcc 
# CFLAGS     := -O0 -qopenmp
# CFLAGS     := -O0 

#gcc O1
# CC         := gcc 
# CFLAGS     := -O1 -qopenmp
# CFLAGS     := -O1 

#gcc O2
# CC         := gcc 
# CFLAGS     := -O2 -qopenmp
# CFLAGS     := -O2 

#gcc O3
# CC         := gcc 
# CFLAGS     := -O3 -fopenmp
# CFLAGS     := -O3 

NRUNS = 10

OBJECT_FILES := helper.o optimised.o
opti: $(OBJECT_FILES)
	$(CC) $(CFLAGS) $(OBJECT_FILES) -o opti

# anti-diagonal, without tiling
# OBJECT_FILES := helper.o without_parallelisation.o
# wparl: $(OBJECT_FILES)
# $(CC) $(CFLAGS) $(OBJECT_FILES) -o wparl

#naive implementation
# OBJECT_FILES := helper.o naive.o 
# naive: $(OBJECT_FILES)
# $(CC) $(CFLAGS) $(OBJECT_FILES) -o naive

# anti-diagonal, with tiling
# OBJECT_FILES := helper.o tiled.o
# tiled: $(OBJECT_FILES)
# $(CC) $(CFLAGS) $(OBJECT_FILES) -o tiled

optimised.o: optimised.c
	$(CC) $(CFLAGS) -c optimised.c
tiled.o: tiled.c
	$(CC) $(CFLAGS) -c tiled.cnaive.o: naive.c
	$(CC) $(CFLAGS) -c naive.c

without_parallelisation.o : without_parallelisation.c
	$(CC) $(CFLAGS) -c without_parallelisation.c

helper.o : helper.c
	$(CC) $(CFLAGS) -c helper.c

clean:
	rm *.o 