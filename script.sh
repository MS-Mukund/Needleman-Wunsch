make clean;
make opti; 
for i in {1000..30000..1000}; 
do 
    # for j in {1..3};
    # do 
    ./opti ${i};
    # done
done