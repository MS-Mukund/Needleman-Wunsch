make clean;
make opti; 
for i in {100..5000..100}; 
do 
    for j in {1..10};
    do 
        ./opti ${i};
    done
done