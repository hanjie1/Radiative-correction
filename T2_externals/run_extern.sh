#! /usr/bin/expect

array set target {H1 D2 He3 H3}
array set kin {0 1 2 3 4 5 7 9 11 13 15}

for i in `seq 2 2`
do
    for j in `seq 0 1`
    do
        ./externals_all
        expect " Enter the input file name (in INP directory)"
        send " ${target[i]}_kin${kin[j]}.inp"
         
    done

done
