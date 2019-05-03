#! /usr/bin/expect

target=("H1" "D2" "He3" "H3")
kin=(0 1 2 3 4 5 7 9 11 13 15)

for i in `seq 0 3`
do
    for j in `seq 0 10`
    do
#	if [ "${kin[$j]}" -lt 7 ] || ["${kin[$j]}" -eq 15 ]
#	then
#           ./externals_all <<< "INP/${target[$i]}_kin${kin[$j]}.inp"
           ./externals_all <<< "${target[$i]}_kin${kin[$j]}.inp"
#	else
#           ./externals_all <<< "INP/ACC/${target[$i]}_kin${kin[$j]}_1st.inp"
#	   ./externals_all <<< "INP/ACC/${target[$i]}_kin${kin[$j]}_2nd.inp"
#	fi
#        expect " Enter the input file name (in INP directory)"
#        send "INP/ACC/${target[$i]}_kin${kin[$j]}.inp"
         
    done

done
