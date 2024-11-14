#!/bin/bash

# Command line arguments for run_sim are R, Dm, k2, alpha

ki=0.0001
tc=5.0
for kp in 0.00001 0.0001 0.001 0.01 0.1 1.0 5.0 10.0
do
        echo "kp ${kp} ki ${ki}"
        if [ -d kp_${kp}_ki_${ki}_tc_${tc} ]; then
                echo "exists, skipping"
        else
                mkdir kp_${kp}_ki_${ki}_tc_${tc}
                cd kp_${kp}_ki_${ki}_tc_${tc}
                cp ../pid .
                cp ../run_sim.sh .
                sbatch run_sim.sh ${kp} ${ki} ${tc}
                cd ../
        fi
done
