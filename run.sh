#!/bin/bash          
echo "Compiling programs.." 

#################################################
# Compile

mpic++ cgm.cpp -o cgm
#################################################
echo "done!" 

echo "run the program with 4 processes!" 
#################################################
# Execute

mpirun -np 4 ./cgm 800 4 2
#################################################
