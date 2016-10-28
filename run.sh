#!/bin/bash          
echo "Compiling programs.." 

#################################################
# Compile

mpic++ cgm_mpsp14057.cpp -o cgm_mpsp14057
#################################################
echo "done!" 

echo "run the program with 4 processes!" 
#################################################
# Execute

mpirun -np 4 ./cgm_mpsp14057 800 4 2
#################################################
