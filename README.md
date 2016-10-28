# cgm
Parallel MPI implementation of the algorithm presented in [1], for determining the kth smallest element of a set of n elements.
KeyWords: Coarse-grained multicomputer, Parallel algorithm, Selection problem.

## Execution Example

```
mpirun -np NUMER_OF_PROCESSES ./cgm NUMBER_OF_ARRAY_ELEMENTS K-th_ELEMENT PARAM_C
```
* NUMER_OF_PROCESSES: The number of processes we wish to use
* NUMBER_OF_ARRAY_ELEMENTS: The number of elements in the array we will be searching
* K-th_ELEMENT: The ordered position of the element in the array 
* PARAM_C: a parameter the tunes the behavior of the algorithm (read the paper)

So for example if we wish to run the program using 4 processes and find the 4th element in an array of 800 elements with c=2:
```
mpirun -np 4 ./cgm 800 4 2
```

[1]: E. L. G. Saukas and S. W. Song. ,  "A Note on Parallel Selection on Coarse-Grained Multicomputers",  Algorithmica (1999)                                          
