/************************************************************************/
/* Christos Baziotis - mpsp14057                                        */
/************************************************************************/
/* Implementation of the algorithm in:                                  */
/* E. L. G. Saukas and S. W. Song. ,                                    */
/* "A Note on Parallel Selection on Coarse-Grained Multicomputers",     */
/* Algorithmica (1999)                                                  */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>

using namespace std;

#define  MASTER     0

// Creates an array of random numbers
int *create_rand_nums(int num_elements, int discrete) {
    int *rand_nums = new int[num_elements];
    assert(rand_nums != NULL);
    for (int i = 0; i < num_elements; i++) {
        int randnum = rand();
        if (discrete == 1 && i > 0) {
            int is_discrete = 1;
            do {
                for (int j = 0; j < i - 1; j++) {
                    if (rand_nums[i] == rand_nums[j]) {
                        is_discrete = 0;
                        break;
                    }
                }
                if (is_discrete == 1) {
                    rand_nums[i] = randnum;
                }
            } while (is_discrete == 0);
        }
        else {
            rand_nums[i] = randnum;
        }
    }
    return rand_nums;
}

void sort_array(int *array, int length) {
    for (int i = length - 1; i > 0; --i) {
        for (int j = 0; j < i; ++j) {
            if (array[j] > array[j + 1]) {
                int temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
            }
        }
    }
}

void print_array(int *array, int length) {
    printf("[");
    for (int i = 0; i < length; i++) {
        if (i == length - 1) {
            printf("%d", array[i]);
        }
        else {
            printf("%d,", array[i]);
        }
    }
    printf("]");
}

void write_array_to_file(int *array, int length, ofstream &file) {
    file << "[";
    for (int i = 0; i < length; i++) {
        if (i == length - 1) {
            file << array[i];
        }
        else {
            file << array[i] << ",";
        }
    }
    file << "]";
}

int get_array_sum(int *array, int length) {
    int sum = 0;
    for (int i = 0; i < length; i++) {
        sum += array[i];
    }
    return sum;
}

int get_array_median(int *array, int length) {

    if (length < 1) return 0;

    // Allocate an array of the same size and sort it.
    int *sorted = new int[length];
    for (int i = 0; i < length; ++i) {
        sorted[i] = array[i];
    }
    sort_array(sorted, length);

    if (length > 1) {
        // Find Median = ceil(n / 2)
        int m = 0;
        m = (int) ceil((double) length / 2);
        int median = sorted[m];
        free(sorted);
        return median;
    }
    else {
        int median = sorted[0];
        free(sorted);
        return median;
    }
}

/************************************************************************
* From R Documentation:Weighted Median Value
* For the n elements x    = c(x[1], x[2], ..., x[n])
* with positive weights w = c(w[1], w[2], ..., w[n]) such that sum(w) = S,
* the weighted median is defined as the element x[k]
* for which the total weight of all elements x[i] < x[k]
* is less or equal to S/2 and for which
* the total weight of all elements x[i] > x[k] is less or equal to S/2
************************************************************************/
int get_array_weighted_median(int *mis, int *nis, int length) {

    int N = get_array_sum(nis, length);

    // sort
    for (int i = length - 1; i > 0; --i) {
        for (int j = 0; j < i; ++j) {
            if (mis[j] > mis[j + 1]) {
                int temp = mis[j];
                mis[j] = mis[j + 1];
                mis[j + 1] = temp;
                temp = nis[j];
                nis[j] = nis[j + 1];
                nis[j + 1] = temp;
            }
        }
    }

    //calculate weights
    double S = 0; // weights sum
    double *weights = new double[length];
    for (int i = 0; i < length; ++i) {
        weights[i] = (double) nis[i] / (double) N;
        S += weights[i];
    }

    int k = 0;
    double sum = S - weights[0]; // sum is the total weight of all `x[i] > x[k]`
    while (sum > S / 2) {
        ++k;
        sum -= weights[k];
    }
    return mis[k];
}


int main(int argc, char **argv) {

    int rank, size;             // for storing this process' rank, and the number of processes
    int *elements = NULL;       // array containing the elements
    int *proc_elements = NULL;  // array containing the elements of each process
    int elements_count;         // elements count
    int k;                      // k-th element
    int c;                      // adjustable arbitrary parameter
    ofstream outfile;           // the file that the output will be written into

    MPI_Status status;

    /************************************************************************/
    /* initialization                                                       */
    /************************************************************************/

    printf("\n");
    outfile.open("output.txt");

    // read arguments
    if (argc != 4) {
        printf("incorrect number of arguments!\n");
        printf("please pass the necessary arguments in the following order:\n");
        printf("mpirun -np NUMER_OF_PROCESSES  ./program NUMBER_OF_ARRAY_ELEMENTS K-th_ELEMENT PARAM_C\n");
        return 0;
    }
    elements_count = atoi(argv[1]);
    k = atoi(argv[2]);
    c = atoi(argv[3]);

    if (k > elements_count) {
        printf("k cant't be larger than the array size!\n");
        return 0;
    }


    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (elements_count % size != 0) {
        printf("number of array elements NOT dividable by number of processes!\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        return 0;
    }


    // ...good to go!

    if (rank == MASTER) {
        srand(time(NULL));
        elements = create_rand_nums(elements_count, 1);
        printf("\n*****************************************************\n");
        printf("Random array\n");
        print_array(elements, elements_count);
        printf("\n*****************************************************\n");
        printf("\n");
        fflush(stdout); //flush stdout output regularly to avoid wrong output order

        outfile << "NUMER_OF_PROCESSES \t \t" << size << endl;
        outfile << "NUMBER_OF_ARRAY_ELEMENTS \t" << elements_count << endl;
        outfile << "K \t \t \t \t" << k << endl;
        outfile << "C \t \t \t \t" << c << endl;
        outfile << endl << endl;
        outfile << "Random array => ";
        write_array_to_file(elements, elements_count, outfile);
    }

    //process elements count
    int ni = elements_count / size;

    //the array chunk of each process
    proc_elements = new int[ni];
    //printf("proc_elements_count %d => %d", rank, ni);

    // Scatter the random numbers to all processes
    MPI_Scatter(elements, ni, MPI_INT, proc_elements, ni, MPI_INT, MASTER, MPI_COMM_WORLD);

    // no longer needed
    free(elements);

    /************************************************************************/
    /* (1) Set N:=n															                                              */
    /************************************************************************/
    int N = elements_count;

    /************************************************************************/
    /* (2) Repeat until N<= n/c*p											                                         */
    /************************************************************************/
    int M;
    int round = 0;

    int *mis = NULL;    // array of medians of each process
    int *nis = NULL;    // array of elements count of each process

    //while (N > (double) elements_count / (double) (c * size))
    while (!(N <= (double) elements_count / (double) (c * size))) {
        
		round++;
		
		
        if (rank == MASTER) {
            printf("\n*****************************************************\n");
            printf("Round %d\n\n", round);
            // printf("Round %d, N = %d, n= %d, (c * size)=%d \n\n", round, N, elements_count, (c * size));
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);


        /************************************************************************/
        /* (2.1) Each processor i computes the median mi of its ni elements     */
        /************************************************************************/
        int mi = get_array_median(proc_elements, ni);    // Median of process elements

        printf("p:%d => ", rank);
        print_array(proc_elements, ni);
        printf(", median:%d\n", mi);
        fflush(stdout);


        /************************************************************************/
        /* (2.2) Each processor i sends mi and ni to processor 1		        */
        /************************************************************************/
        if (rank == MASTER) {
            mis = new int[size];
            nis = new int[size];
        }
        MPI_Gather(&mi, 1, MPI_INT, mis, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
        MPI_Gather(&ni, 1, MPI_INT, nis, 1, MPI_INT, MASTER, MPI_COMM_WORLD);


        /************************************************************************/
        /* (2.3) Processor 1 computes the weighted median M				        */
        /************************************************************************/
        int weighted_median;
        if (rank == MASTER) {
            weighted_median = get_array_weighted_median(mis, nis, size);
        }


        /************************************************************************/
        /* (2.4) Processor 1 broadcasts M to all other processors		        */
        /************************************************************************/
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&weighted_median, 1, MPI_INT, MASTER, MPI_COMM_WORLD);


        /************************************************************************/
        /* (2.5) Each processor i computes li, ei, gi, respectively the numbers */
        /* of its local elements less than, equal to, or greater than M         */
        /************************************************************************/
        int local_leg[3] = {0, 0, 0};//li, ei, gi
        for (int i = 0; i < ni; i++) {
            if (proc_elements[i] < weighted_median)
                local_leg[0]++;
            else if (proc_elements[i] == weighted_median)
                local_leg[1]++;
            else
                local_leg[2]++;
        }
        //printf("\n PRO=%d   L=%d E=%d G=%d \n", rank, local_leg[0], local_leg[1], local_leg[2]);


        /************************************************************************/
        /* (2.6) Each processor i sends li, ei, gi to processor 1               */
        /* (2.7) Processor 1 computes L, E, G, respectively the total numbers   */
        /*		 of elements less than, equal to, or greater than M             */
        /* (2.8) Processor 1 broadcasts L, E, G, to all other processors		*/
        /************************************************************************/
        int global_leg[3] = {0, 0, 0};
        MPI_Allreduce(&local_leg, &global_leg, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        /*MPI_Allreduce
        Many parallel applications will require accessing the reduced results
        across all processes rather than the root process. In a similar complementary
        style of MPI_Allgather to MPI_Gather, MPI_Allreduce will reduce the values
        and distribute the results to all processes.*/


        /************************************************************************/
        /* (2.9) Check conditions										        */
        /************************************************************************/
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == MASTER) {
            printf("\n");
            printf("\n");
            printf("L=%d E=%d G=%d, M=%d \n\n", global_leg[0], global_leg[1], global_leg[2], weighted_median);
            fflush(stdout);

            if (global_leg[0] < k && k <= (global_leg[0] + global_leg[1])) {
                printf("Condition 1: break loop with M:%d\n", weighted_median);
            }
            else if (k < global_leg[0]) {
                printf("Condition 2: each processor discards all but those elements less than M and set N:=L\n");
            }
            else if (k > (global_leg[0] + global_leg[1])) {
                printf("Condition 3: each processor discards all but those elements greater than M\n");
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);

        M = weighted_median;

        if (global_leg[0] < k && k <= (global_leg[0] + global_leg[1])) {
            break;
        }
            // each processor i discards all but those elements less than M and set N:=L
        else if (k < global_leg[0]) {
            int *temp_array = new int[local_leg[0]];
            int j = 0;
            for (int i = 0; i < ni; i++) {
                if (proc_elements[i] < weighted_median) {
                    temp_array[j] = proc_elements[i];
                    j++;
                }
            }
            ni = local_leg[0];
            free(proc_elements);
            proc_elements = temp_array;

            N = global_leg[0];
        }
            // each processor i discards all but those elements greater than M
        else if (k > (global_leg[0] + global_leg[1])) {
            int *temp_array = new int[local_leg[2]];
            int j = 0;
            for (int i = 0; i < ni; i++) {
                if (proc_elements[i] > weighted_median) {
                    temp_array[j] = proc_elements[i];
                    j++;
                }
            }
            ni = local_leg[2];
            free(proc_elements);
            proc_elements = temp_array;

            N = global_leg[2];
            k = k - (global_leg[0] + global_leg[1]);
        }

        printf("p:%d => ", rank);
        print_array(proc_elements, ni);
        printf("\n");
        fflush(stdout);
    }

    if (rank == MASTER) {
        nis = new int[size];
    }
    MPI_Gather(&ni, 1, MPI_INT, nis, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

    int *displs = new int[size];
    int s;
    if (rank == MASTER) {
        displs[0] = 0;
        for (int i = 1; i < size; i++) {
            displs[i] = nis[i - 1] + displs[i - 1];
        }
        s = get_array_sum(nis, size);
        elements = new int[s];
    }

    fflush(stdout);

    MPI_Gatherv(proc_elements, ni, MPI_INT, elements, nis, displs, MPI_INT, MASTER, MPI_COMM_WORLD);

    // wait for other processes
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == MASTER) {
        sort_array(elements, s);
        printf("\n*****************************************************\n");
        /*printf("ELEMENTS ");
        int s = get_array_sum(nis, size);
        print_array(elements, s);*/
        printf("\n");
        printf("\n");
        printf("Solution (%d) => %d !!!", k, elements[k - 1]);
        outfile << endl << "Solution (" << k << ") => " << elements[k - 1];
        printf("\n");
        fflush(stdout);
    }

    free(elements);
    free(mis);
    free(nis);
    free(proc_elements);
    outfile.close();

    MPI_Finalize();
    return 0;
}
