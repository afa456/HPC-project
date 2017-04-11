/**
 * @file    parallel_sort.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Declares the parallel sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef PARALLEL_SORT_H
#define PARALLEL_SORT_H

#include <mpi.h>

/**
 * @brief   Parallel, distributed sorting over all processors in `comm`. Each
 *          processor has the local input [begin, end).
 *
 * Note that `end` is given one element beyond the input. This corresponds to
 * the API of C++ std::sort! You can get the size of the local input with:
 * int local_size = end - begin;
 *
 * @param begin Pointer to the first element in the input sequence.
 * @param end   Pointer to one element past the input sequence. Don't access this!
 * @param comm  The MPI communicator with the processors participating in the
 *              sorting.
 */
void parallel_sort(int * begin, int* end, MPI_Comm comm);


/*********************************************************************
 *              Declare your own helper functions here               *
 *********************************************************************/

/**
 * @brief   Parallel sorting the array on all processors in `comm`. Each
 *          processor has the local input size. However, the size of 
 *          the local array may change. The new size is returned by the function.    
 * @param begin      Pointer to the first element in the input sequence.
 * @param arraysize  Local array length
 * @param comm       The MPI communicator with the processors participating in the
 *                   sorting.
 *
 * @return The new length of local array.
 */ 

int Quick_sorting(int* &begin, int arraysize, MPI_Comm comm);

/**
 * @brief   Partitions the given array about the pivot.
 * @param  begin     		 Pointer to the first element in the input sequence.
 * @param  end       		 Pointer to last element in the input sequence
 * @param  pivot_element     Pivot value.
 */

int Partition_Array(int* begin, int* end, int pivot_element);

/**
 * @brief   Partition the processors.
 * @param procLeft   Number of processor for the left-part array
 * @param procRight  Number of processor for the right-part array
 * @param numProc    Total number of processor  
 * @param numLeft    Number of left part elements
 * @param numRight   Number of right part elements
 * @param comm       The MPI communicator with the processors participating in the
 *                   sorting.
 */

void Partition_Processor(int &procLeft, int &procRight, int numprocs, int numLeft, int numRight, MPI_Comm comm);

/**
 * @brief   Distributes the partitioned data in each processor to the processors assigned to the partition.
 *
 * @param sendcnts     Stores the number of elements received from Processor i in ith entry; 
 * @param senddis      Stores the displacement of the sendbuf for Processor i in ith entry;
 * @param recvcnts     Receive Counts, Similar to sendcnts;
 * @param recvdis      Receive Displacement, Similar to senddis;
 * @param curlen      Current size of the local array
 * @param objectSize   Final local array size at each processor
 * @param comm         The MPI communicator with the processors participating in the
 *                     sorting.
 */

void Distribute_data(int* sendcnts, int* senddis, int* recvcnts, int* recvdis, int curlen, int objlen, MPI_Comm &comm);



#endif // PARALLEL_SORT_H
