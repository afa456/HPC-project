/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "parallel_sort.h"
#include "utils.h"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <ctime>

// implementation of your parallel sorting
void parallel_sort(int* begin, int* end, MPI_Comm comm) {
    
    //Basic Information
    int numProc;
    int localSize = end - begin;
    MPI_Comm_size(comm, &numProc);

    //Copy [begin, end) to 'cache'
    int* cache = new int[localSize];
    for(int i = 0; i < localSize; i++)
      cache[i] = *(begin+i);

    std::srand(time(NULL));
    /*
    int myRank;
    MPI_Comm_rank(comm, &myRank);
    std::cout << "Rank: " << myRank << std::endl;
    for(int i = 0; i < localSize; i++)
      std::cout << begin[i] << "\t";
    std::cout << std::endl;
    */

    //Sort 'cache' with size-changeable sorting function parallel_sort_changeable.
    int newLocalSize = parallel_sort_changeable(cache,localSize,comm);

    /*
    std::cout << "OK" ;
    */

    //Compute the configure for data transformation
    int *sendcnts, *sdispls, *recvcnts, *rdispls;
    sendcnts = new int[numProc];
    sdispls = new int[numProc];
    recvcnts = new int[numProc];
    rdispls = new int[numProc];
    Assign_Data(sendcnts, sdispls, recvcnts, rdispls, newLocalSize, localSize, comm);

    //Use alltoall to reallocate to [begin, end);
    MPI_Alltoallv(cache, sendcnts, sdispls, MPI_INTEGER, begin, recvcnts, rdispls, MPI_INTEGER, comm);

    delete [] cache;
    delete [] sendcnts;
    delete [] sdispls;
    delete [] recvcnts;
    delete [] rdispls;

    return;
}


/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

// ...
