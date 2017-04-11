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
    
    // Get the processor number & local array size
    int numprocs;
    int local_size = end - begin;
    MPI_Comm_size(comm, &numprocs);

    // Use new list to store [begin, end)
    int* sorted_list = new int[local_size];
    for(int i = 0; i < local_size; i++) {
    	sorted_list[i] = *(begin + i);
    }
    std::srand(time(NULL));

    // Sort 'list', get new local size.
    int newlocal_size = Quick_sorting(sorted_list, local_size, comm);

    // Compute the configure for data transformation
    int* sendcnts = new int[numprocs];
    int* senddis = new int[numprocs];
    int* recvcnts = new int[numprocs];
    int* recvdis = new int[numprocs];
    Distribute_data(sendcnts, senddis, recvcnts, recvdis, newlocal_size, local_size, comm);

    // Use alltoall to reallocate to [begin, end);
    MPI_Alltoallv(sorted_list, sendcnts, senddis, MPI_INTEGER, begin, recvcnts, recvdis, MPI_INTEGER, comm);

    delete [] sorted_list;
    delete [] sendcnts;
    delete [] senddis;
    delete [] recvcnts;
    delete [] recvdis;
    return;
}


/*********************************************************************
 *             Implement your own helper functions here:             *
 *********************************************************************/

int Quick_sorting(int* &begin, int arraysize, MPI_Comm comm) {

    int myRank;
    int numprocs;
    int total_Arraysize;
    int local_Arraysize = arraysize;
    MPI_Comm_rank(comm, &myRank);
    MPI_Comm_size(comm, &numprocs);


    // Compute the size of whole array
    MPI_Allreduce(&local_Arraysize, &total_Arraysize, 1, MPI_INTEGER, MPI_SUM, comm);

    // Case 1: Processor > elements
    // Construct new comm with 0,1,..,n-1;
    if (numprocs > total_Arraysize) {
        MPI_Comm newComm;
        if (myRank < total_Arraysize) {
            MPI_Comm_split(comm, 0, myRank, &newComm);
            return Quick_sorting(begin, arraysize, newComm);
        }
        else {
          return 0;
        }
    }

    // Case 2: Processor == 1
    if (numprocs == 1) {
        std::sort(begin, begin + arraysize);
        return arraysize;
    }

    // Case 3: 1< Processor <= n

    // Randomly choose the pivot_element [0, total_Arraysize)
    int pivot_element, pivot_rank;

    // Randomly generate the processor rank where the pivot_element locates.
    if (myRank == 0) {
        pivot_rank = rand() % numprocs;
    }
    MPI_Bcast(&pivot_rank, 1, MPI_INTEGER, 0, comm);

    if (myRank == pivot_rank) {
        pivot_element = *(begin + pivot_rank);
    }
    MPI_Bcast(&pivot_element, 1, MPI_INTEGER, pivot_rank, comm);

    // Partition locally array
    // Count the number of left elements on pivot_element;
    int numLeft;
    numLeft = Partition_Array(begin, begin + arraysize, pivot_element);
    int numTotalLeft;
    MPI_Allreduce(&numLeft, &numTotalLeft, 1, MPI_INTEGER, MPI_SUM, comm);


    // Calculate the distination of sending & receiving
    // For greater
    int* sendcnts = new int[numprocs];
    int* senddis = new int[numprocs];
    int* recvcnts = new int[numprocs];
    int* recvdis = new int[numprocs];
    // For smaller
    int* sendcnts_l = new int[numprocs];
    int* senddis_l = new int[numprocs];
    int* recvcnts_l = new int[numprocs];
    int* recvdis_l = new int[numprocs];

    //Partition the processor to the subproblem
    int numprocsLeft, numprocsRight;
    Partition_Processor(numprocsLeft, numprocsRight, numprocs, numLeft, local_Arraysize - numLeft, comm);


    // Reassigned local array size
    int new_local_size; 
    if (myRank < numprocsLeft) {
  		new_local_size = block_decompose(numTotalLeft, numprocsLeft, myRank);
    }
    else {
  		new_local_size = block_decompose(total_Arraysize - numTotalLeft, numprocsRight, myRank - numprocsLeft);
    }

    // For data transform
    if(myRank < numprocsLeft) {
        // For the smaller or equal to elements
        Distribute_data(sendcnts_l, senddis_l, recvcnts_l, recvdis_l, numLeft, new_local_size, comm);
        // For the larger elements
        Distribute_data(sendcnts, senddis, recvcnts, recvdis, local_Arraysize - numLeft, 0, comm);
    }
    else {
        Distribute_data(sendcnts_l, senddis_l, recvcnts_l, recvdis_l, numLeft, 0, comm);
        Distribute_data(sendcnts, senddis, recvcnts, recvdis, local_Arraysize - numLeft, new_local_size, comm);
    }

    if(myRank < numprocsLeft) {
        for(int i = 0; i < numprocsLeft; i++) {
            sendcnts[i] = sendcnts_l[i];
            recvcnts[i] = recvcnts_l[i];

        }
        for(int i = numprocsLeft; i < numprocs; i++) {
            recvcnts[i] = recvcnts_l[i];
        }
    }
    else {
        for(int i = 0; i < numprocsLeft; i++) {
      		sendcnts[i] = sendcnts_l[i];
        }
    }
    delete [] sendcnts_l;
    delete [] senddis_l;
    delete [] recvcnts_l;
    delete [] recvdis_l;

    senddis[0] = 0; recvdis[0] = 0;
    for(int i = 1; i < numprocs; i++) {
        senddis[i] = senddis[i - 1] + sendcnts[i - 1];
        recvdis[i] = recvdis[i - 1] + recvcnts[i - 1];
    }

    //Exchange data
    int* temp = new int[new_local_size];
    MPI_Alltoallv(begin, sendcnts, senddis, MPI_INTEGER, temp, recvcnts, recvdis, MPI_INTEGER, comm);

    delete [] begin;
    begin = new int[new_local_size];
    for(int i = 0; i < new_local_size; i++) {
  		begin[i] = temp[i];
    }
    delete [] temp;
    // Construct new communicator
    MPI_Comm newcomm;
    MPI_Comm_split(comm, (myRank < numprocsLeft), myRank, &newcomm);

    int ret_size = Quick_sorting(begin, new_local_size, newcomm);


    delete [] sendcnts;
    delete [] senddis;
    delete [] recvcnts;
    delete [] recvdis;
    MPI_Comm_free(&newcomm);

    return ret_size;
}

int Partition_Array(int* begin, int* end, int pivot_element) {
    int cnt = 0;
    int tmp;
    for(int i = 0; i < end - begin; i++) {
        if(begin[i] <= pivot_element) {
        	tmp = begin[cnt];
            begin[cnt] = begin[i];
            begin[i] = tmp;
            cnt++;
        }
    }
    //numLeft = cnt;
    return cnt;
}

void Partition_Processor(int &procLeft, int &procRight, int numprocs, int numLeft, int numRight, MPI_Comm comm) {
    int sum_left, sum_right;
    MPI_Allreduce(&numLeft, &sum_left, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allreduce(&numRight, &sum_right, 1, MPI_INTEGER, MPI_SUM, comm);

    if(sum_left < sum_right) {
  		procLeft = (numprocs * sum_left)/(sum_left + sum_right);
    }
    else {
  		procLeft = (numprocs * sum_left - 1)/(sum_left + sum_right) + 1;
    }
    procRight = numprocs - procLeft;

    if(procLeft == 0 && sum_left > 0) {
        procLeft = 1;
        procRight--;
    }

    if(procRight == 0 && sum_right > 0) {
        procRight = 1;
        procLeft--;
    }
    return;
}



// if(myRank < numprocsLeft) {
//     // For the smaller or equal to elements
//     Distribute_data(tmp_sendcnts, tmp_senddis, tmp_recvcnts, tmp_recvdis, numLeft, new_local_size, comm);
//     // For the larger elements
//     Distribute_data(sendcnts, senddis, recvcnts, recvdis, local_Arraysize - numLeft, 0, comm);
// }
// else {
//     Distribute_data(tmp_sendcnts, tmp_senddis, tmp_recvcnts, tmp_recvdis, numLeft, 0, comm);
//     Distribute_data(sendcnts, senddis, recvcnts, recvdis, local_Arraysize - numLeft, new_local_size, comm);
// }




void Distribute_data(int* sendcnts, int* senddis, int* recvcnts, int* recvdis, int curlen, int objlen, MPI_Comm &comm) {
    int myRank;
    int numprocs;
    MPI_Comm_size(comm, &numprocs);
    MPI_Comm_rank(comm, &myRank);

    int* pre_curlen = new int[numprocs];
    int* pre_objlen = new int[numprocs];
    int tmp_cur, tmp_obj;

    //Collect the prefix sum of current size and allgather them
    MPI_Scan(&curlen, &tmp_cur, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allgather(&tmp_cur, 1, MPI_INTEGER, pre_curlen, 1, MPI_INTEGER, comm);

    //Collect the prefix sum of objective size and allgather them
    MPI_Scan(&objlen, &tmp_obj, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allgather(&tmp_obj, 1, MPI_INTEGER, pre_objlen, 1, MPI_INTEGER, comm);

    int head = 0, tail = 0;

    //Find the index of processors [head, tail] to RECEIVE data from
    if (myRank != 0) {
    	while(pre_curlen[head] < pre_objlen[myRank - 1]) {
    		head++;
   		}
    }

    tail = head;
    while(pre_curlen[tail] < pre_objlen[myRank]) {
    	tail++;
   	}

    // Setup the RECEIVE counts
    for(int i = 0; i < head; i++) {
  		recvcnts[i] = 0;
    }

    if(tail == head) {
  		recvcnts[head] = objlen;
    }
    else {
        if (myRank > 0) {
            recvcnts[head] = pre_curlen[head] - pre_objlen[myRank - 1];
        }
        else {
            recvcnts[head] = pre_curlen[head];
        }
        //recvcnts[head] = pre_curlen[head] - ((myRank > 0) ? pre_objlen[myRank - 1] : 0);
        for(int i = head + 1; i < tail; i++) {
  			recvcnts[i] = pre_curlen[i] - pre_curlen[i - 1];
        }
        recvcnts[tail] = pre_objlen[myRank] - pre_curlen[tail - 1];
    }
    for(int i = tail + 1; i < numprocs; i++) {
  		recvcnts[i] = 0;
    }

    // For SEND data
    head = 0; tail = 0;
    if (myRank != 0) {
  		while(pre_objlen[head] < pre_curlen[myRank - 1]) {
  			head++;
  		}
    }

    tail = head;
    while(pre_objlen[tail] < pre_curlen[myRank]) {
    	tail++;
    }

    // Setup the SEND counts
    for(int i = 0; i < head; i++) {
  		sendcnts[i] = 0;
    }

    if(tail == head) {
  		sendcnts[head] = curlen;
    }
    else {
        if (myRank > 0) {
            sendcnts[head] = pre_objlen[head] - pre_curlen[myRank - 1];
        }
        else {
            sendcnts[head] = pre_objlen[head];
        }
        for(int i = head + 1; i < tail; i++) {
      		sendcnts[i] = pre_objlen[i] - pre_objlen[i - 1];
        }
        sendcnts[tail] = pre_curlen[myRank] - pre_objlen[tail - 1];
    }
    for(int i = tail + 1; i < numprocs; i++) {
  		sendcnts[i] = 0; 
    }

    // Setup Displacement
    senddis[0] = 0; 
    recvdis[0] = 0;
    for(int i = 1; i < numprocs; i++) {
        senddis[i] = senddis[i - 1] + sendcnts[i - 1];
        recvdis[i] = recvdis[i - 1] + recvcnts[i - 1];
    }
    delete [] pre_curlen;
    delete [] pre_objlen;
    return;
}





// ...
