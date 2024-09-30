#pragma once
#include"kernel_functions.cuh"
#include<iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>

#define CUDA_CHECK_ERROR() { \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA error: " << cudaGetErrorString(err) << " at line " << __LINE__ << std::endl; \
        exit(EXIT_FAILURE); \
    } \
}

using namespace std;

namespace tools {
    __device__ 
    void bobbleSortPairs(Graph_D_H::TDrank* array, int left, int right)
    {
        for (int i = left; i <= right; i++) {
            //bool swapped = false;

            for (int j = left; j <= right - 1; j++) {
                if (array[j + 1] < array[j]) {
                    Graph_D_H::TDrank temp = array[j];
                    array[j] = array[j + 1];
                    array[j + 1] = temp;

                    //swapped = true;
                }
            }

            //if (!swapped) {
            //    break;
            //}
        }
    }

    __device__
        void bobbleSortPairs_2(Graph_D_H::TDrank* array, int left, int right)
    {
        for (int i = left; i <= right; i++) {
            //bool swapped = false;

            for (int j = left; j <= right - 1; j++) {
                if (array[j + 1] < array[j]) {
                    Graph_D_H::TDrank temp = array[j];
                    array[j] = array[j + 1];
                    array[j + 1] = temp;

                    //swapped = true;
                }
            }

            //if (!swapped) {
               // break;
            //}
        }
    }

    __device__
        void bobbleSortmyPair(Graph_D_H::myPair<int>* array, int left, int right)
    {
        for (int i = left; i <= right; i++) {
            bool swapped = false;

            for (int j = left; j <= right - 1; j++) {
                if (array[j + 1] < array[j]) {
                    Graph_D_H::myPair<int> temp = array[j];
                    array[j] = array[j + 1];
                    array[j + 1] = temp;

                    swapped = true;
                }
            }

            if (!swapped) {
                break;
            }
        }
    }

    __device__ int partition(Graph_D_H::TDrank* array, int low, int high) {
        Graph_D_H::TDrank pivot = array[high];
        int i = (low - 1);

        for (int j = low; j <= high - 1; j++) {
            if (array[j] <  pivot) {
                i++;
                Graph_D_H::TDrank temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }

        Graph_D_H::TDrank temp = array[i + 1];
        array[i + 1] = array[high];
        array[high] = temp;
        return (i + 1);
    }

    __device__ void quickSort(Graph_D_H::TDrank* array, int low, int high) {
        if (low < high) {
            int pi = partition(array, low, high);
            quickSort(array, low, pi - 1);
            quickSort(array, pi + 1, high);
        }
    }

}

__global__
void sets(int s)
{
    //int index = threadIdx.x;
    printf("cuda kernel in\n");
}

void tests()
{
    sets<<<1,1>>>(1);
    cudaDeviceSynchronize();
}

void checkCudaError(cudaError_t err) {
    if (err != cudaSuccess) {
        std::cerr << "CUDA Error: " << cudaGetErrorString(err) << std::endl;
        exit(err);
    }
}


__global__
void makeCHPerHeight_D(int lowestIndexStart, int lowestIndexEnd,
    Graph_D_H::TDrank* nonAdjcentNode,
    int* ID_hash, 
    Graph_D_H::myPair<int>* CHTreeHash, 
    Graph_D_H::myPair<int>* NANHash,
    Graph_D_H::myPair<int>* CHTree, 
    bool* visited)
{
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tid > lowestIndexEnd - lowestIndexStart)
        return;
    int nonAdjStartIndex = NANHash[tid + lowestIndexStart].first;
    const int TDNodeSize = NANHash[tid + lowestIndexStart].NodeID;

    for (int j = nonAdjStartIndex; j < nonAdjStartIndex + TDNodeSize; j++)
    {
        //tools::quickSort(nonAdjcentNode, j, nonAdjStartIndex + TDNodeSize -1);
        tools::bobbleSortPairs(nonAdjcentNode, j, nonAdjStartIndex + TDNodeSize - 1);
        //for (int it = j; it < nonAdjStartIndex + TDNodeSize; it++) {

        //}
        int TDNodeNow = nonAdjcentNode[j].second;
        visited[TDNodeNow] = true;
        int CHTreeHash_index = ID_hash[TDNodeNow];
        int CHTree_Index = CHTreeHash[CHTreeHash_index].NodeID;
        int CHTree_ActualSize = CHTreeHash[CHTreeHash_index].first;
        //int act_empty = firstempty[CHTreeHash_index];
        int act = 1;
        tools::bobbleSortmyPair(CHTree,CHTree_Index,CHTree_Index + CHTree_ActualSize - 1);
        for (int k = CHTree_Index + 1; k < CHTree_Index + CHTree_ActualSize; k++)
        {
            act++;
            if (CHTree[k].first == CHTree[k - 1].first && CHTree[k].first != INT_MAX)
            {
                if (CHTree[k - 1].second < CHTree[k].second)
                {
                    CHTree[k].second = CHTree[k - 1].second;
                    CHTree[k].NodeID = CHTree[k - 1].NodeID;
                }
                CHTree[k - 1].first = INT_MAX;
                act--;
            }
           // CHTree[k].first == INT_MAX;
        }
        tools::bobbleSortmyPair(CHTree, CHTree_Index, CHTree_Index + CHTree_ActualSize - 1);
        for (int k = CHTree_Index; k < CHTree_Index + CHTree_ActualSize; k++)
        {
            if (CHTree[k].first == INT_MAX)
            {
                break;
            }
            CHTreeHash[CHTreeHash_index].first = k - CHTree_Index + 1;
        }
        CHTree_ActualSize = CHTreeHash[CHTreeHash_index].first;

        for (int k = CHTree_Index; k < CHTree_Index + CHTree_ActualSize; k++)
        {
            int outpoint = CHTree[k].first;
            int outpointCHTreeHash_index = ID_hash[outpoint];
            int outPoint_CHTree_index = CHTreeHash[outpointCHTreeHash_index].NodeID;
            int outPoint_CHTree_actualSize = CHTreeHash[outpointCHTreeHash_index].first;
            for (int q = CHTreeHash[outpointCHTreeHash_index].NodeID; q < CHTreeHash[outpointCHTreeHash_index].NodeID + CHTreeHash[outpointCHTreeHash_index].first; q++)
            {
                if (CHTree[q].first == TDNodeNow)
                {
                    CHTree[q].first = INT_MAX;
                    //tools::bobbleSortmyPair(CHTree, outPoint_CHTree_index, outPoint_CHTree_index + outPoint_CHTree_actualSize - 1);
                    //CHTreeHash[outpointCHTreeHash_index].first--;
                    //break;
                }
            }
            tools::bobbleSortmyPair(CHTree, outPoint_CHTree_index, outPoint_CHTree_index + outPoint_CHTree_actualSize - 1);
            for (int k = CHTreeHash[outpointCHTreeHash_index].NodeID; k < outPoint_CHTree_index + outPoint_CHTree_actualSize; k++)
            {
                if (CHTree[k].first == INT_MAX)
                {
                    break;
                }
                CHTreeHash[outpointCHTreeHash_index].first = k - CHTreeHash[outpointCHTreeHash_index].NodeID + 1;
            }
        }
        //add left->right and right -> left
        for (int k = CHTree_Index; k < CHTree_Index + CHTree_ActualSize; k++)
        {
            //if (CHTree[k].first == INT_MAX)
               // continue;
            int leftID = CHTree[k].first;
            for (int q = k + 1; q < CHTree_Index + CHTree_ActualSize; q++)
            {
               // if (CHTree[q].first == INT_MAX)
                   // continue;
                int rightID = CHTree[q].first;
                if (leftID == rightID)
                    continue;
                int length = CHTree[k].second + CHTree[q].second;
                bool leftHaveRight = false;
               
                for (int p = CHTreeHash[ID_hash[leftID]].NodeID; p < CHTreeHash[ID_hash[leftID]].NodeID + CHTreeHash[ID_hash[leftID]].first; p++)
                {
                    if (CHTree[p].first == rightID)
                    {
                        leftHaveRight = true;
                        if (CHTree[p].second > length)
                        {
                            CHTree[p].second = length;
                            CHTree[p].NodeID = TDNodeNow;
                        }
                        break;
                    }
                }
                for (int p = CHTreeHash[ID_hash[rightID]].NodeID; p < CHTreeHash[ID_hash[rightID]].NodeID + CHTreeHash[ID_hash[rightID]].first; p++)
                {
                    if (CHTree[p].first == leftID)
                    {
                        if (CHTree[p].second > length)
                        {
                            CHTree[p].second = length;
                            CHTree[p].NodeID = TDNodeNow;
                        }
                        break;
                    }
                }
                if (leftHaveRight)
                    continue;
                else
                {
                    //left->right

                    //no such link
                    //left->right
                    //int leftStartIndex = CHTreeHash[ID_hash[leftID]].NodeID;
                    int emptyPlace = CHTreeHash[ID_hash[leftID]].NodeID + CHTreeHash[ID_hash[leftID]].first;
                    CHTree[emptyPlace].NodeID = TDNodeNow;
                    CHTree[emptyPlace].first = rightID;
                    CHTree[emptyPlace].second = length;
                    CHTreeHash[ID_hash[leftID]].first++;

                    //right->left
                    //leftStartIndex = CHTreeHash[ID_hash[rightID]].NodeID;
                    emptyPlace = CHTreeHash[ID_hash[rightID]].NodeID + CHTreeHash[ID_hash[rightID]].first;
                    CHTree[emptyPlace].NodeID = TDNodeNow;
                    CHTree[emptyPlace].first = leftID;
                    CHTree[emptyPlace].second = length;
                    CHTreeHash[ID_hash[rightID]].first++;
                }
                
            }
        }


        //change vertex degree
        for (int k = j + 1; k < nonAdjStartIndex + TDNodeSize; k++)
        {

            nonAdjcentNode[k].first = CHTreeHash[ID_hash[nonAdjcentNode[k].second]].first;
            for (int q = CHTree_Index; q < CHTree_Index + CHTree_ActualSize; q++)
            {
                if (CHTree[q].first == INT_MAX)
                    continue;
                //int outpoint = CHTree[q].first;
                if (nonAdjcentNode[k].second == CHTree[q].first)
                {
                    nonAdjcentNode[k].third = (nonAdjcentNode[j].third + 1 > nonAdjcentNode[k].third) ? nonAdjcentNode[j].third + 1 : nonAdjcentNode[k].third;
                    break;
                }
            }
        }
    }
}


__global__
void makeCHPerHeight_D_2(int lowestIndexStart, int lowestIndexEnd, 
    Graph_D_H::TDrank* nonAdjcentNode,
    int* ID_hash, 
    Graph_D_H::myPair<int>* CHTreeHash, 
    Graph_D_H::myPair<int>* NANHash,
    Graph_D_H::myPair<int>* CHTree, 
    bool* visited, 
    Graph_D_H::pairs* CHInfo)
{
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tid > lowestIndexEnd - lowestIndexStart)
        return;
    int nonAdjStartIndex = NANHash[tid + lowestIndexStart].first;
    const int TDNodeSize = NANHash[tid + lowestIndexStart].NodeID;

    for (int j = nonAdjStartIndex; j < nonAdjStartIndex + TDNodeSize; j++)
    {
        //tools::quickSort(nonAdjcentNode, j, nonAdjStartIndex + TDNodeSize -1);
        tools::bobbleSortPairs(nonAdjcentNode, j, nonAdjStartIndex + TDNodeSize - 1);
        int TDNodeNow = nonAdjcentNode[j].second;

        visited[TDNodeNow] = true;
        int CHTreeHash_index = ID_hash[TDNodeNow];
        int CHTree_Index = CHTreeHash[CHTreeHash_index].NodeID;
        int CHTree_ActualSize = CHTreeHash[CHTreeHash_index].first;
        //int act_empty = firstempty[CHTreeHash_index];

        //reorder
        tools::bobbleSortmyPair(CHTree, CHTree_Index, CHTree_Index + CHTree_ActualSize - 1);
        for (int k = CHTree_Index + 1; k < CHTree_Index + CHTree_ActualSize; k++)
        {

            if (CHTree[k].first == CHTree[k - 1].first && CHTree[k].first != INT_MAX)
            {
                if (CHTree[k - 1].second < CHTree[k].second)
                {
                    CHTree[k].second = CHTree[k - 1].second;
                    CHTree[k].NodeID = CHTree[k - 1].NodeID;
                }
                CHTree[k - 1].first = INT_MAX;
            }
        }
        tools::bobbleSortmyPair(CHTree, CHTree_Index, CHTree_Index + CHTree_ActualSize - 1);
        for (int k = CHTree_Index; k < CHTree_Index + CHTree_ActualSize; k++)
        {
            if (CHTree[k].first == INT_MAX)
            {
                break;
            }
            CHTreeHash[CHTreeHash_index].first = k - CHTree_Index + 1;
        }
        CHTree_ActualSize = CHTreeHash[CHTreeHash_index].first;

        //delete
        for (int k = CHTree_Index; k < CHTree_Index + CHTree_ActualSize; k++)
        {
            int outpoint = CHTree[k].first;
            int outpointCHTreeHash_index = ID_hash[outpoint];
            int outPoint_CHTree_index = CHTreeHash[outpointCHTreeHash_index].NodeID;
            int outPoint_CHTree_actualSize = CHTreeHash[outpointCHTreeHash_index].first;
            for (int q = CHTreeHash[outpointCHTreeHash_index].NodeID; q < CHTreeHash[outpointCHTreeHash_index].NodeID + CHTreeHash[outpointCHTreeHash_index].first; q++)
            {
                if (CHTree[q].first == TDNodeNow)
                {
                    CHTree[q].first = INT_MAX;
                    //tools::bobbleSortmyPair(CHTree, outPoint_CHTree_index, outPoint_CHTree_index + outPoint_CHTree_actualSize - 1);
                    //CHTreeHash[outpointCHTreeHash_index].first--;
                    //break;
                }
            }
            tools::bobbleSortmyPair(CHTree, outPoint_CHTree_index, outPoint_CHTree_index + outPoint_CHTree_actualSize - 1);
            for (int k = CHTreeHash[outpointCHTreeHash_index].NodeID; k < outPoint_CHTree_index + outPoint_CHTree_actualSize; k++)
            {
                if (CHTree[k].first == INT_MAX)
                {
                    break;
                }
                CHTreeHash[outpointCHTreeHash_index].first = k - CHTreeHash[outpointCHTreeHash_index].NodeID + 1;
            }
            
        }
        //add left->right and right -> left
        for (int k = CHTree_Index; k < CHTree_Index + CHTree_ActualSize; k++)
        {
            //if (CHTree[k].first == INT_MAX)
               // continue;
            int leftID = CHTree[k].first;
            for (int q = k + 1; q < CHTree_Index + CHTree_ActualSize; q++)
            {
                // if (CHTree[q].first == INT_MAX)
                    // continue;
                int rightID = CHTree[q].first;
                if (leftID == rightID)
                    continue;
                int length = CHTree[k].second + CHTree[q].second;
                bool leftHaveRight = false;

                for (int p = CHTreeHash[ID_hash[leftID]].NodeID; p < CHTreeHash[ID_hash[leftID]].NodeID + CHTreeHash[ID_hash[leftID]].first; p++)
                {
                    if (CHTree[p].first == rightID)
                    {
                        leftHaveRight = true;
                        if (CHTree[p].second > length)
                        {
                            CHTree[p].second = length;
                            CHTree[p].NodeID = TDNodeNow;
                        }
                        break;
                    }
                }
                for (int p = CHTreeHash[ID_hash[rightID]].NodeID; p < CHTreeHash[ID_hash[rightID]].NodeID + CHTreeHash[ID_hash[rightID]].first; p++)
                {
                    if (CHTree[p].first == leftID)
                    {
                        if (CHTree[p].second > length)
                        {
                            CHTree[p].second = length;
                            CHTree[p].NodeID = TDNodeNow;
                        }
                        break;
                    }
                }
                if (leftHaveRight)
                    continue;
                else
                {
                    //left->right

                    //no such link
                    //left->right
                    //int leftStartIndex = CHTreeHash[ID_hash[leftID]].NodeID;
                    int emptyPlace = CHTreeHash[ID_hash[leftID]].NodeID + CHTreeHash[ID_hash[leftID]].first;
                    CHTree[emptyPlace].NodeID = TDNodeNow;
                    CHTree[emptyPlace].first = rightID;
                    CHTree[emptyPlace].second = length;
                    CHTreeHash[ID_hash[leftID]].first++;

                    //right->left
                    //leftStartIndex = CHTreeHash[ID_hash[rightID]].NodeID;
                    emptyPlace = CHTreeHash[ID_hash[rightID]].NodeID + CHTreeHash[ID_hash[rightID]].first;
                    CHTree[emptyPlace].NodeID = TDNodeNow;
                    CHTree[emptyPlace].first = leftID;
                    CHTree[emptyPlace].second = length;
                    CHTreeHash[ID_hash[rightID]].first++;
                }

            }
        }


        //change vertex degree

        for (int q = CHTree_Index; q < CHTree_Index + CHTree_ActualSize; q++)
        {
            if (CHTree[q].first == INT_MAX)
                continue;
            int outpoint = CHTree[q].first;
            int HeightNow = nonAdjcentNode[j].third;
            CHInfo[outpoint].first = CHTreeHash[ID_hash[outpoint]].first;// global degree
            CHInfo[outpoint].second = (HeightNow + 1 > CHInfo[outpoint].second) ? HeightNow + 1 : CHInfo[outpoint].second;//global height
            
        }
        for (int k = j + 1; k < nonAdjStartIndex + TDNodeSize; k++)
        {
            nonAdjcentNode[k].first = CHTreeHash[ID_hash[nonAdjcentNode[k].second]].first;//degree
            nonAdjcentNode[k].third = CHInfo[nonAdjcentNode[k].second].second; //height
        }

    }
}

void makeCHPerHeight(int TreeHeight,int cudaThreadNum,int TreeHeightTarget,
    thrust::device_vector<Graph_D_H::TDrank>& nonAdjcentNode_D,
    thrust::device_vector<int>& ID_hash_D, 
    thrust::device_vector<Graph_D_H::myPair<int>>& CHTreeHash_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& CHTree_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& NANHash_D,
    thrust::device_vector<bool>& visited_D)
{
    int tempHeight = TreeHeight;
    Graph_D_H::TDrank* nonAdjcentNode = thrust::raw_pointer_cast(nonAdjcentNode_D.data());
    int* ID_hash = thrust::raw_pointer_cast(ID_hash_D.data());
    Graph_D_H::myPair<int>* CHTreeHash = thrust::raw_pointer_cast(CHTreeHash_D.data());
    Graph_D_H::myPair<int>* CHTree = thrust::raw_pointer_cast(CHTree_D.data());
    Graph_D_H::myPair<int>* NANHash = thrust::raw_pointer_cast(NANHash_D.data());
    bool* visited = thrust::raw_pointer_cast(visited_D.data());
   // Graph_D_H::time_Mine time;
    //cout << "TD start " << endl;
    //time.updateStart();
    while (tempHeight > TreeHeightTarget)
    {
        int lowestIndexStart = (int)(std::pow(2, tempHeight - 1) - 1);
        int lowestIndexEnd = (int)std::pow(2, tempHeight) - 2;
       // cout << "height now : " << tempHeight << endl;
       // Graph_D_H::time_Mine time1;

       // time1.updateStart();
        tempHeight--;
        int blockSize = ((lowestIndexStart + 1) / cudaThreadNum) + 1;

        makeCHPerHeight_D<<<blockSize,cudaThreadNum>>>(lowestIndexStart, lowestIndexEnd, nonAdjcentNode, ID_hash, CHTreeHash, NANHash, CHTree,visited);
        cudaDeviceSynchronize();
       //time1.updateEnd();
       // cout << "time using: " << time1.get_microsecond_duration() << endl;
    }
   // time.updateEnd();
   // cout << " using time: " << time.get_microsecond_duration() << endl;
}


void makeCHPerHeight_2(int TreeHeight, int cudaThreadNum, int TreeHeightTarget,
    thrust::device_vector<Graph_D_H::TDrank>& nonAdjcentNode_D,
    thrust::device_vector<int>& ID_hash_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& CHTreeHash_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& CHTree_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& NANHash_D,
    thrust::device_vector<bool>& visited_D,
    thrust::device_vector<Graph_D_H::pairs>& CHInfo_D)
{
    int tempHeight = TreeHeight;
    Graph_D_H::TDrank* nonAdjcentNode = thrust::raw_pointer_cast(nonAdjcentNode_D.data());
    int* ID_hash = thrust::raw_pointer_cast(ID_hash_D.data());
    Graph_D_H::myPair<int>* CHTreeHash = thrust::raw_pointer_cast(CHTreeHash_D.data());
    Graph_D_H::myPair<int>* CHTree = thrust::raw_pointer_cast(CHTree_D.data());
    Graph_D_H::myPair<int>* NANHash = thrust::raw_pointer_cast(NANHash_D.data());
    bool* visited = thrust::raw_pointer_cast(visited_D.data());
    Graph_D_H::pairs* CHInfo = thrust::raw_pointer_cast(CHInfo_D.data());
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    while (tempHeight > TreeHeightTarget)
    {
        int lowestIndexStart = (int)std::pow(2, tempHeight - 1) - 1;
        int lowestIndexEnd = (int)std::pow(2, tempHeight) - 2;
        tempHeight--;
        int blockSize = ((lowestIndexStart + 1) / cudaThreadNum) + 1;
        cudaEventRecord(start);
        makeCHPerHeight_D_2 <<<blockSize, cudaThreadNum >>> (lowestIndexStart, lowestIndexEnd, nonAdjcentNode, ID_hash, CHTreeHash, NANHash, CHTree, visited, CHInfo);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
    }
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}


__global__
void makeH2HLabel_D(Graph_D_H::pairs* H2H_startIndex,
    Graph_D_H::H2HLabel* H2H_label,
    int* frontier, int maxSize
)
{
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tid >= maxSize)
        return;
    int ID = frontier[tid];

    int H2HSize = H2H_startIndex[ID].second - 1;
    int IDNowH2Hindex = H2H_startIndex[ID].first;

    for (int k = IDNowH2Hindex + H2HSize - 1; k >= IDNowH2Hindex; k--)
    {
        if (!H2H_label[k].isTreeNode) {
            continue;
        }
        
        for (int q = k - 1; q >= IDNowH2Hindex; q--)
        {
            if (H2H_label[q].isTreeNode) {
                continue;
            }
            int templength = H2H_label[k].dis + H2H_label[q - IDNowH2Hindex + H2H_startIndex[H2H_label[k].Node].first].dis;
            if (H2H_label[q].dis > templength) {
                H2H_label[q].dis = templength;
                H2H_label[q].Hub = H2H_label[k].Node;
            }
           
        }
    }
}

void makeH2HLabel(int tempHeight,int cudaThreadNum,
    thrust::device_vector<Graph_D_H::pairs>& H2H_startIndex_D, 
    thrust::device_vector<Graph_D_H::H2HLabel>& H2H_label_D, 
    thrust::host_vector<int>& Childs,
    thrust::host_vector<Graph_D_H::pairs>& ChildHash,
    thrust::host_vector<int>& frontier
    ) {
    Graph_D_H::H2HLabel* H2H_label = thrust::raw_pointer_cast(H2H_label_D.data());
    Graph_D_H::pairs* H2H_startIndex = thrust::raw_pointer_cast(H2H_startIndex_D.data());
    thrust::host_vector<int> frontier_now = frontier;
   
   //frontier_now.clear();
    int height = tempHeight;

    //cout << "height now: " << height << endl;


    while (!frontier_now.empty()) {
        thrust::host_vector<int> frontier_next;
        int size = (int)frontier_now.size();
        thrust::device_vector<int> frontier_D = frontier_now;
        int* front = thrust::raw_pointer_cast(frontier_D.data());
        int blockSize = (int)size / cudaThreadNum + 1;

        makeH2HLabel_D<<<blockSize,cudaThreadNum>>>(H2H_startIndex, H2H_label, front, size);
        cudaDeviceSynchronize();

        //cout << " node construct now: ";
        for (auto& ID : frontier_now) {
            //cout << ID << " ";
            for (int j = ChildHash[ID].second; j < ChildHash[ID].second + ChildHash[ID].first; j++)
            {
                frontier_next.push_back(Childs[j]);
            }
        }
        //cout << endl;
        frontier_now.clear();
        frontier_now = frontier_next;
        frontier_next.clear();
        height++;
    }

}


__global__
void makeH2HLabel_noCommunication_D(Graph_D_H::pairs* H2H_startIndex,
    Graph_D_H::H2HLabel* H2H_label,
    int* TreeBFS,int startIndex,int endIndex
)
{
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tid >= endIndex - startIndex)
        return;
    int ID = TreeBFS[tid + startIndex];

    int H2HSize = H2H_startIndex[ID].second - 1;
    int IDNowH2Hindex = H2H_startIndex[ID].first;

    for (int k = IDNowH2Hindex + H2HSize - 1; k >= IDNowH2Hindex; k--)
    {
        if (!H2H_label[k].isTreeNode) {
            continue; 
        }

        for (int q = k - 1; q >= IDNowH2Hindex; q--)
        {
            //if (H2H_label[q].isTreeNode) {
            //	continue;
            //}
            int templength = H2H_label[k].dis + H2H_label[q - IDNowH2Hindex + H2H_startIndex[H2H_label[k].Node].first].dis;
            if (H2H_label[q].dis > templength) {
                H2H_label[q].dis = templength;
                H2H_label[q].Hub = H2H_label[k].Node;
            }
        }
        for (int q = k + 1; q <= IDNowH2Hindex + H2HSize - 1; q++)
        {
            int templength = H2H_label[k].dis + H2H_label[k - IDNowH2Hindex + H2H_startIndex[H2H_label[q].Node].first].dis;
            if (H2H_label[q].dis > templength) {
                H2H_label[q].dis = templength;
                H2H_label[q].Hub = H2H_label[k].Node;
            }
        }
    }
}

void makeH2HLabel_noCommunication(int tempHeight, int cudaThreadNum,
    thrust::device_vector<Graph_D_H::pairs>& H2H_startIndex_D,
    thrust::device_vector<Graph_D_H::H2HLabel>& H2H_label_D,
    thrust::device_vector<int>& TreeBFS_D,
    thrust::host_vector<int>& TreeHash
) {
    Graph_D_H::H2HLabel* H2H_label = thrust::raw_pointer_cast(H2H_label_D.data());
    Graph_D_H::pairs* H2H_startIndex = thrust::raw_pointer_cast(H2H_startIndex_D.data());
    int* TreeBFS = thrust::raw_pointer_cast(TreeBFS_D.data());
    //int* TreeHash = thrust::raw_pointer_cast(TreeHash_D.data());

    size_t maxHeight = TreeHash.size() -1;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    for (int i = tempHeight; i < maxHeight; i++)
    {
        int startIndex = TreeHash[i];
        int endIndex = TreeHash[i + 1];
        int size = endIndex - startIndex;

        int blockSize = size / cudaThreadNum + 1;
        cudaEventRecord(start);
        makeH2HLabel_noCommunication_D <<<blockSize, cudaThreadNum >>>(H2H_startIndex, H2H_label, TreeBFS, startIndex, endIndex);
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
    }
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

}








__global__ void makeH2HLabel_noCommunication_D_noHub_3(
    int64_t* H2H_dis_hash, int* H2H_dis,
    int* TreeBFS_ID, int* TreeBFS_adj, int* TreeBFS_pos, //int* TreeBFS_changeTime,
    int* father,
    int64_t startIndex, int64_t endIndex
) {
    int64_t tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tid >= endIndex - startIndex)
        return;

    int nodeID = TreeBFS_ID[tid + startIndex];
    //int adjID = TreeBFS_adj[tid + startIndex];
    int pos = TreeBFS_pos[tid + startIndex];

    int64_t indexNode = H2H_dis_hash[nodeID];
    int64_t index1_adj = H2H_dis_hash[TreeBFS_adj[tid + startIndex]];
    int tempLength_ = H2H_dis[indexNode + pos];
    for (int64_t j = pos - 1; j >= 0; j--) {
        atomicMin(&H2H_dis[indexNode + j], (tempLength_ + H2H_dis[index1_adj + j]));
    }

    //L4 label
    //int fatherID = nodeID;
    //for (int64_t j = (int)(H2H_dis_hash[nodeID + 1] - H2H_dis_hash[nodeID]) - 2; j > pos; j--) {
    //    fatherID = father[fatherID];
    //    atomicMin(&H2H_dis[indexNode + j], (tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos]));
    //}
}



double makeH2HLabel_noCommunication_noHub_3(int tempHeight, int cudaThreadNum,
    thrust::device_vector<int64_t>& H2H_dis_hash_D,
    thrust::device_vector<int>& H2H_dis_D,
    thrust::device_vector<int>& TreeBFS_ID_D,
    thrust::device_vector<int>& TreeBFS_adj_D,
    thrust::device_vector<int>& TreeBFS_pos_D,
    //thrust::device_vector<int>& TreeBFS_changeTime_D,
    thrust::device_vector<int>& father_D,

    thrust::host_vector<int64_t>& TreeBFS_Hash
) {
    //cout << "H2H label size: \n";
    //cout << "\t pos hash: " << H2H_pos_hash_D.size() << " pos id: " << H2H_pos_ID_D.size() << " pos: " << H2H_pos_POS_D.size() << endl;
    //cout << "\t dis hash: " << H2H_dis_hash_D.size() << " dis: " << H2H_dis_D.size()<<endl;
    //cout << "TreeBFS size: " << TreeBFS_D.size()<<endl ;
    int64_t* H2H_dis_hash = thrust::raw_pointer_cast(H2H_dis_hash_D.data());
    int* H2H_dis = thrust::raw_pointer_cast(H2H_dis_D.data());
    int* TreeBFS_ID = thrust::raw_pointer_cast(TreeBFS_ID_D.data());
    int* TreeBFS_adj = thrust::raw_pointer_cast(TreeBFS_adj_D.data());
    int* TreeBFS_pos = thrust::raw_pointer_cast(TreeBFS_pos_D.data());
    //int* TreeBFS_changeTime = thrust::raw_pointer_cast(TreeBFS_changeTime_D.data());
    int* father = thrust::raw_pointer_cast(father_D.data());
    cudaEvent_t start, stop;
    double Usingtime = 0;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    size_t maxHeight = TreeBFS_Hash.size() - 1;

    for (int i = tempHeight; i < maxHeight; i++)
    {
        int64_t startIndex = TreeBFS_Hash[i];
        int64_t endIndex = TreeBFS_Hash[i + 1];
        int size = endIndex - startIndex;
        int blockSize = size / cudaThreadNum + 1;
        cudaEventRecord(start);

        makeH2HLabel_noCommunication_D_noHub_3 << <blockSize, cudaThreadNum >> > (H2H_dis_hash, H2H_dis, TreeBFS_ID, TreeBFS_adj, TreeBFS_pos, father, startIndex, endIndex);

        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        Usingtime += milliseconds * 1000;
        //std::cout << "At height: " << i << " Kernel execution time: " << milliseconds << " ms" << std::endl;
    }
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    return Usingtime;
}






__global__ void makeH2HLabel_noCommunication_D_noHub(
    int64_t* H2H_pos_hash, int64_t* H2H_dis_hash, int* H2H_pos_ID, int* H2H_pos_POS, int* H2H_dis, int* father, 
    int* TreeBFS, int startIndex, int endIndex
)
{
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tid >= endIndex - startIndex)
        return;
    //printf(" nodeID: %d\n" , tid);
    int nodeID = TreeBFS[tid + startIndex];
    //if (nodeID == 2) {
       // printf("\t %d size is: %d\n", nodeID, H2H_pos_hash[nodeID + 1] - H2H_pos_hash[nodeID]);
    //}
    for (int64_t temp = H2H_pos_hash[nodeID + 1] - 2; temp >= H2H_pos_hash[nodeID]; temp--) {
        int adjID = H2H_pos_ID[temp];
        int pos = H2H_pos_POS[temp];
        int64_t indexNode = H2H_dis_hash[nodeID];
        int64_t index1_adj = H2H_dis_hash[adjID];
        int tempLength_ = H2H_dis[indexNode + pos];
        for (int64_t j = pos - 1; j >= 0; j--) {
            bool check = H2H_dis[indexNode + j] > tempLength_ + H2H_dis[index1_adj + j];
            H2H_dis[indexNode + j] = check * (tempLength_ + H2H_dis[index1_adj + j]) + (!check) * (H2H_dis[indexNode + j]);
            //if (H2H_dis[indexNode + j] > tempLength_ + H2H_dis[index1_adj + j]) {
            //    H2H_dis[indexNode + j] = tempLength_ + H2H_dis[index1_adj + j];
            //}
        }

        //L4 label
        //int fatherID = nodeID;
        //for (int64_t j = (int)(H2H_dis_hash[nodeID + 1] - H2H_dis_hash[nodeID]) - 2; j > pos; j--) {
        //    fatherID = father[fatherID];
        //    bool check = H2H_dis[indexNode + j] > tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos];
        //    H2H_dis[indexNode + j] = check * (tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos]) + (!check) * (H2H_dis[indexNode + j]);

        //    //if (H2H_dis[indexNode + j] > tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos]) {
        //    //    H2H_dis[indexNode + j] = tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos];
        //    //}
        //}
    }
}



double makeH2HLabel_noCommunication_noHub(int tempHeight, int cudaThreadNum,
    thrust::device_vector<int64_t>& H2H_pos_hash_D,
    thrust::device_vector<int64_t>& H2H_dis_hash_D,
    thrust::device_vector<int>& H2H_pos_ID_D,
    thrust::device_vector<int>& H2H_pos_POS_D,
    thrust::device_vector<int>& H2H_dis_D ,
    thrust::device_vector<int>& father_D,
    thrust::device_vector<int>& TreeBFS_D,
    thrust::host_vector<int>& TreeHash
) {
    //cout << "H2H label size: \n";
    //cout << "\t pos hash: " << H2H_pos_hash_D.size() << " pos id: " << H2H_pos_ID_D.size() << " pos: " << H2H_pos_POS_D.size() << endl;
    //cout << "\t dis hash: " << H2H_dis_hash_D.size() << " dis: " << H2H_dis_D.size()<<endl;
    //cout << "TreeBFS size: " << TreeBFS_D.size()<<endl ;

    int64_t* H2H_pos_hash = thrust::raw_pointer_cast(H2H_pos_hash_D.data());
    int64_t* H2H_dis_hash = thrust::raw_pointer_cast(H2H_dis_hash_D.data());
    int* H2H_pos_ID = thrust::raw_pointer_cast(H2H_pos_ID_D.data());
    int* H2H_pos_POS = thrust::raw_pointer_cast(H2H_pos_POS_D.data());
    int* H2H_dis = thrust::raw_pointer_cast(H2H_dis_D.data());
    int* father = thrust::raw_pointer_cast(father_D.data());
    int* TreeBFS = thrust::raw_pointer_cast(TreeBFS_D.data());
    cudaEvent_t start, stop;
    double Usingtime = 0;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    size_t maxHeight = TreeHash.size() - 1;
    for (int i = tempHeight; i < maxHeight; i++)
    {
        int startIndex = TreeHash[i];
        int endIndex = TreeHash[i + 1];
        int size = endIndex - startIndex;
        //cout << "at Height: " << i << " and size: " << size << endl;
        int blockSize = size / cudaThreadNum + 1;
        cudaEventRecord(start);
        makeH2HLabel_noCommunication_D_noHub << <blockSize, cudaThreadNum >> > (H2H_pos_hash, H2H_dis_hash, H2H_pos_ID, H2H_pos_POS, H2H_dis, father, TreeBFS, startIndex, endIndex);
        //cudaDeviceSynchronize();
        cudaEventRecord(stop);
        cudaEventSynchronize(stop);
        float milliseconds = 0;
        cudaEventElapsedTime(&milliseconds, start, stop);
        Usingtime += milliseconds*1000;
        //std::cout <<"At height: "<< i << " Kernel execution time: " << milliseconds << " ms" << std::endl;
    }
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    return Usingtime;
}


__global__ void makeH2HLabel_noCommunication_D_noHub_2(
    int64_t* H2H_pos_hash, int64_t* H2H_dis_hash, int* H2H_pos_ID, int* H2H_pos_POS, int* H2H_dis, int* father,
    int* frontier_Hash, int* frontier_ID, int size
)
{
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (tid >= size)
        return;
    int start = frontier_Hash[tid];
    int end = frontier_Hash[tid + 1];
    for (int i = start; i < end; i++) {

        int nodeID = frontier_ID[i];
        for (int64_t temp = H2H_pos_hash[nodeID + 1] - 2; temp >= H2H_pos_hash[nodeID]; temp--) {
            int adjID = H2H_pos_ID[temp];
            int pos = H2H_pos_POS[temp];
            int64_t indexNode = H2H_dis_hash[nodeID];
            int64_t index1_adj = H2H_dis_hash[adjID];
            int tempLength_ = H2H_dis[indexNode + pos];
            for (int j = pos - 1; j >= 0; j--) {
                if (H2H_dis[indexNode + j] > tempLength_ + H2H_dis[index1_adj + j]) {
                    H2H_dis[indexNode + j] = tempLength_ + H2H_dis[index1_adj + j];
                }
            }

            //L4 label
            //int fatherID = nodeID;
            //for (int j = (int)(H2H_dis_hash[nodeID + 1] - H2H_dis_hash[nodeID]) - 2; j > pos; j--) {
            //    fatherID = father[fatherID];
            //    if (H2H_dis[indexNode + j] > tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos]) {
            //        H2H_dis[indexNode + j] = tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos];
            //    }
            //}
        }

    }
}

void makeH2HLabel_noCommunication_noHub_2(int cudaThreadNum,
    thrust::device_vector<int64_t>& H2H_pos_hash_D,
    thrust::device_vector<int64_t>& H2H_dis_hash_D,
    thrust::device_vector<int>& H2H_pos_ID_D,
    thrust::device_vector<int>& H2H_pos_POS_D,
    thrust::device_vector<int>& H2H_dis_D,
    thrust::device_vector<int>& father_D,

    thrust::device_vector<int> frontier_Hash_D,
    thrust::device_vector<int> frontier_ID_D,
    int64_t size
) {
    //cout << "H2H label size: \n";
    //cout << "\t pos hash: " << H2H_pos_hash_D.size() << " pos id: " << H2H_pos_ID_D.size() << " pos: " << H2H_pos_POS_D.size() << endl;
    //cout << "\t dis hash: " << H2H_dis_hash_D.size() << " dis: " << H2H_dis_D.size()<<endl;
    //cout << "TreeBFS size: " << TreeBFS_D.size()<<endl ;
    int64_t* H2H_pos_hash = thrust::raw_pointer_cast(H2H_pos_hash_D.data());
    int64_t* H2H_dis_hash = thrust::raw_pointer_cast(H2H_dis_hash_D.data());
    int* H2H_pos_ID = thrust::raw_pointer_cast(H2H_pos_ID_D.data());
    int* H2H_pos_POS = thrust::raw_pointer_cast(H2H_pos_POS_D.data());
    int* H2H_dis = thrust::raw_pointer_cast(H2H_dis_D.data());
    int* father = thrust::raw_pointer_cast(father_D.data());

    int* frontier_Hash = thrust::raw_pointer_cast(frontier_Hash_D.data());
    int* frontier_ID = thrust::raw_pointer_cast(frontier_ID_D.data());


    int blockSize = size / cudaThreadNum + 1;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);
    makeH2HLabel_noCommunication_D_noHub_2 << <blockSize, cudaThreadNum >> > (H2H_pos_hash, H2H_dis_hash, H2H_pos_ID, H2H_pos_POS, H2H_dis, father, frontier_Hash, frontier_ID, size);
    //cudaDeviceSynchronize();
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);

        // \BC\C6\CB\E3ʱ\BC\E4
        float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    std::cout << "Kernel execution time: " << milliseconds << " ms" << std::endl;
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}







__global__
void H2Hquery_UsingLCA_D(int* x, int* y, int* result, const int querySize, Graph_D_H::pairs* H2H_startIndex, 
    Graph_D_H::H2HLabel* H2H_label,
    int* firstAppare, Graph_D_H::pairs* RMQ_OneLine, int RMQ_Size, int RMQ_Line_Size) {
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (querySize <= tid)
        return;
    if (x[tid] == y[tid])
    {
        result[tid] = 0;
        return;
    }
    int FirstA = firstAppare[x[tid]];//position  in  OULA_DFS
    int FirstB = firstAppare[y[tid]];
    int firstA = (FirstA < FirstB) * FirstA + (FirstB <= FirstA) * FirstB;
    int firstB = (FirstA < FirstB) * FirstB + (FirstB <= FirstA) * FirstA;
    int k = (int)(log2((double)(firstB - firstA + 1)));//k means log2(gap);

    bool check = (RMQ_OneLine[firstA * RMQ_Line_Size + k].first <= RMQ_OneLine[(firstB - (1 << k) + 1) * RMQ_Line_Size + k].first);
    int LCA = check * RMQ_OneLine[firstA * RMQ_Line_Size + k].second +
        (!check) * RMQ_OneLine[(firstB - (1 << k) + 1) * RMQ_Line_Size + k].second;

    int checkPointX = H2H_startIndex[x[tid]].first;
    int checkPointY = H2H_startIndex[y[tid]].first;
    int LCApoint = H2H_startIndex[LCA].first;
    int frontier = H2H_startIndex[LCA].second;
    int resu = INT_MAX;
    for (int i = 0; i < frontier; i++)
    {
        if (!H2H_label[LCApoint + i].isTreeNode)
            continue;
        resu = ((H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) < resu) ?
            (H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) : resu;
    }
    result[tid] = resu;

}

long long int H2HQuery_UsingLCA(
    thrust::device_vector<int>& x,
    thrust::device_vector<int>& y,
    thrust::device_vector<int>& result,
    thrust::device_vector<int>& firstAppeare_D,
    thrust::device_vector<Graph_D_H::pairs>& RMQ_OneLine_D,
    thrust::device_vector<Graph_D_H::pairs>& H2H_startIndex_D,
    thrust::device_vector<Graph_D_H::H2HLabel>& H2H_label_D, 
    int RMQ_Size,int  RMQ_Line_Size, int cudaThreadNum,int querysize) {
    int querySize = querysize;
    int* X = thrust::raw_pointer_cast(x.data());
    int* Y = thrust::raw_pointer_cast(y.data());
    int* RESULT = thrust::raw_pointer_cast(result.data());
    int* firstAppare = thrust::raw_pointer_cast(firstAppeare_D.data());
    Graph_D_H::pairs* RMQ_Oneline = thrust::raw_pointer_cast(RMQ_OneLine_D.data());
    Graph_D_H::pairs* H2H_startIndex = thrust::raw_pointer_cast(H2H_startIndex_D.data());
    Graph_D_H::H2HLabel* H2H_label = thrust::raw_pointer_cast(H2H_label_D.data());

    int blockSize = querySize / cudaThreadNum + 1;

    Graph_D_H::time_Mine time;
    cout << "H2H_D_UsingLCA start query,querySize=" <<querysize << endl;
    time.updateStart();

    H2Hquery_UsingLCA_D<<<blockSize, cudaThreadNum >>>(X, Y, RESULT, querySize, H2H_startIndex, H2H_label, firstAppare,
        RMQ_Oneline, RMQ_Size, RMQ_Line_Size);
    cudaDeviceSynchronize();

    time.updateEnd();
    cout << "using time: " << time.get_microsecond_duration() << "us" << endl;
    return time.get_microsecond_duration();
}




__global__
void H2Hquery_NoLCA_D(int* x, int* y, int* result, const int querySize, Graph_D_H::pairs* H2H_startIndex, 
    Graph_D_H::H2HLabel* H2H_label) {
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (querySize <= tid)
        return;
    if (x[tid] == y[tid])
    {
        result[tid] = 0;
        return;
    }
    int checkPointX = H2H_startIndex[x[tid]].first;
    int checkPointY = H2H_startIndex[y[tid]].first;
    int frontier = -1;
    if (H2H_startIndex[x[tid]].second > H2H_startIndex[y[tid]].second) {
        frontier = H2H_startIndex[y[tid]].second;
    }
    else {
        frontier = H2H_startIndex[x[tid]].second;
    }
    int resu = INT_MAX;
    for (int i = 0; i < frontier; i++)
    {
        if (H2H_label[checkPointX + i].Node != H2H_label[checkPointY + i].Node)
            break;
        if ((H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) < resu) {
            resu = (H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis);
        }
        //resu = ((H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) < resu) ?
        //    (H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) : resu;
    }
    result[tid] = resu;
}

long long int H2HQuery_NoLCA(
    thrust::device_vector<int>& x,
    thrust::device_vector<int>& y,
    thrust::device_vector<int>& result,
    thrust::device_vector<Graph_D_H::pairs>& H2H_startIndex_D,
    thrust::device_vector<Graph_D_H::H2HLabel>& H2H_label_D, int cudaThreadNum,int querysize) {
    //int querySize = querysize;
    int* X = thrust::raw_pointer_cast(x.data());
    int* Y = thrust::raw_pointer_cast(y.data());
    int* RESULT = thrust::raw_pointer_cast(result.data());
    Graph_D_H::pairs* H2H_startIndex = thrust::raw_pointer_cast(H2H_startIndex_D.data());
    Graph_D_H::H2HLabel* H2H_label = thrust::raw_pointer_cast(H2H_label_D.data());

    int blockSize = querysize / cudaThreadNum + 1;

    Graph_D_H::time_Mine time;
    cout << "H2H_D_noLCA start query,querySize=" << querysize << endl;
    time.updateStart();

    H2Hquery_NoLCA_D << <blockSize, cudaThreadNum >> > (X, Y, RESULT, querysize, H2H_startIndex, H2H_label);
    cudaDeviceSynchronize();

    time.updateEnd();
    cout << "using time: " << time.get_microsecond_duration() << "us" << endl;
    return time.get_microsecond_duration();
}






















__global__
void H2HQuery_UsingLCA_D_noHub(int* x, int* y, int* result, const int querySize, 

    int64_t* H2H_pos_hash,
    int* H2H_pos_POS,
    int64_t* H2H_dis_hash,
    int* H2H_dis,

    int* firstAppeare,
    int64_t* RMQHash,
    int* RMQ_ID,
    int* H2H_TreeHeight_Hash
) 
{
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (querySize <= tid)
        return;
    if (x[tid] == y[tid])
    {
        result[tid] = 0;
        return;
    }
    int FirstA = firstAppeare[x[tid]];//position  in  OULA_DFS
    int FirstB = firstAppeare[y[tid]];
    int firstA = (FirstA < FirstB) * FirstA + (FirstB <= FirstA) * FirstB;
    int firstB = (FirstA < FirstB) * FirstB + (FirstB <= FirstA) * FirstA;
    int k = (int)(log2((double)(firstB - firstA + 1)));//k means log2(gap);

    int ID_1 = RMQ_ID[RMQHash[firstA] + (int64_t)k];
    int ID_2 = RMQ_ID[RMQHash[firstB - (1 << k) + 1] + (int64_t)k];
    int LCA = ID_2;
    if (H2H_TreeHeight_Hash[ID_1] < H2H_TreeHeight_Hash[ID_2]) {
      LCA = ID_1;
    }

    int resu = INT_MAX;
    for (int64_t i = H2H_pos_hash[LCA]; i < H2H_pos_hash[LCA + 1]; i++) {
       
        int templength = H2H_dis[H2H_dis_hash[x[tid]] + (int64_t)H2H_pos_POS[i]] + 
            H2H_dis[H2H_dis_hash[y[tid]] + (int64_t)H2H_pos_POS[i]];
        if (templength < resu) {
            resu = templength;
        }
        //atomicMin(&resu, templength);
        //atomicMin(&result[tid], templength);
    }
    result[tid] = resu;
}

long long int H2HQuery_UsingLCA_noHub(
    thrust::device_vector<int>& x,
    thrust::device_vector<int>& y,
    thrust::device_vector<int>& result,

    thrust::device_vector<int64_t>& H2H_pos_hash_D,
    thrust::device_vector<int64_t>& H2H_dis_hash_D,
    thrust::device_vector<int>& H2H_pos_POS_D,
    thrust::device_vector<int>& H2H_dis_D,

    thrust::device_vector<int>& firstAppeare_D,
    thrust::device_vector<int64_t>& RMQHash_D,
    thrust::device_vector<int>& RMQ_ID_D, 
    thrust::device_vector<int>& H2H_TreeHeight_Hash_D,

    int cudaThreadNum, int querysize
) {
    //int querySize = querysize;
    thrust::fill(result.begin(), result.end(), INT_MAX);
    int* X = thrust::raw_pointer_cast(x.data());
    int* Y = thrust::raw_pointer_cast(y.data());
    int* RESULT = thrust::raw_pointer_cast(result.data());


    int64_t* H2H_pos_hash = thrust::raw_pointer_cast(H2H_pos_hash_D.data());
    int64_t* H2H_dis_hash = thrust::raw_pointer_cast(H2H_dis_hash_D.data());
    int* H2H_pos_POS = thrust::raw_pointer_cast(H2H_pos_POS_D.data());
    int* H2H_dis = thrust::raw_pointer_cast(H2H_dis_D.data());

    int* firstAppare = thrust::raw_pointer_cast(firstAppeare_D.data());
    int64_t* RMQHash = thrust::raw_pointer_cast(RMQHash_D.data());
    int* RMQ_ID = thrust::raw_pointer_cast(RMQ_ID_D.data());
    int* H2H_TreeHeight_Hash = thrust::raw_pointer_cast(H2H_TreeHeight_Hash_D.data());

    int blockSize = querysize / cudaThreadNum + 1;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    Graph_D_H::time_Mine time;
    cout << "H2H_D_UsingLCA start query,querySize=" << querysize << endl;
    time.updateStart();
    cudaEventRecord(start);
    H2HQuery_UsingLCA_D_noHub << <blockSize, cudaThreadNum >> > (X, Y, RESULT, querysize, H2H_pos_hash, H2H_pos_POS, H2H_dis_hash, H2H_dis, firstAppare, RMQHash, RMQ_ID, H2H_TreeHeight_Hash);
    //cudaDeviceSynchronize();
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    long long int micro = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    micro = (long long int)(milliseconds * 1000);
    //time.updateEnd();

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    //CUDA_CHECK_ERROR();

    cout << "using time: " << micro << "us" << endl;
    return micro;
}


__global__
void H2HQuery_noLCA_D_noHub(int* x, int* y, int* result, const int querySize,
    int64_t* H2H_dis_hash,
    int* H2H_dis,
    int* firstAppeare,
    int64_t* RMQHash,
    int* RMQ_ID,
    int* H2H_TreeHeight_Hash)
{
    unsigned int tid = threadIdx.x + (blockDim.x * blockIdx.x);
    if (querySize <= tid)
        return;
    if (x[tid] == y[tid])
    {
        result[tid] = 0;
        return;
    }
    int FirstA = firstAppeare[x[tid]];//position  in  OULA_DFS
    int FirstB = firstAppeare[y[tid]];
    int firstA = (FirstA < FirstB) * FirstA + (FirstB <= FirstA) * FirstB;
    int firstB = (FirstA < FirstB) * FirstB + (FirstB <= FirstA) * FirstA;
    int k = (int)(log2((double)(firstB - firstA + 1)));//k means log2(gap);

    int ID_1 = RMQ_ID[RMQHash[firstA] + (int64_t)k];
    int ID_2 = RMQ_ID[RMQHash[firstB - (1 << k) + 1] + (int64_t)k];
    int LCASize = H2H_TreeHeight_Hash[ID_2];
    //atomicMin(&LCASize, H2H_TreeHeight_Hash[ID_1]);
    if (H2H_TreeHeight_Hash[ID_1] < H2H_TreeHeight_Hash[ID_2]) {
        LCASize = H2H_TreeHeight_Hash[ID_1];
    }

    int resu = INT_MAX;

    for (int64_t i = 0; i < LCASize; i++) {

        int templength = H2H_dis[H2H_dis_hash[x[tid]] + i] + H2H_dis[H2H_dis_hash[y[tid]] + i];
        if (templength < resu) {
            resu = templength;
        }
        //atomicMin(&resu, templength);
        //atomicMin(&result[tid], templength);
    }
    result[tid] = resu;
}

long long int H2HQuery_noLCA_noHub(
    thrust::device_vector<int>& x,
    thrust::device_vector<int>& y,
    thrust::device_vector<int>& result,

    thrust::device_vector<int64_t>& H2H_dis_hash_D,
    thrust::device_vector<int>& H2H_dis_D,

    thrust::device_vector<int>& firstAppeare_D,
    thrust::device_vector<int64_t>& RMQHash_D,
    thrust::device_vector<int>& RMQ_ID_D,
    thrust::device_vector<int>& H2H_TreeHeight_Hash_D,
    
    int cudaThreadNum, int querysize) {

    thrust::fill(result.begin(), result.end(), INT_MAX);
    //int querySize = querysize;
    int* X = thrust::raw_pointer_cast(x.data());
    int* Y = thrust::raw_pointer_cast(y.data());
    int* RESULT = thrust::raw_pointer_cast(result.data());


    int64_t* H2H_dis_hash = thrust::raw_pointer_cast(H2H_dis_hash_D.data());
    int* H2H_dis = thrust::raw_pointer_cast(H2H_dis_D.data());

    int* firstAppare = thrust::raw_pointer_cast(firstAppeare_D.data());
    int64_t* RMQHash = thrust::raw_pointer_cast(RMQHash_D.data());
    int* RMQ_ID = thrust::raw_pointer_cast(RMQ_ID_D.data());
    int* H2H_TreeHeight_Hash = thrust::raw_pointer_cast(H2H_TreeHeight_Hash_D.data());


    int blockSize = querysize / cudaThreadNum + 1;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //Graph_D_H::time_Mine time;
    cout << "H2H_D_noLCA start query,querySize = " << querysize << endl;
    //time.updateStart();
    cudaEventRecord(start);
    H2HQuery_noLCA_D_noHub <<<blockSize, cudaThreadNum>>>(X, Y, RESULT, querysize, H2H_dis_hash, H2H_dis, firstAppare, RMQHash, RMQ_ID, H2H_TreeHeight_Hash);
    //cudaDeviceSynchronize();
    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    long long int micro = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);
    micro = (long long int)(milliseconds * 1000);
    //time.updateEnd();

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    //CUDA_CHECK_ERROR();

     cout << "using time: " << micro << "us" << endl;
     return micro;
}
