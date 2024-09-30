#pragma once
#include "Graph_D.cuh"

void tests();

void makeCHPerHeight(int TreeHeight, int cudaThreadNum, int TreeHeightTarget,
    thrust::device_vector<Graph_D_H::TDrank>& nonAdjcentNode_D,
    thrust::device_vector<int>& ID_hash_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& CHTreeHash_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& CHTree_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& NANHash_D,
    thrust::device_vector<bool>& visited_D
);

void makeCHPerHeight_2(int TreeHeight, int cudaThreadNum, int TreeHeightTarget,
    thrust::device_vector<Graph_D_H::TDrank>& nonAdjcentNode_D,
    thrust::device_vector<int>& ID_hash_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& CHTreeHash_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& CHTree_D,
    thrust::device_vector<Graph_D_H::myPair<int>>& NANHash_D,
    thrust::device_vector<bool>& visited_D,
    thrust::device_vector<Graph_D_H::pairs>& CHInfo_D);


void makeH2HLabel(int tempHeight, int cudaThreadNum,
    thrust::device_vector<Graph_D_H::pairs>& H2H_startIndex_D,
    thrust::device_vector<Graph_D_H::H2HLabel>& H2H_label_D,
    thrust::host_vector<int>& Childs,
    thrust::host_vector<Graph_D_H::pairs>& ChildHash,
    thrust::host_vector<int>& frontier
);

void makeH2HLabel_noCommunication(int tempHeight, int cudaThreadNum,
    thrust::device_vector<Graph_D_H::pairs>& H2H_startIndex_D,
    thrust::device_vector<Graph_D_H::H2HLabel>& H2H_label_D,
    thrust::device_vector<int>& TreeBFS_D,
    thrust::host_vector<int>& TreeHash_D
);

double makeH2HLabel_noCommunication_noHub(int tempHeight, int cudaThreadNum,
    thrust::device_vector<int64_t>& H2H_pos_hash_D,
    thrust::device_vector<int64_t>& H2H_dis_hash_D,
    thrust::device_vector<int>& H2H_pos_ID_D,
    thrust::device_vector<int>& H2H_pos_POS_D,
    thrust::device_vector<int>& H2H_dis_D,
    thrust::device_vector<int>& father_D,
    thrust::device_vector<int>& TreeBFS_D,
    thrust::host_vector<int>& TreeHash
);

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
);

double makeH2HLabel_noCommunication_noHub_3(int tempHeight, int cudaThreadNum,
    thrust::device_vector<int64_t>& H2H_dis_hash_D,
    thrust::device_vector<int>& H2H_dis_D,
    thrust::device_vector<int>& TreeBFS_ID_D,
    thrust::device_vector<int>& TreeBFS_adj_D,
    thrust::device_vector<int>& TreeBFS_pos_D,
    //thrust::device_vector<int>& TreeBFS_changeTime_D,
    thrust::device_vector<int>& father_D,

    thrust::host_vector<int64_t>& TreeBFS_Hash
);

long long int H2HQuery_UsingLCA(
    thrust::device_vector<int>& x,
    thrust::device_vector<int>& y,
    thrust::device_vector<int>& result,
    thrust::device_vector<int>& firstAppeare_D,
    thrust::device_vector<Graph_D_H::pairs>& RMQ_OneLine_D,
    thrust::device_vector<Graph_D_H::pairs>& H2H_startIndex_D,
    thrust::device_vector<Graph_D_H::H2HLabel>& H2H_label_D,
    int RMQ_Size, int  RMQ_Line_Size, int cudaThreadNum, int querysize);

long long int H2HQuery_NoLCA(
    thrust::device_vector<int>& x,
    thrust::device_vector<int>& y,
    thrust::device_vector<int>& result,
    thrust::device_vector<Graph_D_H::pairs>& H2H_startIndex_D,
    thrust::device_vector<Graph_D_H::H2HLabel>& H2H_label_D,
     int cudaThreadNum, int querysize);


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
);


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

    int cudaThreadNum, int querysize
);