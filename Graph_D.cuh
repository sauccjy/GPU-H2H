#pragma once
//#pragma auto_inline(off)
#include<iostream>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include<thrust/count.h>
#include<bitset>
#include<fstream>
#include<thrust/sort.h>
#include<thrust/shuffle.h>
#include<thrust/random.h>
#include<thrust/execution_policy.h>
#include<stack>
#include<queue>
#include<algorithm>
#include<chrono>
#include<random>
#include<vector>
#include<map>
#include<set>
#include<unordered_map>
#include <cstdint>
#include <thread>
#include <functional>
#include <future>
#include <mutex>
#include <condition_variable>



using namespace std;
using namespace std::chrono;

namespace Graph_D_H
{
	struct time_Mine {
		std::chrono::high_resolution_clock::time_point start;
		std::chrono::high_resolution_clock::time_point end;
		void updateStart()
		{
			start = std::chrono::high_resolution_clock::now();
		}
		void updateEnd()
		{
			end = std::chrono::high_resolution_clock::now();
		}
		long long int get_microsecond_duration()
		{
			return std::chrono::duration_cast<microseconds>(end - start).count();
		}

	};

	struct H2HLabel
	{
	public:
		int Node = -1;
		bool isTreeNode = false;
		int Hub = -1;
		int dis = INT_MAX;
		__host__ __device__
		H2HLabel() = default;
		__host__ __device__
		H2HLabel(int Node, bool isTreeNode, int Hub, int _dis) : Node(Node), isTreeNode(isTreeNode), Hub(Hub), dis(_dis) { }
		__host__ __device__
		void setH2HLabel(int _Node, bool _isTreeNode, int _Hub, int _dis)
		{
			Node = _Node;
			isTreeNode = _isTreeNode;
			Hub = _Hub;
			dis = _dis;
		}
		__host__ __device__
		void setH2HLabel(H2HLabel& others)
		{
			Node = others.Node;
			isTreeNode = others.isTreeNode;
			Hub = others.Hub;
			dis = others.dis;
		}
	};

	struct pairs
	{
		int first = -1;
		int second = -1;

		__host__ __device__
		pairs() = default;
		__host__ __device__
		pairs(int first, int second)
			: first(first), second(second)
		{
		}
		__host__ __device__
			pairs(const pairs& other)
			: first(other.first), second(other.second)
		{
		}
		__host__ __device__
		void pairsCopy(const pairs& others)
		{
			first = others.first;
			second = others.second;
		}
		__host__ __device__
		void pairsReset(int first_, int second_)
		{
			first = first_;
			second = second_;
		}
		__host__ __device__
		bool operator<(const pairs& other) const {
			if (first < other.first)
				return true;
			else if (first == other.first)
				return second < other.second;
			else
				return false;
		}
		__host__ __device__
			bool operator>(const pairs& other) const {
			if (first > other.first)
				return true;
			else if (first == other.first)
				return second > other.second;
			else
				return false;
		}
	};

	struct TDrank {
		int first = -1;//Node degree
		int second = -1;//Node ID
		int third = 0;  //continus degree
		int childSize = 0;
		__host__ __device__
			TDrank() = default;
		__host__ __device__
			TDrank(int first, int second, int third) : first(first), second(second), third(third) {}
		__host__ __device__
			TDrank(const TDrank& others) : first(others.first), second(others.second), third(others.third) {}
		__host__ __device__
			void setDetal(int f, int s, int t)
		{
			first = f;
			second = s;
			third = t;
		}
		__host__ __device__
			bool operator<(const TDrank& other) const
		{
			//if (third + first + childSize != other.third + other.first + other.childSize) {
				//return third + first + childSize < other.third + other.first + other.childSize;
			//}
			//else {
			//if (third != other.third ) {
				//return third < other.third;
			//}
				if (third + first != other.third + other.first) {
					return third + first < other.third + other.first;
				}
				if (first != other.first) {
					return first < other.first;
				}
				return second < other.second;

				
			//}
		}
		__host__ __device__
			bool operator>(const TDrank& other) const
		{
			//if (third + first + childSize != other.third + other.first + other.childSize) {
				//return third + first + childSize > other.third + other.first + other.childSize;
			//}
			//else {
				if (third + first != other.third + other.first) {
					return third + first > other.third + other.first;
				}
				if (first != other.first) {
					return first > other.first;
				}
				return second > other.second;
					
				
			//}
		}
		__host__ __device__
			void operator=(const TDrank& other)
		{
			this->first = other.first;
			second = other.second;
			third = other.third;

		}
		__host__ __device__
			bool operator==(const TDrank& other)const
		{
			return second == other.second;
		}


	};

	static vector<pairs> CHInfo;
	static vector<int> CHChild;
	static vector<int64_t> DBC;
	struct HostRank {
		int x;
		bool operator<(const HostRank& other) const
		{
			//if ( CHInfo[x].second !=  CHInfo[other.x].second)
				//return CHInfo[x].second < CHInfo[other.x].second;
			if (CHInfo[x].first + CHInfo[x].second != CHInfo[other.x].first + CHInfo[other.x].second)
				return CHInfo[x].first + CHInfo[x].second < CHInfo[other.x].first + CHInfo[other.x].second;
			if (CHInfo[x].first != CHInfo[other.x].first)
				return CHInfo[x].first < CHInfo[other.x].first;
			return x < other.x;
		}

		HostRank(int x) : x(x) {
		}

		HostRank() = default;

	};

	template<typename T>
	struct myPair {
		int NodeID = -1;
		T first = 0;
		T second = 0;
		//__host__ __device__
		myPair() = default;
		__host__ __device__
		myPair(const T first, const T second) : first(first), second(second) {}
		__host__ __device__
		myPair(const T first, const T second, int ID) : first(first), second(second) ,NodeID(ID){}
		__host__ __device__
			myPair(const myPair<T>& other ) : first(other.first), second(other.second), NodeID(other.NodeID) {}
		__host__ __device__
		void setPairs( int ID,const T& _first, const T& _second) {
			NodeID = ID;
			first = _first;
			second = _second;
		}

		__host__ __device__
			bool operator<(const myPair<T>& others)const {
			return first < others.first;
		}
		__host__ __device__
			void operator=(const myPair<T>& other)
		{
			this->first = other.first;
			second = other.second;
			NodeID = other.NodeID;

		}
	};

	struct myPair_latitude_less_than
	{

		__host__ __device__
			bool operator()(const myPair<double> &other, const myPair<double> &other1) const
		{
			return other.second < other1.second;
		}
	};

	struct myPair_longitude_less_than
	{
		__host__ __device__
			bool operator()(const myPair<double> &other, const myPair<double> &other1) const
		{
			return other.first < other1.first;
		}
	};

	class Graph {

	public:

		//main info
		int NodeNumber = 0;
		int EdgeNumber = 0;
		int threadNumber = 16;
		string NodeFile = "";
		string EdgeFile = "";
		string graphName = "";
		std::string partitionCut = "";
		std::string partitionID = "";
		int mainTreeHeight = 17;

		//input info
		int contractDevice = -1;
		int parallelFine = -1;
		int constructDevice = -1;

		//adjList and  NE(longitude and latitude)
		thrust::host_vector<thrust::host_vector<pairs> > adjList = {};
		thrust::host_vector<myPair<double> > NE_P = {};//first : N ; second :E ; nodeID : ID;

		thrust::host_vector<int> ID_hash = {};//hash NE_P
		thrust::host_vector<int> ranks = {};
		int head = -1; //headID,
		//partitionTree
		thrust::host_vector<pairs> partition_Tree = {}; // partition (left,right) index in NE_P , binary tree structure
		int TreeHeight = 1;
		int lowestMaxSize = 0;
		int CHTreeHeight = 0;
		//order and partition
		//Partition Rank Tree
		int PartitionMethod = 1;
		thrust::host_vector<TDrank> nonAdjcentNode = {};// ID = nonAdjcentNode[rank].second, from rank to ID
		thrust::host_vector<myPair<int>> NANHash = {}; //first : nonAdjcentNode index; second maxSize; nodeID : nonadjcent Size
		thrust::host_vector<int> father = {};
		// first: childNumber, second: childsIndex; 
		thrust::host_vector<pairs> ChildHash = {};
		//childs ID; childs[ChildHash[ID].second ] to childs[ChildHash[ID].second + ChildHash[ID].first - 1] 
		thrust::host_vector<int> Childs = {};

		//first:actual size; second:maxSize; nodeID:CHTree index; Following ID_hash;
		thrust::host_vector<myPair<int>> CHTreeHash = {};
		thrust::host_vector<bool> visited = {};
		thrust::host_vector<myPair<int>> CHTree = {}; //first: outPoint; second: weight; NodeID: hub;
		long long int chSize = 0;//in byte
		//oula sequence and RMQ
		thrust::host_vector<pairs > OULA_DFS_ = {};//first: nodeID ; second : depth
		thrust::host_vector<int> firstAppeare = {};
		thrust::host_vector<pairs> RMQ_OneLine = {};

		long long int RMQ_Size = -1;
		int RMQ_Line_Size = -1;
		thrust::host_vector<int> OULA_Only = {};
		bool oulaOnly = false;
		//H2H
		long long int H2H_Size = 0;
		thrust::host_vector<pairs> CHInfo_H = {};
		//H2H with hub for Path Backtracking
		thrust::host_vector<pairs> H2H_startIndex = {}; //label start in H2H_label, size
		thrust::host_vector<H2HLabel> H2H_label = {}; //
		thrust::host_vector<int> TreeBFS = {};
		thrust::host_vector<int> TreeHash = {};

		//Device CSR
		thrust::device_vector<pairs> CHInfo_D = {};
		thrust::device_vector<int> ID_hash_D;//
		thrust::device_vector<int> ranks_D ;//
		thrust::device_vector<TDrank> nonAdjcentNode_D;//
		thrust::device_vector<myPair<int>> NANHash_D;//
		thrust::device_vector<myPair<int>> CHTreeHash_D;//
		thrust::device_vector<bool> visited_D;
		thrust::device_vector<myPair<int>> CHTree_D;//
		thrust::device_vector<int> firstAppeare_D ;
		thrust::device_vector<pairs > RMQ_OneLine_D ;
		thrust::device_vector<int> OULA_Only_D;
		thrust::device_vector<pairs> H2H_startIndex_D;
		thrust::device_vector<H2HLabel> H2H_label_D;
		thrust::device_vector<int> TreeBFS_D;
		thrust::device_vector<int> father_D;

		//no hub label
		thrust::host_vector<int64_t> H2H_pos_hash = {};
		thrust::host_vector<int64_t> H2H_dis_hash = {};
		thrust::host_vector<int> H2H_pos_ID = {};
		thrust::host_vector<int> H2H_pos_POS = {};
		thrust::host_vector<int> H2H_dis = {};

		thrust::host_vector<int64_t> TreeBFS_Hash = {};
		thrust::host_vector<int> TreeBFS_ID = {};
		thrust::host_vector<int> TreeBFS_adj = {};
		thrust::host_vector<int> TreeBFS_pos = {};
		thrust::host_vector<int> TreeBFS_changeTime = {};

		thrust::host_vector<int64_t> RMQHash = {};
		//thrust::host_vector<int> RMQ_Height = {};
		thrust::host_vector<int> RMQ_ID = {};
		thrust::host_vector<int> H2H_TreeHeight_Hash = {};

		thrust::device_vector<int64_t> H2H_pos_hash_D = {};
		thrust::device_vector<int64_t> H2H_dis_hash_D = {};
		thrust::device_vector<int> H2H_pos_ID_D = {};
		thrust::device_vector<int> H2H_pos_POS_D = {};
		thrust::device_vector<int> H2H_dis_D = {};

		//thrust::device_vector<int64_t> TreeBFS_Hash_D = {};
		thrust::device_vector<int> TreeBFS_ID_D = {};
		thrust::device_vector<int> TreeBFS_adj_D = {};
		thrust::device_vector<int> TreeBFS_pos_D = {};
		thrust::device_vector<int> TreeBFS_changeTime_D = {};

		thrust::device_vector<int64_t> RMQHash_D;
		//thrust::device_vector<int> RMQ_Height_D;
		thrust::device_vector<int> RMQ_ID_D;
		thrust::device_vector<int> H2H_TreeHeight_Hash_D = {};

		int changeToGPUHeight = -1;
		thrust::host_vector<int> frontier_Hash = {};
		thrust::host_vector<int> frontier_ID = {};

		thrust::device_vector<int> frontier_Hash_D;
		thrust::device_vector<int> frontier_ID_D;

		void beforeH2H_noHub();
		void inConstructH2H_noHub();
		void inConstructH2H_noHub_multiThread();
		void inConstructH2H_noHub_multiThread_2();
		void inConstructH2H_noHub_D();//GPU + BFS
		void inConstructH2H_noHub_D_2(); //GPU + TD
		void mallocH2HLabel_noHub();
		void makeH2HLabel_noHub_serial();
		void makeH2HLabel_noHub_multiThred();
		void makeH2HLabel_noHub_multiThred_2();
		void makeH2HLabel_noHub_noComm();  //GPU + BFS
		void makeH2HLabel_noHub_noComm_2();
		void makeH2HLabel_noHub_noComm_3();//GPU + TD
		void translateH2H_noHub();
		void translateH2HBack_noHub();
		void displayH2H_noHub();
		void translateRMQ_noHub();
		//no Hub query
		int findLCA_noHub(int x, int y);
		int H2HDistancQuery_UsingLCA_noHub(int x, int y);
		int H2HDistancQuery_noLCA_noHub(int x, int y);
		long long int H2Hquery_bunch_noHub(int size);
		long long int H2Hquery_bunch_MultiThread_noHub(int size);
		long long int H2Hquery_bunch_UsingLCA_D_noHub(int size);
		long long int H2Hquery_bunch_noLCA_D_noHub(int size);
		void answerRandomQuery_noHub();
		void answerRealNYQuery_noHub();


		int cudaThreadNum = 1024;


		//construct data and time info
		long long int size1 = 0;
		int parallelFine_Construct = -1;
		int ConstructMethod = 1;
		int changeHeight = 0;
		long long int partitionTime = 0;
		long long int makePartitionRankTreeTime = 0;
		long long int AllocateLUBTime = 0;
		double LUBSize = 0;
		long long int translateLUBToGPUTime = 0;
		long long int TDInGPUTime = 0;
		long long int translateBackAndMallocToCPUTime = 0;
		int vertexesConstructedInGPU = 0;
		double DecompositionTreeSize = 0;
		long long int TDInCPUTime = 0;
		int TDTREEHeight = 0;
		double RMQsize = 0;
		double H2HLabelSize = 0;
		int changeConstructDeviceSize = 64;

		long long int H2HmallocTime = 0;
		long long int H2HConstructBFSTreeTime = 0;
		long long int H2HTranslateBFSTreeTime = 0;
		long long int H2HUsingTime_CPU = 0;
		long long int H2HTranslateTime = 0;
		long long int H2HUsingTime_GPU = 0;

		long long int translateRMQTime = 0;

		//query info
		long long int queryTranslateTime = 0;
		int queryStep = 1000;
		vector<long long int> QueryTime_seq = {};
		vector<long long int> QueryTime_MT_4 = {};
		vector<long long int> QueryTime_MT_6 = {};
		vector<long long int> QueryTime_MT_8 = {};
		vector<long long int> QueryTime_MT_10 = {};
		vector<long long int> QueryTime_MT_12 = {};
		vector<long long int> QueryTime_MT_14 = {};
		vector<long long int> QueryTime_MT_16 = {};
		vector<long long int> QueryTime_noLCA_64 = {};
		vector<long long int> QueryTime_LCA_64 = {};
		vector<long long int> QueryTime_noLCA_128 = {};
		vector<long long int> QueryTime_LCA_128 = {};
		vector<long long int> QueryTime_noLCA_256 = {};
		vector<long long int> QueryTime_LCA_256 = {};
		vector<long long int> QueryTime_noLCA_512 = {};
		vector<long long int> QueryTime_LCA_512 = {};
		vector<long long int> QueryTime_noLCA_1024 = {};
		vector<long long int> QueryTime_LCA_1024 = {};

		//CH using map
		set<HostRank> heap;
		//set<TDrank> Heap = {};
		vector<unordered_map<int, pairs>> CHAdjlist = {};
		int tempnonAdjcentNodeIndex = INT_MAX;

		//query
		int queryMaxSize = 1000000;
		thrust::host_vector<int> x = {};
		thrust::host_vector<int> y = { };
		thrust::host_vector<int> result = {};

		thrust::device_vector<int> x_D = {};
		thrust::device_vector<int> y_D = { };
		thrust::device_vector<int> result_D = {};

		//fundamental
		Graph(const string graphName, const int TreeHeight, const string Nodefile, const string EdgeFile, int threadPoolSize, int changeHeight, string PartitionCutFile, string PartitionIDFile);
		~Graph();
		void checklink();
		
		//LUB allocate and Tree-Decomposition
		void makeAdjcentNode(int goalHeight);
		//GPU-CPU hybrid contraction
		void partition_TD_Parallel_D(int height);
		void partition_TD_Parallel_D_2(int height);
		//traditional CH after GPU 
		void translateInTD();
		void TDInCPU(int tempheight); //if there is no vertices contracted on GPU, such function becomes only -DBC optimized H2H
		void uploadCHset(int tempheight);
		void CHUsingHeap();
		void CHUsingHeap_2();
		//CPU-multi-thread contraction
		void partition_TD_T(int height);
		void TDInCPUMultiThread(int tempheight);


		//order and tree structure
		bool useRMQ = true;
		void beforeH2H();
		void getrank();
		void getTreeStructure();
		void onlyOULA();
		void onlyOULADFS(const int headNow, const int depth, thrust::host_vector<bool>& visited);
		void makeOULAAndRMQ();
		void OULADFS(const int headNow, const int depth, thrust::host_vector<bool>& visited);



		//construct H2H
		void inConstructH2H_MultiThread();
		void inConstructH2H_D();
		void inConstructH2H_mix();
		void inConstructH2H_noComm();
		void mallocH2HLabel();
		void makeH2Hlabel();
		void makeH2Hlabel_MultiThread();
		void makeH2Hlabel_D();
		void makeH2Hlabel_noComm();
		void makeH2Hlabel_mix();
		int findLCA(int x,int y);

		void constructBFSTree();


		//data translate
		void translateBeforeCH();
		void translateCHTree();
		void translateOULARMQ();
		void translateH2HLabel();
		void translateH2HLabelBack();
		void translateBeforeQuery();
		void translateQuery();
		void translateQueryResultBack();

		//data maintainance
		void cleanAfterTreeStruct();
		void cleanBeforeQuery();
		void cleanResult();
		void cleanBeforeConstruct();
		//display
		void displayPartition();//NE_P partition
		void displayOrderAndCHTree();
		void displayCHTree();
		void displayOULARMQ();
		void displayH2HLabel();
		void displaySingalLabel(int id);
		void displayAdjList();
		void displayCHT();
		void displayCHTree_mapVersion();
		void displayGraph(std::map<int, std::map<int, int>>& standard_graph_in_partition);
		void displayCHAdjList();

		//query
		void generateQuery();
		void uploadRealQuery();
		int H2HDistancQuery_UsingLCA(int x, int y);
		int H2HDistancQuery_UsingOULA(int x, int y);
		int H2HDistancQuery_AncestorToPost(int x, int y);
		long long int H2Hquery_bunch(int size);
		long long int H2Hquery_bunch_MultiThread(int size);
		long long int H2Hquery_bunch_UsingLCA_D(int size);
		long long int H2Hquery_bunch_noLCA_D(int size);

		void answerRandomQuery();
		void answerRealNYQuery();

		//others
		void startGPU();
		void printConstructInfo();


		//void cleanBeforeQuery();
	};

}