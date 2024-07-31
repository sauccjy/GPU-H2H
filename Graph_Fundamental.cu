#pragma once
#include"Graph_D.cuh"

Graph_D_H::Graph::Graph(const string graphName , int TreeHeight, const string Nodefile, const string EdgeFile, int threadPoolSize, int changeHeight, string PartitionCutFile, string PartitionIDFile)
{
	this->graphName = graphName;
	this->NodeFile = Nodefile;
	this->EdgeFile = EdgeFile;
	this->changeHeight = changeHeight;
	mainTreeHeight = TreeHeight;
	//this->TreeHeight = TreeHeight;
	this->partitionCut = PartitionCutFile;
	this->partitionID = PartitionIDFile;

	threadNumber = threadPoolSize;
	time_Mine time;

	cout << "graph node start loading" << endl;
	time.updateStart();
	//load Node
	fstream partitionNode(PartitionIDFile, ios::in | ios::out);
	partitionNode >> NodeNumber;
	NE_P.assign(NodeNumber, myPair<double>());
	ID_hash.assign(NodeNumber, -1);
	int Hash = -1, ID = -1;
	double la = 0, gr = 0;
	for (int i = 0; i < NodeNumber; i++) {
		partitionNode >> Hash >> la >> gr >> ID;
		NE_P[Hash].setPairs(ID, la, gr);
		ID_hash[ID] = Hash;
	}
	partitionNode.close();
	time.updateEnd();
	cout << "\t graph node end and using time:" << time.get_microsecond_duration() << endl;

	cout << "graph edge start loading" << endl;
	time.updateStart();
	fstream Edge(EdgeFile, ios::in | ios::out);
	Edge >>NodeNumber>> EdgeNumber;
	int NodeID1 = 0, NodeID2 = 0, EdgeWeight = 0;
	char s;
	//int tempNodeID = 0;
	int actualEdgeNumber = EdgeNumber;
	adjList.assign( NodeNumber , thrust::host_vector<pairs>() );

	for (int i = 0; i < EdgeNumber; i++)
	{
		Edge >>s >> NodeID1 >> NodeID2 >> EdgeWeight;
		if (NodeID1 == NodeID2)
		{
			actualEdgeNumber--;
			continue;
		}
		adjList[NodeID1 - 1].push_back(pairs(NodeID2 - 1, EdgeWeight));
	}
	EdgeNumber = actualEdgeNumber;
	Edge.close();
	time.updateEnd();
	cout << "\t graph edge finished and using time:" << time.get_microsecond_duration() << endl;

	cout << "Partition Cut loading" << endl;
	time.updateStart();

	fstream Partition1(PartitionCutFile, ios::in | ios::out);
	int maxHeight;
	Partition1 >> maxHeight;
	maxHeight = min(maxHeight, mainTreeHeight);
	mainTreeHeight = maxHeight;
	this->TreeHeight = mainTreeHeight;
	long long int fullPar = std::pow(2, maxHeight) - 1;

	partition_Tree.assign(fullPar, pairs());
	int h = 0, left = 0, right = 0;
	for (long long int i = 0; i < fullPar; i++) {
		Partition1 >> h >> left >> right;
		partition_Tree[i].pairsReset(left, right);
	}
	Partition1.close();
	time.updateEnd();
	cout<<"\t Partition Cut loading end and using time: "<< time.get_microsecond_duration() << endl;

	//displayPartition();

	
	CHInfo_H.assign(NodeNumber, pairs(-1,0));
	//Graph_D_H::CHInfo.assign(NodeNumber, pairs(0, 0));
	
	for (int i = 0; i < NodeNumber; i++)
	{
		CHInfo_H[i].first = adjList[i].size();
	}
	//CHInfo_D = CHInfo_H;
	visited.assign(NodeNumber, false);
	CHAdjlist.assign(NodeNumber, unordered_map<int, pairs>());
}


Graph_D_H::Graph::~Graph()
{
}

void Graph_D_H::Graph::checklink()
{
	vector<bool> visit(NodeNumber,false);
	priority_queue<int> que;
	int size = 0;
	que.push(0);
	visit[0] = true;
	while (!que.empty()) {
		int ID = que.top();
		que.pop();
		size++;
		for (int i = 0; i < adjList[ID].size(); i++)
		{
			if (visit[adjList[ID][i].first])
				continue;
			que.push(adjList[ID][i].first);
			visit[adjList[ID][i].first] = true;
		}

	}
	if (size == NodeNumber) {
		cout << "connect!"<<endl;
	}
	else {
		cout << "disconnect!"<<endl;
	}
}



void Graph_D_H::Graph::printConstructInfo()
{
	string constructData = "./ConstructInfo/"+graphName+"/";
	if(ContractM == 1 && ConstructMethod == 1){
		constructData += "CPU-CPU.csv";
	}
	if(ContractM == 2 && ConstructMethod == 1){
		constructData += "GPU-CPU.csv";
	}
	if(ContractM == 1 && ConstructMethod == 2){
		constructData += "CPU-GPU.csv";
	}
	if(ContractM == 2 && ConstructMethod == 2){
		constructData += "GPU-GPU.csv";
	}
	// if (ConstructMethod == 1) {
	// 	constructData += "MultiThread.csv";
	// }
	// else {
	// 	constructData += "GPUHybrid.csv";
	// }

	//string constructData = "constructInfo.csv";
	std::ofstream outputFile(constructData, std::ios::app);
	if (outputFile.is_open()) {

		outputFile << "graphName,PartitionTreeHeight,changeHeight,ParititionTime,makePartitionRankTreeTime,AlllocateLUBTime,LUBSize,translateLUBToGPUTime,TDInGPUTime,translateBackAndMallocToCPUTime,vertexesConstructedInGPU,DecompositionTreeSize,TDInCPUTime,TDTREEHeight,,RMQsize,translateRMQTime,H2HLabelSize,,changeConstructDeviceSize,H2HmallocTime,H2HConstructBFSTreeTime,H2HTranslateBFSTreeTime,H2HUsingTime_CPU,H2HTranslateTime,H2HUsingTime_GPU,PartitionMethod,PartitionMethod\n";

		outputFile << graphName << "," << TreeHeight << "," << changeHeight << "," << partitionTime 
			<< "," << makePartitionRankTreeTime << "," << AllocateLUBTime << ","
			<< LUBSize << "," << translateLUBToGPUTime << "," << TDInGPUTime << "," 
			<< translateBackAndMallocToCPUTime << "," 
			<< vertexesConstructedInGPU << "," << DecompositionTreeSize 
			<< "," << TDInCPUTime << "," << TDTREEHeight << ",," << RMQsize << "," 
			<< translateRMQTime << "," << H2HLabelSize << ",," << changeConstructDeviceSize << "," 
			<< H2HmallocTime << "," << H2HConstructBFSTreeTime << "," << H2HTranslateBFSTreeTime 
			<< "," << H2HUsingTime_CPU << "," << H2HTranslateTime 
			<< "," << H2HUsingTime_GPU <<","<< PartitionMethod << "\n";


		outputFile.close();
		//std::cout << "Data has been written to output.csv\n";
	}
	else {
		std::cerr << "Unable to open file for writing.\n";
	}
}


