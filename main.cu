
#include"Graph_D.cuh"
using namespace std;

bool usingGPUConstruct = true;

void TDAndH2HConstruct_Tradition_plus_G_2(Graph_D_H::Graph* Agraph, int goalHeight, int changeHeight)
{
	Agraph->ConstructMethod = 2;
	Agraph->startGPU();
	//Agraph->displayAdjList();
	Graph_D_H::time_Mine time;
	//Agraph->displayPartition();
	Agraph->makeAdjcentNode(0);
	//Agraph->displayOrderAndCHTree();

	cout << endl; cout << endl;
	cout << "________________TD start!_____________________" << endl;
	if (goalHeight != changeHeight)
	{
		cout << "translate before CH start " << endl;
		time.updateStart();
		Agraph->translateBeforeCH();
		time.updateEnd();
		Agraph->translateLUBToGPUTime = time.get_microsecond_duration();
		cout << "\t translate before CH end, using time: " << time.get_microsecond_duration() << endl;

		cout << "TD in GPU start " << endl;
		time.updateStart();
		Agraph->partition_TD_Parallel_D_2(changeHeight);
		time.updateEnd();
		Agraph->TDInGPUTime = time.get_microsecond_duration();
		cout << "\t TD in GPU using time: " << time.get_microsecond_duration() << endl;
		Agraph->translateInTD();
	}
	//Agraph->displayCHTree();
	Agraph->TDInCPU(changeHeight);
	cout << endl;
	//Agraph->displayCHAdjList();
	cout << "start OULARMQ " << endl;
	time.updateStart();
	Agraph->beforeH2H();
	time.updateEnd();
	cout << "\t OULARMQ end, using time: " << time.get_microsecond_duration() << endl;
	//Agraph->displayOULARMQ();
	cout << endl; cout << endl;
	cout << "_________________start construct H2H______________________" << endl;

	//Agraph->inConstructH2H_mix();
	Agraph->inConstructH2H_noComm();
}

void Partition_CH_Construct(Graph_D_H::Graph* Agraph, int goalHeight, int changeHeight) {
	Agraph->startGPU();
	//Agraph->displayAdjList();
	Graph_D_H::time_Mine time;

	Agraph->makeAdjcentNode(0);
	//Agraph->displayOrderAndCHTree();
	cout << endl; cout << endl;
	cout << "________________TD start!_____________________" << endl;
	Agraph->TDInCPUMultiThread(changeHeight);
	//Agraph->displayCHAdjList();

	cout << endl;
	//Agraph->displayCHAdjList();
	//time.updateEnd();
	//cout << " using time: " << time.get_microsecond_duration() << endl;
	cout << "start OULARMQ " << endl;
	time.updateStart();
	Agraph->beforeH2H();
	time.updateEnd();
	cout << "\t OULARMQ end, using time: " << time.get_microsecond_duration() << endl;
	//Agraph->displayOULARMQ();
	cout << endl; cout << endl;
	cout << "_________________start construct H2H______________________" << endl;

	Agraph->inConstructH2H_MultiThread();

	//Agraph->displayH2HLabel();
}


void TDAndH2HConstruct_Tradition_plus_G_2_noHub(Graph_D_H::Graph* Agraph, int goalHeight, int changeHeight)
{

	Agraph->startGPU();
	//Agraph->displayAdjList();
	Graph_D_H::time_Mine time;
	//Agraph->displayPartition();
	Agraph->makeAdjcentNode(0);
	//Agraph->displayOrderAndCHTree();

	cout << endl; cout << endl;
	cout << "________________TD start!_____________________" << endl;
	if (goalHeight != changeHeight)
	{
		cout << "translate before CH start " << endl;
		time.updateStart();
		Agraph->translateBeforeCH();
		time.updateEnd();
		Agraph->translateLUBToGPUTime = time.get_microsecond_duration();
		cout << "\t translate before CH end, using time: " << time.get_microsecond_duration() << endl;

		cout << "TD in GPU start " << endl;
		time.updateStart();
		Agraph->partition_TD_Parallel_D_2(changeHeight);
		time.updateEnd();
		Agraph->TDInGPUTime = time.get_microsecond_duration();
		cout << "\t TD in GPU using time: " << time.get_microsecond_duration() << endl;
		Agraph->translateInTD();
	}
	//Agraph->displayCHTree();
	Agraph->TDInCPU(changeHeight);
	cout << endl;
	//Agraph->displayCHAdjList();
	cout << "start OULARMQ " << endl;
	time.updateStart();
	Agraph->beforeH2H_noHub();
	time.updateEnd();
	cout << "\t OULARMQ end, using time: " << time.get_microsecond_duration() << endl;
	//Agraph->displayOULARMQ();
	cout << endl; cout << endl;
	cout << "_________________start construct H2H______________________" << endl;

	//Agraph->inConstructH2H_noHub();
	//Agraph->inConstructH2H_noHub_multiThread();
	if(	usingGPUConstruct ){
		Agraph->inConstructH2H_noHub_D();
			Agraph->ConstructMethod = 2;
	} else {
		Agraph->inConstructH2H_noHub_multiThread();
			Agraph->ConstructMethod = 1;
	}
	//Agraph->inConstructH2H_noHub_D();

	//Agraph->displayH2H_noHub();
}

void Partition_CH_Construct_noHub(Graph_D_H::Graph* Agraph, int goalHeight, int changeHeight) {
	Agraph->startGPU();
	//Agraph->displayAdjList();
	Graph_D_H::time_Mine time;

	Agraph->makeAdjcentNode(0);
	//Agraph->displayOrderAndCHTree();
	cout << endl; cout << endl;
	cout << "________________TD start!_____________________" << endl;
	Agraph->TDInCPUMultiThread(changeHeight);
	//Agraph->displayCHAdjList();

	cout << endl;
	//Agraph->displayCHAdjList();
	//time.updateEnd();
	//cout << " using time: " << time.get_microsecond_duration() << endl;
	cout << "start OULARMQ " << endl;
	time.updateStart();
	Agraph->beforeH2H_noHub();
	time.updateEnd();
	cout << "\t OULARMQ end, using time: " << time.get_microsecond_duration() << endl;
	//Agraph->displayOULARMQ();
	cout << endl; cout << endl;
	cout << "_________________start construct H2H______________________" << endl;

	//Agraph->inConstructH2H_noHub();
	//Agraph->inConstructH2H_noHub_multiThread();
	//Agraph->inConstructH2H_noHub_D();
	if(	usingGPUConstruct ){
		Agraph->inConstructH2H_noHub_D();
			Agraph->ConstructMethod = 2;
	} else {
		Agraph->inConstructH2H_noHub_multiThread();
			Agraph->ConstructMethod = 1;
	}


	//Agraph->displayH2H_noHub();
}



void printInfo(Graph_D_H::Graph* Agraph) {
	Agraph->printConstructInfo();
}

void Query(Graph_D_H::Graph* Agraph) {
	Agraph->cleanBeforeQuery();
	Graph_D_H::time_Mine time;
	cout << "translate RMQ" << endl;
	time.updateStart();
	Agraph->translateRMQ_noHub();
	time.updateEnd();
	Agraph->translateRMQTime = time.get_microsecond_duration();
	cout << "\t translate RMQ using time: " << time.get_microsecond_duration() << endl;
	if (Agraph->graphName == "NY") {
		Agraph->answerRealNYQuery();
	}
	else {
		Agraph->answerRandomQuery();
	}
	cout << "translate query" << endl;
	time.updateStart();
	Agraph->translateQueryResultBack();
	time.updateEnd();
	cout << "using time: " << time.get_microsecond_duration() << "us" << endl;

}

void Query_noHub(Graph_D_H::Graph* Agraph) {
	Agraph->cleanBeforeQuery();
	Graph_D_H::time_Mine time;

	cout << "translate RMQ" << endl;
	time.updateStart();
	Agraph->translateRMQ_noHub();
	time.updateEnd();
	Agraph->translateRMQTime = time.get_microsecond_duration();
	cout << "\t translate RMQ using time: " << time.get_microsecond_duration() << endl;
	//Agraph->translateQuery();
	if (Agraph->graphName == "NY") {
		Agraph->answerRealNYQuery_noHub();
	}
	else {
		Agraph->answerRandomQuery_noHub();
	}
	cout << "translate query" << endl;
	time.updateStart();
	Agraph->translateQueryResultBack();
	time.updateEnd();
	cout << "using time: " << time.get_microsecond_duration() << "us" << endl;
}



int main(int argc, char** argv)
{
	if(argc != 7)
		return;
	string gname(argv[1]);
	char *endptr;
	
    int secondArgLong = (int)strtol(argv[2], &endptr, 10);//tree height
	int thirdArgLong = (int)strtol(argv[3], &endptr, 10); //change device height(in contract, from GPU to CPU)
	int forthArgLong = (int)strtol(argv[4], &endptr, 10); //device to contract 0:CPU, 1:GPU;
	int fifthArgLong = (int)strtol(argv[5], &endptr, 10); //device to construct 0:CPU, 1:GPU;
	int sixthArgLong = (int)strtol(argv[6], &endptr, 10); //is Query?  0:no query, 1: query;

	string NodeFile, EdgeFile;
	NodeFile = "./Graph/USA/";
	EdgeFile = "./Graph/USA/";

	int threadPoolSize = 16;
	int TreeHeight = secondArgLong;                                                               
	int changeHeight =thirdArgLong;
	string graphName = gname;

	int partitionMethod = 1;

	bool printConInfo = 1;

	usingGPUConstruct = fifthArgLong;
	int ConstructM = forthArgLong + 1;
	bool isQuery = sixthArgLong; 


	bool noHub = true;

	NodeFile += graphName;
	EdgeFile += graphName;
	NodeFile += "Node.txt";
	EdgeFile += "Edge.txt";
	cout << NodeFile << endl;
	cout << EdgeFile << endl;

	string PartitionCutFile, PartitionIDFile;
	PartitionCutFile = "./Partition/" + graphName + "/";
	PartitionIDFile = PartitionCutFile;
	if (partitionMethod == 2) {
		PartitionCutFile += "MinimumCut.txt";
		PartitionIDFile += "MinimumID.txt";
	}
	else if (partitionMethod == 1) {
		PartitionCutFile += "LatitudeCut.txt";
		PartitionIDFile += "LatitudeID.txt";
	}

	cout << PartitionCutFile << endl;
	cout << PartitionIDFile << endl;
	//cout << "myPair size: " << sizeof(Graph_D_H::myPair<int>) << endl;

	Graph_D_H::Graph* Agraph = new Graph_D_H::Graph(graphName, TreeHeight, NodeFile, EdgeFile, threadPoolSize, changeHeight, PartitionCutFile, PartitionIDFile);
	Agraph->PartitionMethod = partitionMethod;
	Agraph->ContractM = ConstructM;
	Agraph->checklink();
	if(!noHub)
	{
		if (ConstructM == 1) {
			Partition_CH_Construct(Agraph, TreeHeight, changeHeight);
			//Partition_CH_Construct_noHub(Agraph, TreeHeight, changeHeight);
		}
		else if (ConstructM == 2) {
			TDAndH2HConstruct_Tradition_plus_G_2(Agraph, TreeHeight, changeHeight);
		}
		if (printConInfo) {
			printInfo(Agraph);
		}
		if (isQuery) {
			Query(Agraph);
		}
	}
	else {
		if (ConstructM == 1) {
			//Partition_CH_Construct(Agraph, TreeHeight, changeHeight);
			Partition_CH_Construct_noHub(Agraph, TreeHeight, changeHeight);
		}
		else if (ConstructM == 2) {
			//TDAndH2HConstruct_Tradition_plus_G_2(Agraph, TreeHeight, changeHeight);
			TDAndH2HConstruct_Tradition_plus_G_2_noHub(Agraph, TreeHeight, changeHeight);
		}
		if (printConInfo) {
			printInfo(Agraph);
		}
		if (isQuery) {
			Query_noHub(Agraph);
		}
	}

}