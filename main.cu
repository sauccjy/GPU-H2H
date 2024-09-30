
#include"Graph_D.cuh"
using namespace std;



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
	Agraph->beforeH2H_noHub();
	time.updateEnd();
	cout << "\t OULARMQ end, using time: " << time.get_microsecond_duration() << endl;
	//Agraph->displayOULARMQ();
	cout << endl; cout << endl;
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
	cout << "start OULARMQ " << endl;
	time.updateStart();
	Agraph->beforeH2H_noHub();
	time.updateEnd();
	cout << "\t OULARMQ end, using time: " << time.get_microsecond_duration() << endl;
	//Agraph->displayOULARMQ();
	cout << endl; cout << endl;
}

void contraction_noHub(Graph_D_H::Graph* Agraph, const int contractDevice, const int goalHeight, const int changeHeight) {
	if (contractDevice == 0) { //CPU-only (CPU parallel contraction + global contraction)
		Partition_CH_Construct_noHub(Agraph, goalHeight, changeHeight);
		return;
	}
	if (contractDevice == 1) { //GPU
		TDAndH2HConstruct_Tradition_plus_G_2_noHub(Agraph, goalHeight, changeHeight);
		return;
	}
}

void constructLabel_noHub(Graph_D_H::Graph* Agraph, const int parallelFine, const int constructDevice) {

	cout << "_________________start construct H2H______________________" << endl;
	if (parallelFine == 0) { //serial
		Agraph->inConstructH2H_noHub();
		return;
	}
	if (parallelFine == 1) { //BFS fine parallel
		if (constructDevice == 0) { // CPU
			Agraph->inConstructH2H_noHub_multiThread();
			return;
		}
		if (constructDevice == 1) { // GPU
			Agraph->inConstructH2H_noHub_D_2();
			return;
		}
		
	}
	if (parallelFine == 2) { //TD fine parallel
		if (constructDevice == 0) { // CPU
			Agraph->inConstructH2H_noHub_multiThread_2();
			return;
		}
		if (constructDevice == 1) { // GPU
			Agraph->inConstructH2H_noHub_D();
			return;
		}
	}
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
	//if (Agraph->graphName == "NY") {
		//Agraph->answerRealNYQuery_noHub();
	//}
	//else {
		Agraph->answerRandomQuery_noHub();
	//}
	cout << "translate query" << endl;
	time.updateStart();
	Agraph->translateQueryResultBack();
	time.updateEnd();
	std::cout << "using time: " << time.get_microsecond_duration() << "us" << endl;
}



int main(int argc, char** argv)
{
	if (argc != 9)
		return 0;

	//string gname(argv[1]);
	char* endptr;

	string graphName(argv[1]);
	int TreeHeight = (int)strtol(argv[2], &endptr, 10);			//tree height 
	int changeHeight = (int)strtol(argv[3], &endptr, 10);		//change contract order base, from partition order to global order
	int threadPoolSize = (int)strtol(argv[4], &endptr, 10);		//max thread number on CPU
	int partitionMethod = (int)strtol(argv[5], &endptr, 10);	//[1,2,3,4]: 1: XYCoord partition; 2: HC2L-like partition; 3: metis partition; 4: SCOTCH partition

	int contractDevice = (int)strtol(argv[6], &endptr, 10);		//[0,1]: partition order contraction device: 0: CPU; 1: GPU
	int parallelFine = (int)strtol(argv[7], &endptr, 10);		//[0,1,2]: 0: serial; 1: BFS fine; 2: TD fine
	int constructDevice = (int)strtol(argv[8], &endptr, 10);	//[0,1]: construction device: 0: CPU; 1: GPU

	string NodeFile, EdgeFile;
	NodeFile = "./Graph/USA/";
	EdgeFile = "./Graph/USA/";
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
	else if (partitionMethod == 3) {
		PartitionCutFile += "metisCut.txt";
		PartitionIDFile += "metisID.txt";
	}
	else if (partitionMethod == 4) {
		PartitionCutFile += "scotchCut.txt";
		PartitionIDFile += "scotchID.txt";
	}
	cout << PartitionCutFile << endl;
	cout << PartitionIDFile << endl;
	//cout << "myPair size: " << sizeof(Graph_D_H::myPair<int>) << endl;

	Graph_D_H::Graph* Agraph = new Graph_D_H::Graph(graphName, TreeHeight, NodeFile, EdgeFile, threadPoolSize, changeHeight, PartitionCutFile, PartitionIDFile);
	Agraph->PartitionMethod = partitionMethod;
	Agraph->parallelFine_Construct = parallelFine;
	Agraph->contractDevice = contractDevice;
	Agraph->parallelFine = parallelFine;
	Agraph->constructDevice = constructDevice;
	Agraph->checklink();

	contraction_noHub(Agraph, contractDevice, TreeHeight, changeHeight);
	constructLabel_noHub(Agraph, parallelFine, constructDevice);
	printInfo(Agraph);
	Query_noHub(Agraph);
	
}