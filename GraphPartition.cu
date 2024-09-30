#pragma once
#include"Graph_D.cuh"

namespace Graph_D_H
{
	void Graph::makeAdjcentNode(int goalHeight)
	{
		Graph_D_H::time_Mine time;

		cout << "start construct Partition Rank tree " << endl;
		time.updateStart();

		int tempHeight = TreeHeight;
		thrust::host_vector<bool> isAdjcent(NodeNumber, true);
		//thrust::host_vector<pairs> range(NodeNumber, pairs(INT_MAX, INT_MIN));
		thrust::host_vector<int> placeOccupy(NodeNumber, 0);

		for (int i = 0; i < NodeNumber; i++) {
			placeOccupy[i] = adjList[i].size();
		}

		NANHash.assign(partition_Tree.size(), myPair<int>());
		nonAdjcentNode.clear();
		vector<double> afs(TreeHeight + 1, 0);
		while (tempHeight > goalHeight)
		{
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 2;
			thrust::host_vector<bool> isAdjcent_cp = isAdjcent;
			for (int i = lowestIndexEnd; i >= lowestIndexStart; i--)
			{
				NANHash[i].first = nonAdjcentNode.size();//index
				NANHash[i].second = 0; //maxSize
				NANHash[i].NodeID = 0;//nonAdjcent Size
				int left = partition_Tree[i].first;
				int right = partition_Tree[i].second;
				//vector<int> outPartitionSize(right - left, 0);
				map<int, int> subgraph;
				for (int j = left; j < right; j++)//construct sub graph
				{
					if (!isAdjcent[NE_P[j].NodeID])
						continue;
					int nodeID = NE_P[j].NodeID;
					subgraph.emplace(nodeID, 0);
				}
				NANHash[i].second = subgraph.size();//record max size
				for (auto& it : subgraph) {//define outsize
					it.second = subgraph.size();
					for (auto& adj : adjList[it.first]) {
						if (!isAdjcent[adj.first]) {
							continue;
						}
						if (subgraph.find(adj.first) == subgraph.end()) {
							it.second++;
							//isAdjcent_cp[it.first] = false;
						}
					}
				}

				for (auto& it : subgraph) {
					if (it.second == NANHash[i].second) {// mark non-adjNode
						isAdjcent_cp[it.first] = false;
						NANHash[i].NodeID++;
						nonAdjcentNode.push_back(TDrank(adjList[it.first].size(), it.first, 0));
					}
					placeOccupy[it.first] = max(placeOccupy[it.first], it.second);

				}
			}
			isAdjcent = isAdjcent_cp;
			afs[tempHeight] = nonAdjcentNode.size();
			tempHeight--;
		}



		time.updateEnd();
		makePartitionRankTreeTime = time.get_microsecond_duration();
		cout << "\t construct Partition Rank tree end,  using time: " << time.get_microsecond_duration() << endl;
		string latitudePartition = "LatitudePartitionRate.csv";
		if (PartitionMethod == 2) latitudePartition = "MinimumPartitionRate.csv";
		if (PartitionMethod == 3) latitudePartition = "metisPartitionRate.csv";
		if (PartitionMethod == 4) latitudePartition = "scotchPartitionRate.csv";
		std::fstream heightdata(latitudePartition, ios::in | ios::out | ios::app);

		heightdata << graphName << ",";
		////cout << "NodeNumber per layer: \n";
		for (int i = 1; i < afs.size(); i++) {
			heightdata << afs[i] / NodeNumber << ",";
			cout << "\t layer : " << i << " size: " << std::fixed << afs[i] << endl;
		}
		heightdata << "\n";
		heightdata.close();

		//allowcate CHTree
		cout << "start Allocate LUB " << endl;
		time.updateStart();
		CHTree.clear();
		CHTreeHash.assign(NodeNumber, myPair<int>());
		
		for (int i = 0; i < NodeNumber; i++)
		{
			int nodeID = NE_P[i].NodeID;
			CHTreeHash[i].first = 0;//actual size
			CHTreeHash[i].second = placeOccupy[nodeID];//maxSize
			//cout << placeOccupy[i] << endl;
			CHTreeHash[i].NodeID = CHTree.size();//CHTree index
			CHTree.insert(CHTree.end(), CHTreeHash[i].second, myPair<int>(INT_MAX,0,-1)); //insert maxSize label(outpoint,weight,hub)
			
			for (auto& it : adjList[nodeID])
			{
				CHTree[CHTreeHash[i].NodeID + CHTreeHash[i].first].setPairs(-1,it.first,it.second);
				CHTreeHash[i].first++;
				if (CHTreeHash[i].first > CHTreeHash[i].second) {
					cout << "LUB error~! :::" << CHTreeHash[i].first - CHTreeHash[i].second << endl;
					cout << "\t actual adjcent size: " << adjList[nodeID].size() << " lub size: " << CHTreeHash[i].second << endl;
					throw("LUB allocate Error!!");
				}
			}
		}
		
		chSize = CHTree.size();

		time.updateEnd();
		AllocateLUBTime = time.get_microsecond_duration();
		cout << "\t Malloc LUB end,  using time: " << time.get_microsecond_duration() << endl;
		//cout << "LUB Size: " << ((double)(chSize) * sizeof(Graph_D_H::myPair<int>)) / (1024 * 1024) << "MB" << endl;
		LUBSize = ((double)(chSize) * sizeof(Graph_D_H::myPair<int>)) / (1024 * 1024);
		cout << "LUB Size: " << ((double)(chSize) * sizeof(Graph_D_H::myPair<int>)) / (1024 * 1024) << "MB" << endl;
		cout << "LUB length: " << chSize << endl;

		NE_P.clear();
		NE_P.shrink_to_fit();
		partition_Tree.clear();
		partition_Tree.shrink_to_fit();
	}

}
