#pragma once
#include"Graph_D.cuh"

namespace Graph_D_H
{
	void Graph_D_H::Graph::TDInCPU(int tempheight)
	{
		Graph_D_H::time_Mine time;
		std::cout << "translate and malloc back to CPU start " << endl;
		time.updateStart();
		uploadCHset(tempheight);

		time.updateEnd();
		translateBackAndMallocToCPUTime = time.get_microsecond_duration();
		std::cout << "\t translate and malloc back to CPU end, using time: " << time.get_microsecond_duration() << endl;

		std::cout << "\n" << "TD in CPU start " << endl;
		time.updateStart();
		CHUsingHeap_2();
		time.updateEnd();
		TDInCPUTime = time.get_microsecond_duration();
		std::cout << "\t TD in CPU end, using time: " << time.get_microsecond_duration() << endl;


	}


	void Graph::uploadCHset(int tempheight)
	{
		for (auto& it : adjList) {
			it.clear();
			it.shrink_to_fit();
		}
		adjList.clear();
		adjList.shrink_to_fit();
		int tempHeight = tempheight;
		CHInfo.assign(NodeNumber, pairs());
		CHChild.assign(NodeNumber, 0);
		DBC.assign(NodeNumber, 1);
		//bool haveless = false;
		//int wrongSize = 0;
		cout << "clear finished!" << endl;
		while (tempHeight > 0)
		{
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 2;
			tempHeight--;

			for (int i = lowestIndexEnd; i >= lowestIndexStart; i--)
			{
				int startIndex = NANHash[i].first;
				int TDNodeSize = NANHash[i].NodeID;
				tempnonAdjcentNodeIndex = min(startIndex, tempnonAdjcentNodeIndex);
				for (int j = startIndex; j < startIndex + TDNodeSize; j++)
				{
					int TDNodeNow = nonAdjcentNode[j].second;
					int CHTreeHash_index = ID_hash[TDNodeNow];
					int CHTree_Index = CHTreeHash[CHTreeHash_index].NodeID;
					int CHTree_ActualSize = CHTreeHash[CHTreeHash_index].first;
					int size = 0;
					for (int k = CHTree_Index; k < CHTree_Index + CHTree_ActualSize; k++)
					{
						if (CHTree[k].first == INT_MAX)
							continue;
						size++;
						CHAdjlist[TDNodeNow].emplace(CHTree[k].first, pairs(CHTree[k].second, CHTree[k].NodeID));
						//CHTree[k].first = INT_MAX;
					}

					Graph_D_H::CHInfo[TDNodeNow].pairsReset(size, CHInfo_H[TDNodeNow].second);
					heap.insert(HostRank(TDNodeNow));

				}
			}
		}
		tempHeight = TreeHeight;
		while (tempHeight > tempheight)
		{
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 2;
			tempHeight--;
			for (int i = lowestIndexEnd; i >= lowestIndexStart; i--)
			{
				int startIndex = NANHash[i].first;
				int TDNodeSize = NANHash[i].NodeID;
				for (int j = startIndex; j < startIndex + TDNodeSize; j++)
				{
					int TDNodeNow = nonAdjcentNode[j].second;
					int CHTreeHash_index = ID_hash[TDNodeNow];
					int CHTree_Index = CHTreeHash[CHTreeHash_index].NodeID;
					int CHTree_ActualSize = CHTreeHash[CHTreeHash_index].first;
					//int size = 0;
					for (int k = CHTree_Index; k < CHTree_Index + CHTree_ActualSize; k++)
					{
						if (CHTree[k].first == INT_MAX)
							continue;
						//size++;
						CHAdjlist[TDNodeNow].emplace(CHTree[k].first, pairs(CHTree[k].second, CHTree[k].NodeID));
						//CHTree[k].first = INT_MAX;
					}
					//Graph_D_H::CHInfo[TDNodeNow].pairsReset(size, 0);
				}
			}
		}
		CHTree.clear();
		CHTree.shrink_to_fit();
		CHTreeHash.clear();
		CHTreeHash.shrink_to_fit();
		CHTree_D.clear();
		CHTree_D.shrink_to_fit();
		CHTreeHash_D.clear();
		CHTreeHash_D.shrink_to_fit();
		CHInfo_H.clear();
		CHInfo_H.shrink_to_fit();
		CHInfo_D.clear();
		CHInfo_D.shrink_to_fit();
		std::cout << " \t node in heap now: " << heap.size() << endl;
	}


	void Graph::CHUsingHeap()
	{
		long long int tempNANIndex = tempnonAdjcentNodeIndex;
		vertexesConstructedInGPU = tempnonAdjcentNodeIndex;
		std::cout<<"\t vertexes Constructed in lower layer: " << tempNANIndex << endl;
		int total = 0;
		while (!heap.empty())
		{
			int TDNodeNow = heap.begin()->x;
			heap.erase(heap.begin());

			if (heap.size() < 10000) {
				std::cout << "file in " << endl;
			}

			total++;
			vector<pairs> info;
			vector<int> nei;
			nonAdjcentNode[tempNANIndex++].setDetal(CHAdjlist[TDNodeNow].size(), TDNodeNow, Graph_D_H::CHInfo[TDNodeNow].second);
			for (auto i = CHAdjlist[TDNodeNow].begin(); i != CHAdjlist[TDNodeNow].end(); i++)
			{
				nei.push_back(i->first);
				info.push_back(i->second);
			}
			for (int i = 0; i < nei.size(); i++)
			{
				int neiID = nei[i];
				auto it = CHAdjlist[neiID].find(TDNodeNow);
				if (it != CHAdjlist[neiID].end())
				{
					CHAdjlist[neiID].erase(it);
				}
			}
			for (int i = 0; i < nei.size(); i++)
			{
				int leftNode = nei[i];
				for (int j = i + 1; j < nei.size(); j++)
				{
					int rightNode = nei[j];
					int weight = info[i].first + info[j].first;
					if (CHAdjlist[leftNode].find(rightNode) == CHAdjlist[leftNode].end())
					{
						CHAdjlist[leftNode].insert(make_pair(rightNode, pairs(weight, TDNodeNow)));
					}
					else
					{
						if (CHAdjlist[leftNode][rightNode].first > weight)
						{
							CHAdjlist[leftNode][rightNode].first = weight;
							CHAdjlist[leftNode][rightNode].second = TDNodeNow;
						}
					}
					if (CHAdjlist[rightNode].find(leftNode) == CHAdjlist[rightNode].end())
					{
						CHAdjlist[rightNode].insert(make_pair(leftNode, pairs(weight, TDNodeNow)));
					}
					else
					{
						if (CHAdjlist[rightNode][leftNode].first > weight)
						{
							CHAdjlist[rightNode][leftNode].first = weight;
							CHAdjlist[rightNode][leftNode].second = TDNodeNow;
						}
					}
				}
			}
			
			for (int i = 0; i < nei.size(); i++)
			{
				int neiID = nei[i];
				int height = max(Graph_D_H::CHInfo[TDNodeNow].second + 1, Graph_D_H::CHInfo[neiID].second);
				if (height == Graph_D_H::CHInfo[neiID].second) {
					CHChild[neiID]++;
				}
				else {
					CHChild[neiID] = 1;
				}
				Graph_D_H::CHInfo[neiID].pairsReset(CHAdjlist[neiID].size(), height);
			}
			
		}

		int head = nonAdjcentNode[nonAdjcentNode.size() - 1].second;
		std::cout << "head: "<< head <<"and  CH adjlist size: " << CHAdjlist[head].size() << endl;
		long long int size1 = 0;
		for (auto& it : CHAdjlist) {
			size1 += it.size() + 1;
		}
		std::cout << "vertexes constructed in CPU: " << total << endl;
		DecompositionTreeSize = ((double)(size1 * sizeof(Graph_D_H::myPair<int>))) / (1024 * 1024);
		std::cout << "\t actual Decomposition Tree size = " << ((double)(size1* sizeof(Graph_D_H::myPair<int>)))/(1024*1024)<<"MB" << endl;
		std::cout << endl;
	}


	void Graph::CHUsingHeap_2()
	{
		long long int tempNANIndex = tempnonAdjcentNodeIndex;
		vertexesConstructedInGPU = tempnonAdjcentNodeIndex;
		std::cout << "\t vertexes Constructed in GPU: " << tempNANIndex << endl;
		int total = 0;
		while (!heap.empty())
		{

			int TDNodeNow = heap.begin()->x;
			heap.erase(heap.begin());
			if (visited[TDNodeNow])
			{
				std::cout << "constructNode: " << TDNodeNow << " ~" << endl;
			}
			total++;
			visited[TDNodeNow] = true;

			vector<pairs> info;
			vector<int> nei;
			nonAdjcentNode[tempNANIndex++].setDetal(CHAdjlist[TDNodeNow].size(), TDNodeNow, Graph_D_H::CHInfo[TDNodeNow].second);
			//std::cout <<"nodeID: "<<TDNodeNow << " adjlist size: " << CHAdjlist[TDNodeNow].size() ;
			//std::cout << "heap size: " << heap.size();
				//int CHTreeHash_index = ID_hash[TDNodeNow];
				//int CHTree_Index = CHTreeHash[CHTreeHash_index].NodeID;

			for (auto i = CHAdjlist[TDNodeNow].begin(); i != CHAdjlist[TDNodeNow].end(); i++)
			{
				//if (visited[i->first])
				//{
				//	cout << "stopedID: " << TDNodeNow<<"adjNode: "<<i->first << endl;
				//	cin.get();
				//	continue;
				//}
				nei.push_back(i->first);
				info.push_back(i->second);
				//CHTree[CHTree_Index++].setPairs((i->second).second, i->first, (i->second).first);
			}

			//firstempty[CHTreeHash_index] = CHAdjlist[TDNodeNow].size();
			//CHTreeHash[CHTreeHash_index].first = CHAdjlist[TDNodeNow].size();
			//std::cout << "erase id: ";
			for (int i = 0; i < nei.size(); i++)
			{
				int neiID = nei[i];
				auto it = CHAdjlist[neiID].find(TDNodeNow);
				if (it != CHAdjlist[neiID].end())
				{
					CHAdjlist[neiID].erase(it);
				}
				//CHAdjlist[neiID].erase(TDNodeNow);
				//std::cout << " " << neiID;
				heap.erase(HostRank(neiID));
			}
			//std::cout << "; erase adj , heapSize: " << heap.size() << "";
			for (int i = 0; i < nei.size(); i++)
			{
				int leftNode = nei[i];
				for (int j = i + 1; j < nei.size(); j++)
				{
					int rightNode = nei[j];
					int weight = info[i].first + info[j].first;
					if (CHAdjlist[leftNode].find(rightNode) == CHAdjlist[leftNode].end())
					{
						CHAdjlist[leftNode].insert(make_pair(rightNode, pairs(weight, TDNodeNow)));
					}
					else
					{
						if (CHAdjlist[leftNode][rightNode].first > weight)
						{
							CHAdjlist[leftNode][rightNode].first = weight;
							CHAdjlist[leftNode][rightNode].second = TDNodeNow;
						}
					}
					if (CHAdjlist[rightNode].find(leftNode) == CHAdjlist[rightNode].end())
					{
						CHAdjlist[rightNode].insert(make_pair(leftNode, pairs(weight, TDNodeNow)));
					}
					else
					{
						if (CHAdjlist[rightNode][leftNode].first > weight)
						{
							CHAdjlist[rightNode][leftNode].first = weight;
							CHAdjlist[rightNode][leftNode].second = TDNodeNow;
						}
					}
				}
			}

			for (int i = 0; i < nei.size(); i++)
			{
				int neiID = nei[i];
				int height = max(Graph_D_H::CHInfo[TDNodeNow].second + 1, Graph_D_H::CHInfo[neiID].second);
				if (height == Graph_D_H::CHInfo[neiID].second) {
					//CHChild[neiID]++;
					//DBC[neiID] = DBC[neiID] * (CHChild[TDNodeNow] + 1);
					CHChild[neiID] += CHChild[TDNodeNow] + 1;
				}
				else {
					CHChild[neiID] = 1 + CHChild[TDNodeNow];
					//DBC[neiID] =1 + CHChild[TDNodeNow];
				}
				Graph_D_H::CHInfo[neiID].pairsReset(CHAdjlist[neiID].size(), height);
				//if(!visited[neiID])
				heap.insert(HostRank(neiID));
			}
			//std::cout << "; add adj , heapSize: " << heap.size() << ""<<endl;
		}

		int head = nonAdjcentNode[nonAdjcentNode.size() - 1].second;
		std::cout << "head: " << head << "and  CH adjlist size: " << CHAdjlist[head].size() << endl;

		for (auto& it : CHAdjlist) {
			size1 += it.size() + 1;
		}


		std::cout << "vertexes constructed in CPU: " << total << endl;
		DecompositionTreeSize = ((double)(size1 * sizeof(Graph_D_H::myPair<int>))) / (1024 * 1024);
		std::cout << "\t actual Decomposition Tree size = " << ((double)(size1 * sizeof(Graph_D_H::myPair<int>))) / (1024 * 1024) << "MB" << endl;
		std::cout << endl;
		//set<int> mai;
		//for (auto& i : nonAdjcentNode)
		//{
		//	mai.insert(i.second);
		//}
		//std::cout << "Number:" << NodeNumber << "  Actual : " << mai.size() << endl;
		//int size = 0;
		//for (int i = 0; i < NodeNumber; i++)
		//{
		//	if (!visited[i]) {
		//		size++;
		//	}
		//}
		//cout <<"left"<<  size << endl;
	}


	void Graph_D_H::Graph::partition_TD_T(int height)
	{
		CHTree.clear();
		CHTree.shrink_to_fit();
		CHTreeHash.clear();
		CHTreeHash.shrink_to_fit();
		CHTree_D.clear();
		CHTree_D.shrink_to_fit();
		CHTreeHash_D.clear();
		CHTreeHash_D.shrink_to_fit();

		CHInfo.assign(NodeNumber, pairs());
		CHChild.assign(NodeNumber, 0);
		DBC.assign(NodeNumber, 1);
		//loading CHAdjList
		for (int i = 0; i < NodeNumber; i++) {
			for (auto& it : adjList[i]) {
				CHAdjlist[i].emplace(it.first, pairs(it.second, -1));
			}
			CHInfo[i].pairsReset(adjList[i].size(), 0);
			adjList[i].clear();
			adjList[i].shrink_to_fit();
		}
		adjList.clear();
		adjList.shrink_to_fit();
		int tempHeight = TreeHeight;
		cout << "start TD" << endl;
		while (tempHeight > height) {
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 1;

			int totalSubGraph = lowestIndexEnd - lowestIndexStart;
			int maxSubGraphPerThread = totalSubGraph / threadNumber + 1;
			vector<vector<HostRank>> candidateHeap(threadNumber, vector<HostRank>());
			vector<int> range(threadNumber, INT_MAX);
			int it = 0;
			int temp = 0;
			//allocate contract vertexes on different thread
			for (int i = lowestIndexStart; i < lowestIndexEnd; i++) {
				int left = NANHash[i].first;
				int right = NANHash[i].NodeID;

				for (int j = left; j < right + left; j++) {
					candidateHeap[it].push_back(HostRank(nonAdjcentNode[j].second));
				}
				range[it] = min(range[it], left);
				temp++;
				if (temp == maxSubGraphPerThread) {
					it++;
					temp = 0;
				}
			}
			std::vector<std::thread> threads;
			for (int i = 0; i < threadNumber; i++) {
				threads.emplace_back(
					[this, &candidateHeap, &range, i]() {

						while (!candidateHeap[i].empty()) {
							thrust::stable_sort(candidateHeap[i].begin(), candidateHeap[i].end());
							int TDNodeNow = candidateHeap[i].begin()->x;
							vector<pairs> info;
							vector<int> nei;
							nonAdjcentNode[range[i]].setDetal(CHAdjlist[TDNodeNow].size(), TDNodeNow, Graph_D_H::CHInfo[TDNodeNow].second);
							range[i]++;
							for (auto it = CHAdjlist[TDNodeNow].begin(); it != CHAdjlist[TDNodeNow].end(); it++) {
								nei.push_back(it->first);
								info.push_back(it->second);
							}
							for (int its = 0; its < nei.size(); its++) {
								int neiID = nei[its];
								if (CHAdjlist[neiID].find(TDNodeNow) != CHAdjlist[neiID].end()) {
									CHAdjlist[neiID].erase(TDNodeNow);
								}
							}
							for (int it = 0; it < nei.size(); it++) {
								int leftNode = nei[it];
								for (int jt = it + 1; jt < nei.size(); jt++) {
									int rightNode = nei[jt];
									int weight = info[it].first + info[jt].first;
									if (CHAdjlist[leftNode].find(rightNode) == CHAdjlist[leftNode].end()) {
										CHAdjlist[leftNode].insert(make_pair(rightNode, pairs(weight, TDNodeNow)));
									}
									else {
										if (CHAdjlist[leftNode][rightNode].first > weight) {
											CHAdjlist[leftNode][rightNode].first = weight;
											CHAdjlist[leftNode][rightNode].second = TDNodeNow;
										}
									}
									if (CHAdjlist[rightNode].find(leftNode) == CHAdjlist[rightNode].end()) {
										CHAdjlist[rightNode].insert(make_pair(leftNode, pairs(weight, TDNodeNow)));
									}
									else {
										if (CHAdjlist[rightNode][leftNode].first > weight) {
											CHAdjlist[rightNode][leftNode].first = weight;
											CHAdjlist[rightNode][leftNode].second = TDNodeNow;
										}
									}
								}
							}
							for (int it = 0; it < nei.size(); it++) {
								int neiID = nei[it];
								int height = max(Graph_D_H::CHInfo[TDNodeNow].second + 1, Graph_D_H::CHInfo[neiID].second);
								if (height == Graph_D_H::CHInfo[neiID].second) {
									//CHChild[neiID]++;
									CHChild[neiID] += CHChild[TDNodeNow];
								}
								else {
									CHChild[neiID] = 1 + CHChild[TDNodeNow];
								}
								Graph_D_H::CHInfo[neiID].pairsReset(CHAdjlist[neiID].size(), height);

							}
							candidateHeap[i].erase(candidateHeap[i].begin());
						}
					}
				);
			}

			for (auto& its : threads) {
				its.join();
			}
			tempHeight--;
		}





		while (tempHeight > 0) {
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 1;
			tempHeight--;
			for (int i = lowestIndexStart; i < lowestIndexEnd; i++) {
				int left = NANHash[i].first;
				int right = NANHash[i].NodeID;
				tempnonAdjcentNodeIndex = min(left, tempnonAdjcentNodeIndex);
				for (int j = left; j < right + left; j++) {
					heap.insert(HostRank(nonAdjcentNode[j].second));
				}
			}
		}

	}


	void Graph_D_H::Graph::TDInCPUMultiThread(int tempheight)
	{
		Graph_D_H::time_Mine time;
		std::cout << "CPU partition contract start: " << endl;
		time.updateStart();
		partition_TD_T(tempheight);

		time.updateEnd();
		translateBackAndMallocToCPUTime = time.get_microsecond_duration();
		std::cout << "\t CPU multithread  partition contract end, using time: " << time.get_microsecond_duration() << endl;

		std::cout << "\n" << "TD in CPU serial start " << endl;
		time.updateStart();
		CHUsingHeap_2();
		time.updateEnd();
		TDInCPUTime = time.get_microsecond_duration() + translateBackAndMallocToCPUTime;
		std::cout << "\t TD in CPU serial end, using time: " << time.get_microsecond_duration() << endl;
	}
}