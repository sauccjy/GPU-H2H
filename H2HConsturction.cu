#pragma once
#include"Graph_D.cuh"

namespace Graph_D_H
{


	void Graph_D_H::Graph::inConstructH2H_MultiThread()
	{
		Graph_D_H::time_Mine time;
		std::cout << "malloc H2H label start " << std::endl;
		time.updateStart();
		mallocH2HLabel();
		time.updateEnd();
		H2HmallocTime = time.get_microsecond_duration();
		std::cout << "\t malloc end, using time: " << time.get_microsecond_duration() << std::endl;

		std::cout << "construct bfs tree start " << endl;
		time.updateStart();
		constructBFSTree();
		time.updateEnd();
		H2HConstructBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs construct end, using time: " << time.get_microsecond_duration() << std::endl;

		makeH2Hlabel_MultiThread();

	}
	void Graph_D_H::Graph::inConstructH2H_D()
	{
		mallocH2HLabel();
		constructBFSTree();
		translateH2HLabel();
		//displayH2HLabel();
		makeH2Hlabel_D();
		translateH2HLabelBack();
		//displayH2HLabel();
	}

	void Graph_D_H::Graph::inConstructH2H_mix()
	{
		mallocH2HLabel();
		makeH2Hlabel_mix();
		translateH2HLabelBack();
	}

	void Graph_D_H::Graph::inConstructH2H_noComm()
	{
		Graph_D_H::time_Mine time;
		std::cout << "malloc H2H label start " << std::endl;
		time.updateStart();
		mallocH2HLabel();
		time.updateEnd();
		H2HmallocTime = time.get_microsecond_duration();
		std::cout << "\t malloc end, using time: " << time.get_microsecond_duration() << std::endl;

		std::cout << "construct bfs tree start " << std::endl;
		time.updateStart();
		constructBFSTree();
		time.updateEnd();
		H2HConstructBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs construct end, using time: " << time.get_microsecond_duration() << std::endl;

		std::cout << "translate bfs tree start " << std::endl;
		time.updateStart();
		TreeBFS_D = TreeBFS;
		time.updateEnd();
		H2HTranslateBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs translate end, using time: " << time.get_microsecond_duration() << std::endl;


		makeH2Hlabel_noComm();


		translateH2HLabelBack();
	}

	void Graph_D_H::Graph::cleanAfterTreeStruct() {
	}

	void Graph_D_H::Graph::cleanBeforeQuery()
	{
		TreeBFS.clear();
		TreeBFS.shrink_to_fit();
		TreeBFS_D.clear();
		TreeBFS_D.shrink_to_fit();
		TreeHash.clear();
		TreeHash.shrink_to_fit();
		father.clear();
		father.shrink_to_fit();
		father_D.clear();
		father_D.shrink_to_fit();
	}

	void Graph_D_H::Graph::mallocH2HLabel() {

		H2H_startIndex.assign(NodeNumber, pairs(-1, 0));//first: label start index , second :label size
		H2H_label.assign(H2H_Size, H2HLabel());
		int temp = 0;
		queue<int > heap;
		heap.push(head);
		if (oulaOnly) {
			while (!heap.empty())
			{
				int nodeID = heap.front();
				heap.pop();
				H2H_startIndex[nodeID].second = OULA_Only[nodeID];
				H2H_startIndex[nodeID].first = temp;
				temp += H2H_startIndex[nodeID].second;
				H2H_label[temp - 1].Node = nodeID;
				H2H_label[temp - 1].isTreeNode = true;
				H2H_label[temp - 1].Hub = -1;
				H2H_label[temp - 1].dis = 0;
				if (nodeID != head)
				{
					for (int i = 0; i < H2H_startIndex[father[nodeID]].second; i++)
					{
						int sameID = H2H_label[i + H2H_startIndex[father[nodeID]].first].Node;
						H2H_label[i + H2H_startIndex[nodeID].first].Node = sameID;

						if (CHAdjlist[nodeID].find(sameID) != CHAdjlist[nodeID].end()) {
							H2H_label[i + H2H_startIndex[nodeID].first].isTreeNode = true;
							H2H_label[i + H2H_startIndex[nodeID].first].dis = CHAdjlist[nodeID][sameID].first;
							H2H_label[i + H2H_startIndex[nodeID].first].Hub = CHAdjlist[nodeID][sameID].second;
						}
					}
				}


				for (int i = ChildHash[nodeID].second; i < ChildHash[nodeID].second + ChildHash[nodeID].first; i++)
				{
					heap.push(Childs[i]);
				}
			}


		}
		else {
			while (!heap.empty())
			{
				int nodeID = heap.front();
				heap.pop();
				H2H_startIndex[nodeID].second = OULA_DFS_[firstAppeare[nodeID]].second;
				H2H_startIndex[nodeID].first = temp;
				temp += H2H_startIndex[nodeID].second;
				H2H_label[temp - 1].Node = nodeID;
				H2H_label[temp - 1].isTreeNode = true;
				H2H_label[temp - 1].Hub = -1;
				H2H_label[temp - 1].dis = 0;
				if (nodeID != head)
				{
					for (int i = 0; i < H2H_startIndex[father[nodeID]].second; i++)
					{
						int sameID = H2H_label[i + H2H_startIndex[father[nodeID]].first].Node;
						H2H_label[i + H2H_startIndex[nodeID].first].Node = sameID;

						if (CHAdjlist[nodeID].find(sameID) != CHAdjlist[nodeID].end()) {
							H2H_label[i + H2H_startIndex[nodeID].first].isTreeNode = true;
							H2H_label[i + H2H_startIndex[nodeID].first].dis = CHAdjlist[nodeID][sameID].first;
							H2H_label[i + H2H_startIndex[nodeID].first].Hub = CHAdjlist[nodeID][sameID].second;
						}
					}
				}
				for (int i = ChildHash[nodeID].second; i < ChildHash[nodeID].second + ChildHash[nodeID].first; i++)
				{
					heap.push(Childs[i]);
				}
			}

		}

		for (auto& it : CHAdjlist)
		{
			it.clear();
		}
		CHAdjlist.clear();

		H2H_startIndex_D = H2H_startIndex;
		H2H_label_D = H2H_label;
	}

	int Graph::findLCA(int x, int y)
	{
		if (x == y)
			return x;
		int firstA = firstAppeare[x];//position  in  OULA_DFS
		int firstB = firstAppeare[y];
		if (firstA > firstB)
			swap(firstA, firstB);
		int k = (int)(log2(firstB - firstA + 1));//k means log2(gap);
		//int k = (int)(log2(abs(firstA - firstB) + 1));
		// min(RMQ[firstA][k], RMQ[firstB - (1 << k) + 1][k]);//(1<<k - 1 ) means davicate index;
		if (RMQ_OneLine[firstA * RMQ_Line_Size + k].first <= RMQ_OneLine[(firstB - (1 << k) + 1) * RMQ_Line_Size + k].first)
			return RMQ_OneLine[firstA * RMQ_Line_Size + k].second;
		else
			return RMQ_OneLine[(firstB - (1 << k) + 1) * RMQ_Line_Size + k].second;
	}

	int Graph_D_H::Graph::findLCA_noHub(int x, int y)
	{
		if (x == y)
			return x;
		int firstA = firstAppeare[x];//position  in  OULA_DFS
		int firstB = firstAppeare[y];
		if (firstA > firstB)
			swap(firstA, firstB);
		int k = (int)(log2(firstB - firstA + 1));//k means log2(gap);

		//int64_t firstARMQPlace = RMQHash[firstA] + k;
		//int64_t firstBRMQPlace = RMQHash[firstB - (1 << k) + 1] + k;

		int ID_1 = RMQ_ID[RMQHash[firstA] + (int64_t)k];
		int ID_2 = RMQ_ID[RMQHash[firstB - (1 << k) + 1] + (int64_t)k];

		if (H2H_TreeHeight_Hash[ID_1] <= H2H_TreeHeight_Hash[ID_2]) {
			return ID_1;
		}
		else {
			return ID_2;
		}

	}

	int Graph::H2HDistancQuery_UsingLCA(int x, int y)
	{
		if (x == y)
			return 0;
		int LCA = findLCA(x, y);
		int checkPointX = H2H_startIndex[x].first;
		int checkPointY = H2H_startIndex[y].first;
		int LCApoint = H2H_startIndex[LCA].first;
		int frontier = H2H_startIndex[LCA].second;
		int result = INT_MAX;
		for (int i = 0; i < frontier; i++)
		{
			if (!H2H_label[LCApoint + i].isTreeNode)
				continue;
			result = ((H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) < result) ?
				(H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) : result;
		}
		return result;
	}

	int Graph::H2HDistancQuery_UsingOULA(int x, int y)
	{
		if (x == y)
			return 0;
		int checkPointX = H2H_startIndex[x].first;
		int checkPointY = H2H_startIndex[y].first;

		int frontier = min(H2H_startIndex[x].second, H2H_startIndex[y].second);
		int result = INT_MAX;
		for (int i = 0; i < frontier; i++)
		{
			if (H2H_label[checkPointX + i].Node != H2H_label[checkPointY + i].Node)
				break;
			result = ((H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) < result) ?
				(H2H_label[checkPointX + i].dis + H2H_label[checkPointY + i].dis) : result;
		}
		return result;
	}

	int Graph::H2HDistancQuery_AncestorToPost(int x, int y)
	{
		//int checkPointX = H2H_startIndex[x].first;
		//int checkPointY = H2H_startIndex[y].first;
		//int minSize = H2H_startIndex[x].second;
		return H2H_label[H2H_startIndex[y].first + H2H_startIndex[x].second - 1].dis;
	}

	void Graph_D_H::Graph::makeH2Hlabel() {
		queue<int> frontier;
		for (int i = ChildHash[head].second; i < ChildHash[head].second + ChildHash[head].first; i++)
		{
			frontier.push(Childs[i]);
		}
		//int height = 1;
		while (!frontier.empty())
		{
			int size = frontier.size();
			//cout << "height: " << height++ << " size : " << size <<" total: "<<total << endl;
			for (int i = 0; i < size; i++)
			{
				int ID = frontier.front();
				frontier.pop();
				int H2HSize = H2H_startIndex[ID].second - 1;
				int IDNowH2Hindex = H2H_startIndex[ID].first;
				//displaySingalLabel(ID);
				//construct others
				for (int k = IDNowH2Hindex + H2HSize - 1; k >= IDNowH2Hindex; k--)
				{
					if (!H2H_label[k].isTreeNode) {
						continue;
					}
					//int adjH2HStartIndex = H2H_startIndex[H2H_label[k].Node].first;
					int q = k - 1;
					for (; q >= IDNowH2Hindex; q--)
					{
						int templength = H2H_label[k].dis + H2H_label[q - IDNowH2Hindex + H2H_startIndex[H2H_label[k].Node].first].dis;
						if (H2H_label[q].dis > templength) {
							H2H_label[q].dis = templength;
							H2H_label[q].Hub = H2H_label[k].Node;
						}
						//H2H_label[q].dis = min(H2H_label[q].dis, H2H_label[k].dis + H2H_label[q - IDNowH2Hindex + H2H_startIndex[H2H_label[k].Node].first].dis);
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

				for (int j = ChildHash[ID].second; j < ChildHash[ID].second + ChildHash[ID].first; j++)
				{
					frontier.push(Childs[j]);
				}
			}
		}
	}


	void Graph_D_H::Graph::constructBFSTree()
	{
		queue<int> frontier;
		int maxLabelWidth = 1;
		int tempHeight = 0;
		int height = 0;
		set<int> vertexInMaxWitdh = {};
		TreeHash.push_back(0);
		TreeBFS.push_back(head);
		TreeHash.push_back(TreeBFS.size());

		TreeBFS_Hash.push_back(0);
		TreeBFS_ID.push_back(head);
		TreeBFS_adj.push_back(head);
		TreeBFS_pos.push_back(0);
		//TreeBFS_changeTime.push_back(0);
		TreeBFS_Hash.push_back(1);

		for (int i = ChildHash[head].second; i < ChildHash[head].second + ChildHash[head].first; i++)
		{
			frontier.push(Childs[i]);
		}
		//int height = 1;
		while (!frontier.empty())
		{
			int size = frontier.size();
			//cout << "height: " << height++ << " size : " << size <<" total: "<<total << endl;
			for (int i = 0; i < size; i++)
			{
				int ID = frontier.front();
				frontier.pop();
				TreeBFS.push_back(ID);
				for (int64_t index = H2H_pos_hash[ID]; index < H2H_pos_hash[ID + 1]; index++) {
					TreeBFS_ID.push_back(ID);
					TreeBFS_adj.push_back(H2H_pos_ID[index]);
					TreeBFS_pos.push_back(H2H_pos_POS[index]);
					//TreeBFS_changeTime.push_back(H2H_pos_hash[ID + 1] - 2 - index);
				}
				for (int j = ChildHash[ID].second; j < ChildHash[ID].second + ChildHash[ID].first; j++)
				{
					frontier.push(Childs[j]);
				}
			}
			TreeHash.push_back(TreeBFS.size());
			TreeBFS_Hash.push_back(TreeBFS_ID.size());
		}


		Childs.clear();
		Childs.shrink_to_fit();
		ChildHash.clear();
		ChildHash.shrink_to_fit();
		//TreeBFS_D = TreeBFS;
	}









	void Graph_D_H::Graph::makeH2Hlabel_MultiThread()
	{
		Graph_D_H::time_Mine time;

		std::cout << "construct H2H on CPU start " << std::endl;
		time.updateStart();

		for (int tempHeight = 0; tempHeight < TreeHash.size() - 1; tempHeight++)
		{
			int startIndex = TreeHash[tempHeight];
			int endIndex = TreeHash[tempHeight + 1];

			int tempThreadNum = min(threadNumber, endIndex - startIndex);

			int tempSize = (endIndex - startIndex) / tempThreadNum + 1;
			vector<vector<int>> candidateHeap(tempThreadNum, vector<int>());
			int it = 0, temp = 0;
			for (int i = startIndex; i < endIndex; i++) {
				int ID = TreeBFS[i];
				candidateHeap[it].emplace_back(ID);
				temp++;
				if (temp == tempSize) {
					it++;
					temp = 0;
				}
			}
			std::vector<std::thread> threads;

			for (int i = 0; i < tempThreadNum; i++) {

				threads.emplace_back(
					[this, &candidateHeap, i]() {
						for (auto& ID : candidateHeap[i]) {
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
					}
				);

			}

			for (auto& its : threads) {
				its.join();
			}

		}
		translateBeforeQuery();
		time.updateEnd();
		H2HUsingTime_CPU = time.get_microsecond_duration();
		std::cout << "\t CPU construct end, using time: " << time.get_microsecond_duration() << std::endl;
	}



	void Graph_D_H::Graph::inConstructH2H_noHub()
	{
		//malloc
		//mallocH2HLabel_noHub();

		Graph_D_H::time_Mine time;
		std::cout << "malloc H2H label start " << std::endl;
		time.updateStart();
		mallocH2HLabel_noHub();
		time.updateEnd();
		H2HmallocTime = time.get_microsecond_duration();
		std::cout << "\t malloc end, using time: " << time.get_microsecond_duration() << std::endl;

		//constructBFSTree
		std::cout << "construct bfs tree start " << endl;
		time.updateStart();
		constructBFSTree();
		time.updateEnd();
		H2HConstructBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs construct end, using time: " << time.get_microsecond_duration() << std::endl;
		//make
		makeH2HLabel_noHub_serial();
	}

	void Graph_D_H::Graph::inConstructH2H_noHub_multiThread()
	{
		//malloc
		Graph_D_H::time_Mine time;
		std::cout << "malloc H2H label start " << std::endl;
		time.updateStart();
		mallocH2HLabel_noHub();
		time.updateEnd();
		H2HmallocTime = time.get_microsecond_duration();
		std::cout << "\t malloc end, using time: " << time.get_microsecond_duration() << std::endl;

		//constructBFSTree
		std::cout << "construct bfs tree start " << endl;
		time.updateStart();
		constructBFSTree();
		time.updateEnd();
		H2HConstructBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs construct end, using time: " << time.get_microsecond_duration() << std::endl;
		//make
		makeH2HLabel_noHub_multiThred();
	}

	void Graph_D_H::Graph::inConstructH2H_noHub_D()
	{
		//malloc
		Graph_D_H::time_Mine time;
		std::cout << "malloc H2H label start " << std::endl;
		time.updateStart();
		mallocH2HLabel_noHub();
		time.updateEnd();
		H2HmallocTime = time.get_microsecond_duration();
		std::cout << "\t malloc end, using time: " << time.get_microsecond_duration() << std::endl;


		//constructBFSTree
		std::cout << "construct bfs tree start " << endl;
		time.updateStart();
		constructBFSTree();
		time.updateEnd();
		H2HConstructBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs construct end, using time: " << time.get_microsecond_duration() << std::endl;

		std::cout << "translate bfs tree start " << std::endl;
		time.updateStart();
		frontier_Hash_D = frontier_Hash;
		frontier_ID_D = frontier_ID;
		TreeBFS_D = TreeBFS;
		//TreeBFS_changeTime_D = TreeBFS_changeTime;
		time.updateEnd();
		H2HTranslateBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs translate end, using time: " << time.get_microsecond_duration() << std::endl;
		//make
		//displayH2H_noHub();

		//makeH2HLabel_noHub_noComm();
		//makeH2HLabel_noHub_noComm_2();
		makeH2HLabel_noHub_noComm_3();
		H2H_dis = H2H_dis_D;

		//H2H_pos_ID_D.clear();
		//H2H_pos_ID_D.shrink_to_fit();
		//H2H_pos_ID.clear();
		//H2H_pos_ID.shrink_to_fit();

		frontier_Hash.clear();
		frontier_Hash.shrink_to_fit();
		frontier_Hash_D.clear();
		frontier_Hash_D.shrink_to_fit();
		frontier_ID.clear();
		frontier_ID.shrink_to_fit();
		frontier_ID_D.clear();
		frontier_ID_D.shrink_to_fit();
	}

	void Graph_D_H::Graph::inConstructH2H_noHub_D_2()
	{
		//malloc
		Graph_D_H::time_Mine time;
		std::cout << "malloc H2H label start " << std::endl;
		time.updateStart();
		mallocH2HLabel_noHub();
		time.updateEnd();
		H2HmallocTime = time.get_microsecond_duration();
		std::cout << "\t malloc end, using time: " << time.get_microsecond_duration() << std::endl;


		//constructBFSTree
		std::cout << "construct bfs tree start " << endl;
		time.updateStart();
		constructBFSTree();
		time.updateEnd();
		H2HConstructBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs construct end, using time: " << time.get_microsecond_duration() << std::endl;

		std::cout << "translate bfs tree start " << std::endl;
		time.updateStart();
		frontier_Hash_D = frontier_Hash;
		frontier_ID_D = frontier_ID;
		TreeBFS_D = TreeBFS;
		//TreeBFS_changeTime_D = TreeBFS_changeTime;
		time.updateEnd();
		H2HTranslateBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs translate end, using time: " << time.get_microsecond_duration() << std::endl;
		//make
		//displayH2H_noHub();

		makeH2HLabel_noHub_noComm_2();
		//makeH2HLabel_noHub_noComm_3();
		H2H_dis = H2H_dis_D;

		//H2H_pos_ID_D.clear();
		//H2H_pos_ID_D.shrink_to_fit();
		//H2H_pos_ID.clear();
		//H2H_pos_ID.shrink_to_fit();

		frontier_Hash.clear();
		frontier_Hash.shrink_to_fit();
		frontier_Hash_D.clear();
		frontier_Hash_D.shrink_to_fit();
		frontier_ID.clear();
		frontier_ID.shrink_to_fit();
		frontier_ID_D.clear();
		frontier_ID_D.shrink_to_fit();
	}


	void Graph_D_H::Graph::mallocH2HLabel_noHub()
	{
		int totalSize = 0;
		for (auto& it : CHAdjlist) {
			totalSize += it.size() + 1;
		}
		cout << " \t pos size: " << totalSize << endl;
		cout << " \t dis size: " << H2H_Size << endl;
		//H2H_pos_hash.assign(NodeNumber + 1, -1);
		//H2H_dis_hash.assign(NodeNumber + 1, -1);
		//H2H_dis.assign(H2H_Size, INT_MAX);
		//H2H_pos_ID.assign(totalSize, -1);
		//H2H_pos_POS.assign(totalSize, -1);
		cout << "\t CPU allocate start! " << endl;

		try{
		H2H_pos_hash.push_back(0);
		H2H_dis_hash.push_back(0);
		for (int i = 0; i < NodeNumber; i++) {
			H2H_dis.insert(H2H_dis.end(), OULA_DFS_[firstAppeare[i]].second, INT_MAX);

			for (auto& adj : CHAdjlist[i]) {
				int adjID = adj.first;
				int adjPos = OULA_DFS_[firstAppeare[adjID]].second - 1;
				H2H_pos_ID.push_back(adjID);
				H2H_pos_POS.push_back(adjPos);
				H2H_dis[H2H_dis_hash[i] + (int64_t)adjPos] = adj.second.first;
			}
			H2H_pos_ID.push_back(i);
			H2H_pos_POS.push_back(OULA_DFS_[firstAppeare[i]].second - 1);
			H2H_dis[H2H_dis.size() - 1] = 0;
			H2H_pos_hash.push_back(H2H_pos_ID.size());
			H2H_dis_hash.push_back(H2H_dis.size());

			thrust::sort_by_key(H2H_pos_POS.begin() + H2H_pos_hash[i], H2H_pos_POS.end(),
				H2H_pos_ID.begin() + H2H_pos_hash[i]);
		}
		}
		catch(const std::bad_alloc& e){
				std::cerr << "bad allocate : " << e.what() <<" at RAM-label-allocation"<< '\n';  
				std::exit(1);
		}
		//int64_t posTemp = 0;
		//int64_t distemp = 0;
		//for (auto nodeID = 0; nodeID < NodeNumber; nodeID++) {
		//	H2H_pos_hash[nodeID] = posTemp;
		//	H2H_dis_hash[nodeID] = distemp;
		//	//H2H_dis.insert(H2H_dis.end(), OULA_DFS_[firstAppeare[nodeID]].second, INT_MAX);



		//	for (auto& adj : CHAdjlist[nodeID]) {
		//		int adjID = adj.first;
		//		int adjPos = OULA_DFS_[firstAppeare[adjID]].second - 1;
		//		H2H_pos_POS[posTemp] = adjPos;
		//		H2H_pos_ID[posTemp++] = adjID;
		//		//H2H_pos_ID.push_back(adjID);
		//		//H2H_pos_POS.push_back(adjPos);

		//		H2H_dis[H2H_dis_hash[nodeID] + adjPos] = adj.second.first;
		//	}

		//	//posTemp += CHAdjlist[nodeID].size();
		//	distemp += OULA_DFS_[firstAppeare[nodeID]].second;
		//	H2H_pos_POS[posTemp] = nodeID;
		//	H2H_pos_ID[posTemp++] = OULA_DFS_[firstAppeare[nodeID]].second - 1;
		//	H2H_dis[H2H_dis_hash[nodeID] + OULA_DFS_[firstAppeare[nodeID]].second - 1] = 0;

		//	thrust::sort_by_key(H2H_pos_POS.begin() + H2H_pos_hash[nodeID], H2H_pos_POS.begin() + posTemp,
		//		H2H_pos_ID.begin() + H2H_pos_hash[nodeID]);
		//}


		//H2H_pos_hash[NodeNumber] = H2H_pos_ID.size();
		//H2H_dis_hash[NodeNumber] = H2H_dis.size();
		cout << "\t CPU allocate end" << endl;
		//translate
		//H2H_pos_hash_D = H2H_pos_hash;
		//H2H_dis_hash_D = H2H_dis_hash;
		//H2H_pos_ID_D = H2H_pos_ID;
		//H2H_pos_POS_D = H2H_pos_POS;
		//H2H_dis_D = H2H_dis;

		double H2HLabelSize_T;
		cout << "\t CHTreeHeight: " << CHTreeHeight << endl;
		TDTREEHeight = CHTreeHeight;
		double rmq = RMQHash.size() * sizeof(int64_t) + RMQ_ID.size() * sizeof(int) + H2H_TreeHeight_Hash.size() * sizeof(int);
		RMQsize = rmq / (1024 * 1024);
		H2HLabelSize = ((double)(H2H_pos_hash.size() * sizeof(int64_t) + H2H_dis_hash.size() * sizeof(int64_t)  +
			H2H_pos_POS.size() * sizeof(int) + H2H_dis.size() * sizeof(int))) / (1024 * 1024);
		H2HLabelSize_T = (H2H_pos_POS.size() * sizeof(int) + H2H_dis.size() * sizeof(int)) / (1024 * 1024);
		cout << "\t RMQ size : " << RMQsize
			<< "MB, H2H Label size: " << H2HLabelSize_T << "MB" << endl;

		OULA_DFS_.clear();
		OULA_DFS_.shrink_to_fit();
		for (auto& it : CHAdjlist) {
			it.clear();
		}
		CHAdjlist.clear();
		CHAdjlist.shrink_to_fit();

		//displayH2H_noHub();

	}




	void Graph_D_H::Graph::makeH2HLabel_noHub_serial()
	{
		Graph_D_H::time_Mine time;

		std::cout << "construct H2H on CPU-serial start " << std::endl;
		time.updateStart();
		//for (int tempHeight = 0; tempHeight < TreeHash.size() - 1; tempHeight++)
		//{
		//	int startIndex = TreeHash[tempHeight];
		//	int endIndex = TreeHash[tempHeight + 1];
		//	for (int i = startIndex; i < endIndex; i++) {
		//		int nodeID = TreeBFS[i];

		//		for (size_t temp = H2H_pos_hash[(size_t)nodeID + 1] - 1; temp >= H2H_pos_hash[(size_t)nodeID]; temp--) {
		//			int adjID = H2H_pos_ID[temp];
		//			int pos = H2H_pos_POS[temp];
		//			size_t indexNode = H2H_dis_hash[nodeID];
		//			size_t index1_adj = H2H_dis_hash[adjID];


		//			int tempLength_ = H2H_dis[indexNode + (size_t)pos];
		//			for (size_t j = (size_t)pos - 1; j >= 0; j--) {
		//				H2H_dis[indexNode + j] = min(H2H_dis[indexNode + j],
		//					tempLength_ + H2H_dis[index1_adj + j]);
		//			}

		//			//L4 label
		//			int fatherID = nodeID;
		//			for (size_t j = H2H_dis_hash[(size_t)nodeID + 1] - H2H_dis_hash[(size_t)nodeID] - 2; j > pos; j--) {
		//				fatherID = father[fatherID];
		//				
		//				H2H_dis[indexNode + j] = min(H2H_dis[indexNode + j],
		//					tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + (size_t)pos]);
		//			}
		//		}
		//	}
		//}

		for (int tempHeight = 0; tempHeight < TreeHash.size() - 1; tempHeight++) {
			int startIndex = TreeHash[tempHeight];
			int endIndex = TreeHash[tempHeight + 1];
			//cout << "at height: " << tempHeight << endl;
			for (int i = startIndex; i < endIndex; i++) {
				int nodeID = TreeBFS[i];
				//cout << "\tNodeID: " << nodeID <<"start at: "<< H2H_pos_hash[(size_t)nodeID] << "size: " << H2H_pos_hash[(size_t)nodeID + 1] - H2H_pos_hash[(size_t)nodeID] << endl;
				for (int64_t temp = H2H_pos_hash[(int64_t)nodeID + 1] - 2; temp >= H2H_pos_hash[(int64_t)nodeID]; temp--) {
					int adjID = H2H_pos_ID[temp];
					int pos = H2H_pos_POS[temp];
					int64_t indexNode = H2H_dis_hash[nodeID];
					int64_t index1_adj = H2H_dis_hash[adjID];
					int tempLength_ = H2H_dis[indexNode + (int64_t)pos];
					for (int j = pos - 1; j >= 0; j--) {
						H2H_dis[indexNode + (int64_t)j] = std::min(H2H_dis[indexNode + (int64_t)j],
							tempLength_ + H2H_dis[index1_adj + (int64_t)j]);
					}
					// L4 label
					//int64_t fatherID = nodeID;
					//for (int j = (int)(H2H_dis_hash[nodeID + 1] - H2H_dis_hash[nodeID]) - 2; j > pos; j--) {
					//	fatherID = (int64_t)father[fatherID];

					//	H2H_dis[indexNode + (int64_t)j] = std::min(H2H_dis[indexNode + (int64_t)j],
					//		tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + (int64_t)pos]);
					//}
				}
			}
		}
		cout << "finish calculate" << endl;
		//translateBeforeQuery();
		translateH2H_noHub();

		time.updateEnd();
		H2HUsingTime_CPU = time.get_microsecond_duration();
		std::cout << "\t CPU construct end, using time: " << time.get_microsecond_duration() << std::endl;
	}





	void Graph::makeH2HLabel_noHub_multiThred()
	{
		Graph_D_H::time_Mine time;

		std::cout << "construct H2H on CPU multi-thread start " << std::endl;
		time.updateStart();

		for (int64_t tempHeight = 0; tempHeight < (int64_t)TreeHash.size() - 1; tempHeight++)
		{
			int64_t startIndex = TreeHash[tempHeight];
			int64_t endIndex = TreeHash[tempHeight + 1];
			Graph_D_H::time_Mine time1;
			time1.updateStart();
			int64_t tempThreadNum = min((int64_t)threadNumber, endIndex - startIndex);

			int64_t tempSize = (endIndex - startIndex) / tempThreadNum + 1;
			vector<vector<int64_t>> candidateHeap(tempThreadNum, vector<int64_t>());
			int64_t it = 0, temp = 0;
			for (int64_t i = startIndex; i < endIndex; i++) {
				int64_t ID = TreeBFS[i];
				candidateHeap[it].emplace_back(ID);
				temp++;
				if (temp == tempSize) {
					it++;
					temp = 0;
				}
			}
			std::vector<std::thread> threads;

			for (int i = 0; i < tempThreadNum; i++) {

				threads.emplace_back(
					[this, &candidateHeap, i]() {
						for (auto& nodeID : candidateHeap[i]) {
							for (int64_t temp = H2H_pos_hash[(int64_t)nodeID + 1] - 2; temp >= H2H_pos_hash[(int64_t)nodeID]; temp--) {
								int adjID = H2H_pos_ID[temp];
								int pos = H2H_pos_POS[temp];
								int64_t indexNode = H2H_dis_hash[nodeID];
								int64_t index1_adj = H2H_dis_hash[adjID];
								int tempLength_ = H2H_dis[indexNode + (int64_t)pos];
								for (int j = pos - 1; j >= 0; j--) {
									H2H_dis[indexNode + (int64_t)j] = std::min(H2H_dis[indexNode + (int64_t)j],
										tempLength_ + H2H_dis[index1_adj + (int64_t)j]);
								}
								// L4 label
								//int64_t fatherID = nodeID;
								//for (int j = (int)(H2H_dis_hash[nodeID + 1] - H2H_dis_hash[nodeID]) - 2; j > pos; j--) {
								//	fatherID = (int64_t)father[fatherID];

								//	H2H_dis[indexNode + (int64_t)j] = std::min(H2H_dis[indexNode + (int64_t)j],
								//		tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + (int64_t)pos]);
								//}
							}
						}
					}
				);

			}

			for (auto& its : threads) {
				its.join();
			}
			time1.updateEnd();
			//std::cout << "At height: " << tempHeight << " Kernel execution time: " << time1.get_microsecond_duration() / 1000 << " ms" << std::endl;

		}
		//translateBeforeQuery();
		translateH2H_noHub();

		time.updateEnd();
		H2HUsingTime_CPU = time.get_microsecond_duration();
		std::cout << "\t CPU construct end, using time: " << time.get_microsecond_duration() << std::endl;
	}


	void Graph::makeH2HLabel_noHub_multiThred_2()
	{
		Graph_D_H::time_Mine time;

		std::cout << "construct H2H on CPU multi-thread start " << std::endl;
		time.updateStart();

		for (int64_t tempHeight = 0; tempHeight < (int64_t)TreeBFS_Hash.size() - 1; tempHeight++)
		{
			int64_t startIndex = TreeBFS_Hash[tempHeight];
			int64_t endIndex = TreeBFS_Hash[tempHeight + 1];
			Graph_D_H::time_Mine time1;
			time1.updateStart();
			int64_t tempThreadNum = min((int64_t)threadNumber, endIndex - startIndex);

			int64_t tempSize = (endIndex - startIndex) / tempThreadNum + 1;
			vector<vector<int64_t>> candidateHeap(tempThreadNum, vector<int64_t>());
			int64_t it = 0, temp = 0;
			for (int64_t i = startIndex; i < endIndex; i++) {
				//int64_t ID = TreeBFS[i];
				candidateHeap[it].emplace_back(i);
				temp++;
				if (temp == tempSize) {
					it++;
					temp = 0;
				}
			}
			std::vector<std::thread> threads;

			for (int i = 0; i < tempThreadNum; i++) {

				threads.emplace_back(
					[this, &candidateHeap, i]() {
						for (auto& nodeIDindex : candidateHeap[i]) {

							std::mutex mt;
							int nodeID = TreeBFS_ID[nodeIDindex];
							int adjID = TreeBFS_adj[nodeIDindex];
							int pos = TreeBFS_pos[nodeIDindex];

							int64_t indexNode = H2H_dis_hash[nodeID]; //current vertex's label position
							int64_t index1_adj = H2H_dis_hash[adjID]; //ancestor's label position

							int tempLength_ = H2H_dis[indexNode + pos];
							for (int64_t j = pos - 1; j >= 0; j--) {

								std::lock_guard<std::mutex> lock(mt);
								H2H_dis[indexNode + j] = std::min(H2H_dis[indexNode + j], (tempLength_ + H2H_dis[index1_adj + j]));
								//atomicMin(&H2H_dis[indexNode + j], (tempLength_ + H2H_dis[index1_adj + j]));
							}

							//L4_label
							//int fatherID = nodeID;
							//for (int64_t j = (int)(H2H_dis_hash[nodeID + 1] - H2H_dis_hash[nodeID]) - 2; j > pos; j--) {
							//	fatherID = father[fatherID];
							//	//atomicMin(&H2H_dis[indexNode + j], (tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos]));
							//	std::lock_guard<std::mutex> lock(mt);
							//	H2H_dis[indexNode + j] = std::min(H2H_dis[indexNode + j], (tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + pos]));
							//}

						}
					}
				);
			}

			for (auto& its : threads) {
				its.join();
			}
			time1.updateEnd();
			//std::cout << "At height: " << tempHeight << " Kernel execution time: " << time1.get_microsecond_duration() / 1000 << " ms" << std::endl;

		}
		//translateBeforeQuery();
		translateH2H_noHub();

		time.updateEnd();
		H2HUsingTime_CPU = time.get_microsecond_duration();
		std::cout << "\t CPU construct end, using time: " << time.get_microsecond_duration() << std::endl;
	}

	void Graph_D_H::Graph::inConstructH2H_noHub_multiThread_2()
	{
		//malloc
		Graph_D_H::time_Mine time;
		std::cout << "malloc H2H label start " << std::endl;
		time.updateStart();
		mallocH2HLabel_noHub();
		time.updateEnd();
		H2HmallocTime = time.get_microsecond_duration();
		std::cout << "\t malloc end, using time: " << time.get_microsecond_duration() << std::endl;

		//constructBFSTree
		std::cout << "construct bfs tree start " << endl;
		time.updateStart();
		constructBFSTree();
		time.updateEnd();
		H2HConstructBFSTreeTime = time.get_microsecond_duration();
		std::cout << "\t bfs construct end, using time: " << time.get_microsecond_duration() << std::endl;
		//make
		makeH2HLabel_noHub_multiThred_2();
	}

}
