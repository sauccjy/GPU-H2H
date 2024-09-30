#pragma once
#include"Graph_D.cuh"
#include "kernel_functions.cuh"

namespace Graph_D_H
{
	void Graph_D_H::Graph::startGPU()
	{
		tests();
	}

	void Graph_D_H::Graph::translateBeforeCH()
	{
		ID_hash_D = ID_hash;
		//ranks_D = ranks;
		//partition_Tree_D = partition_Tree;
		nonAdjcentNode_D = nonAdjcentNode;
		NANHash_D = NANHash;
		//father_D = father;
		CHTreeHash_D = CHTreeHash;
		CHTree_D = CHTree;
		CHInfo_D = CHInfo_H;
		visited_D = visited;
	}

	void Graph_D_H::Graph::partition_TD_Parallel_D(int height)
	{

		int heightTarget = height;
		makeCHPerHeight(TreeHeight, cudaThreadNum, heightTarget, nonAdjcentNode_D, ID_hash_D, CHTreeHash_D,
			CHTree_D, NANHash_D,visited_D);

	}

	void Graph_D_H::Graph::partition_TD_Parallel_D_2(int height)
	{

		int heightTarget = height;
		makeCHPerHeight_2(TreeHeight, cudaThreadNum, heightTarget, nonAdjcentNode_D, ID_hash_D, CHTreeHash_D,
			CHTree_D, NANHash_D, visited_D, CHInfo_D);

	}
	void Graph_D_H::Graph::translateInTD()
	{
		CHTreeHash = CHTreeHash_D;
		CHTree = CHTree_D;
		nonAdjcentNode = nonAdjcentNode_D;
		visited = visited_D;
		CHInfo_H = CHInfo_D;
	}
	void Graph::translateCHTree()
	{
		//father = father_D;
		CHTreeHash = CHTreeHash_D;
		CHTree = CHTree_D;
		nonAdjcentNode = nonAdjcentNode_D;
	}


	void Graph::translateOULARMQ()
	{
		firstAppeare_D = firstAppeare;
		RMQ_OneLine_D = RMQ_OneLine;
		//OULA_DFS_D = OULA_DFS_;
	}

	void Graph_D_H::Graph::translateRMQ_noHub()
	{
		firstAppeare_D = firstAppeare;
		H2H_TreeHeight_Hash_D = H2H_TreeHeight_Hash;
		RMQHash_D = RMQHash;
		//RMQ_Height_D = RMQ_Height;
		RMQ_ID_D = RMQ_ID;
	}

	void Graph::translateH2HLabel() {

		H2H_startIndex_D = H2H_startIndex;
		H2H_label_D = H2H_label;
		TreeBFS_D = TreeBFS;
	}


	void Graph::translateH2HLabelBack() {

		H2H_startIndex = H2H_startIndex_D;
		H2H_label = H2H_label_D;
	}

	void Graph_D_H::Graph::translateQuery()
	{
		x_D = x;
		y_D = y;
		result_D = result;
		//translateOULARMQ();
	}

	void Graph_D_H::Graph::translateQueryResultBack()
	{
		result = result_D;
	}

	void Graph_D_H::Graph::makeH2Hlabel_D() {
		thrust::host_vector<int> frontier = {};
		frontier.push_back(head);
		makeH2HLabel(1, cudaThreadNum, H2H_startIndex_D, H2H_label_D, Childs, ChildHash, frontier);
	}


	void Graph_D_H::Graph::makeH2Hlabel_mix() {
		queue<int> frontier;
		for (int i = ChildHash[head].second; i < ChildHash[head].second + ChildHash[head].first; i++)
		{
			frontier.push(Childs[i]);
		}
		int height = 1;
		while (!frontier.empty())
		{
			int size = frontier.size();
			if (size > 128)
				break;
			height++;
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
						if (H2H_label[q].isTreeNode) {
							//k = q + 1;
							continue;
						}
						int templength = H2H_label[k].dis + H2H_label[q - IDNowH2Hindex + H2H_startIndex[H2H_label[k].Node].first].dis;
						if (H2H_label[q].dis > templength) {
							H2H_label[q].dis = templength;
							H2H_label[q].Hub = H2H_label[k].Node;
						}
						//H2H_label[q].dis = min(H2H_label[q].dis, H2H_label[k].dis + H2H_label[q - IDNowH2Hindex + H2H_startIndex[H2H_label[k].Node].first].dis);
					}
				}

				for (int j = ChildHash[ID].second; j < ChildHash[ID].second + ChildHash[ID].first; j++)
				{
					frontier.push(Childs[j]);
				}
			}
		}

		thrust::host_vector<int> frontier1 = {};

		while (!frontier.empty()) {
			frontier1.push_back(frontier.front());
			frontier.pop();
		}
		H2H_startIndex_D = H2H_startIndex;
		H2H_label_D = H2H_label;
		makeH2HLabel(height, cudaThreadNum, H2H_startIndex_D, H2H_label_D, Childs, ChildHash, frontier1);
	}

	void Graph_D_H::Graph::makeH2Hlabel_noComm()
	{
		Graph_D_H::time_Mine time;

		//cout << "construct H2H on CPU start " << endl;
		//time.updateStart();

		int tempHeight = 0;
		//for (; tempHeight < TreeHash.size() - 1; tempHeight++)
		//{
		//	int startIndex = TreeHash[tempHeight];
		//	int endIndex = TreeHash[tempHeight + 1];
		//	if (endIndex - startIndex > 256) {
		//		break;
		//	}
		//	int tempThreadNum = min(threadNumber, endIndex - startIndex);

		//	int tempSize = (endIndex - startIndex) / tempThreadNum + 1;
		//	vector<vector<int>> candidateHeap(tempThreadNum, vector<int>());
		//	int it = 0, temp = 0;
		//	for (int i = startIndex; i < endIndex; i++) {
		//		int ID = TreeBFS[i];
		//		candidateHeap[it].emplace_back(ID);
		//		temp++;
		//		if (temp == tempSize) {
		//			it++;
		//			temp = 0;
		//		}
		//	}
		//	std::vector<std::thread> threads;

		//	for (int i = 0; i < tempThreadNum; i++) {

		//		threads.emplace_back(
		//			[this, &candidateHeap, i]() {
		//				for (auto& ID : candidateHeap[i]) {
		//					int H2HSize = H2H_startIndex[ID].second - 1;
		//					int IDNowH2Hindex = H2H_startIndex[ID].first;

		//					for (int k = IDNowH2Hindex + H2HSize - 1; k >= IDNowH2Hindex; k--)
		//					{
		//						if (!H2H_label[k].isTreeNode) {
		//							continue;
		//						}

		//						for (int q = k - 1; q >= IDNowH2Hindex; q--)
		//						{
		//							//if (H2H_label[q].isTreeNode) {
		//							//	continue;
		//							//}
		//							int templength = H2H_label[k].dis + H2H_label[q - IDNowH2Hindex + H2H_startIndex[H2H_label[k].Node].first].dis;
		//							if (H2H_label[q].dis > templength) {
		//								H2H_label[q].dis = templength;
		//								H2H_label[q].Hub = H2H_label[k].Node;
		//							}
		//						}
		//						for (int q = k + 1; q <= IDNowH2Hindex + H2HSize - 1; q++)
		//						{
		//							int templength = H2H_label[k].dis + H2H_label[k - IDNowH2Hindex + H2H_startIndex[H2H_label[q].Node].first].dis;
		//							if (H2H_label[q].dis > templength) {
		//								H2H_label[q].dis = templength;
		//								H2H_label[q].Hub = H2H_label[k].Node;
		//							}
		//						}
		//					}
		//				}
		//			}
		//		);

		//	}

		//	for (auto& its : threads) {
		//		its.join();
		//	}

		//}
		//time.updateEnd();
		//H2HUsingTime_CPU = time.get_microsecond_duration();
		//cout << "\t CPU construct end, using time: " << time.get_microsecond_duration() << endl;


		//cout << "translate H2H start " << endl;
		//time.updateStart();
		//H2H_startIndex_D = H2H_startIndex;
		//H2H_label_D = H2H_label;
		//time.updateEnd();
		//H2HTranslateTime = time.get_microsecond_duration();
		//cout << "\t translate H2H end, using time: " << time.get_microsecond_duration() << endl;


		cout << "H2H construct in GPU start " << endl;
		time.updateStart();
		makeH2HLabel_noCommunication(tempHeight, cudaThreadNum, H2H_startIndex_D, H2H_label_D, TreeBFS_D, TreeHash);
		time.updateEnd();
		H2HUsingTime_GPU = time.get_microsecond_duration();
		cout << "\t translate H2H end, using time: " << time.get_microsecond_duration() << endl;
	}



	void Graph_D_H::Graph::translateBeforeQuery()
	{

		H2H_startIndex_D = H2H_startIndex;
		H2H_label_D = H2H_label;
	}

	void Graph_D_H::Graph::translateH2H_noHub()
	{
		H2H_pos_hash_D = H2H_pos_hash;
		H2H_dis_hash_D = H2H_dis_hash;
		//H2H_pos_ID_D = H2H_pos_ID;
		H2H_pos_POS_D = H2H_pos_POS;
		H2H_dis_D = H2H_dis;
	}

	void Graph_D_H::Graph::translateH2HBack_noHub()
	{
		H2H_dis = H2H_dis_D;
	}



	long long int Graph::H2Hquery_bunch_UsingLCA_D(int size)
	{
		result_D = result;
		return H2HQuery_UsingLCA(x_D, y_D, result_D, firstAppeare_D, RMQ_OneLine_D, H2H_startIndex_D,
			H2H_label_D, RMQ_Size, RMQ_Line_Size, cudaThreadNum, size);
	}

	long long int Graph::H2Hquery_bunch_noLCA_D(int size)
	{
		result_D = result;
		return H2HQuery_NoLCA(x_D,y_D,result_D,H2H_startIndex_D,H2H_label_D,cudaThreadNum,size);
	}



	long long int Graph_D_H::Graph::H2Hquery_bunch_UsingLCA_D_noHub(int size)
	{
		//cout << "\t using before result size: " << result.size() << " result_D size: " << result_D.size() << endl;
		//result_D = result;
		//cout << "\t using after result size: " << result.size() << " result_D size: " << result_D.size() << endl;
		return H2HQuery_UsingLCA_noHub(x_D, y_D, result_D, H2H_pos_hash_D, H2H_dis_hash_D, 
				H2H_pos_POS_D, H2H_dis_D, firstAppeare_D, RMQHash_D, RMQ_ID_D, H2H_TreeHeight_Hash_D, cudaThreadNum, size);

	}


	long long int Graph_D_H::Graph::H2Hquery_bunch_noLCA_D_noHub(int size)
	{
		//cout << "\t no before result size: " << result.size() << " result_D size: " << result_D.size() << endl;
		//result_D = result; 
		//cout << "no after result size: " << result.size() << " result_D size: " << result_D.size() << endl;
		return H2HQuery_noLCA_noHub(x_D, y_D, result_D, H2H_dis_hash_D,
			H2H_dis_D, firstAppeare_D, RMQHash_D, RMQ_ID_D, H2H_TreeHeight_Hash_D, cudaThreadNum, size);

	}




	void Graph_D_H::Graph::makeH2HLabel_noHub_noComm()
	{
		Graph_D_H::time_Mine time;


		std::cout << "construct H2H on CPU multi-thread start " << std::endl;
		time.updateStart();

		for (int64_t tempHeight = 0; tempHeight < changeToGPUHeight; tempHeight++)
		{
			int64_t startIndex = TreeHash[tempHeight];
			int64_t endIndex = TreeHash[tempHeight + 1];

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
								int64_t fatherID = nodeID;
								for (int j = (int)(H2H_dis_hash[nodeID + 1] - H2H_dis_hash[nodeID]) - 2; j > pos; j--) {
									fatherID = (int64_t)father[fatherID];

									H2H_dis[indexNode + (int64_t)j] = std::min(H2H_dis[indexNode + (int64_t)j],
										tempLength_ + H2H_dis[H2H_dis_hash[fatherID] + (int64_t)pos]);
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
		time.updateEnd();
		H2HUsingTime_CPU = time.get_microsecond_duration();
		std::cout << "\t CPU construct end, using time: " << time.get_microsecond_duration() << std::endl;

		cout << "translate H2H start " << endl;
		time.updateStart();
		H2H_pos_hash_D = H2H_pos_hash;
		H2H_dis_hash_D = H2H_dis_hash;
		H2H_pos_ID_D = H2H_pos_ID;
		H2H_pos_POS_D = H2H_pos_POS;
		H2H_dis_D = H2H_dis;
		time.updateEnd();
		H2HTranslateTime = time.get_microsecond_duration();
		cout << "\t translate H2H end, using time: " << time.get_microsecond_duration() << endl;

		cout << "H2H construct in GPU start " << endl;
		time.updateStart();


		makeH2HLabel_noCommunication_noHub_2(cudaThreadNum, H2H_pos_hash_D, H2H_dis_hash_D, H2H_pos_ID_D, H2H_pos_POS_D,
			H2H_dis_D, father_D, frontier_Hash_D, frontier_ID_D, frontier_Hash.size() - 1);
		time.updateEnd();
		H2HUsingTime_GPU = time.get_microsecond_duration();
		cout << "\t construct H2H end, using time: " << time.get_microsecond_duration() << endl;
	}


	void Graph_D_H::Graph::makeH2HLabel_noHub_noComm_2()
	{
		Graph_D_H::time_Mine time;

		int tempHeight11111 = 2;
		std::cout << "construct H2H on CPU multi-thread start " << std::endl;
		time.updateStart();

		for (int64_t tempHeight = 0; tempHeight < tempHeight11111; tempHeight++)
		{
			int64_t startIndex = TreeHash[tempHeight];
			int64_t endIndex = TreeHash[tempHeight + 1];

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

		}
		time.updateEnd();
		H2HUsingTime_CPU = time.get_microsecond_duration();
		std::cout << "\t CPU construct end, using time: " << time.get_microsecond_duration() << std::endl;

		cout << "translate H2H start " << endl;
		time.updateStart();
		H2H_pos_hash_D = H2H_pos_hash;
		H2H_dis_hash_D = H2H_dis_hash;
		H2H_pos_ID_D = H2H_pos_ID;
		H2H_pos_POS_D = H2H_pos_POS;
		H2H_dis_D = H2H_dis;
		time.updateEnd();
		H2HTranslateTime = time.get_microsecond_duration();
		cout << "\t translate H2H end, using time: " << time.get_microsecond_duration() << endl;

		cout << "H2H construct in GPU start " << endl;
		H2HUsingTime_GPU = makeH2HLabel_noCommunication_noHub(tempHeight11111, cudaThreadNum, H2H_pos_hash_D, H2H_dis_hash_D, H2H_pos_ID_D, H2H_pos_POS_D,
			H2H_dis_D, father_D, TreeBFS_D, TreeHash);
		//H2HUsingTime_GPU = time.get_microsecond_duration();
		cout << "\t construct H2H end, using time: " << H2HUsingTime_GPU << endl;
	}

	void Graph_D_H::Graph::makeH2HLabel_noHub_noComm_3()
	{
		Graph_D_H::time_Mine time;

		int tempHeight11111 = 2;
		std::cout << "construct H2H on CPU multi-thread start " << std::endl;
		time.updateStart();

		for (int64_t tempHeight = 0; tempHeight < tempHeight11111; tempHeight++)
		{
			int64_t startIndex = TreeHash[tempHeight];
			int64_t endIndex = TreeHash[tempHeight + 1];

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

		}
		time.updateEnd();
		H2HUsingTime_CPU = time.get_microsecond_duration();
		std::cout << "\t CPU construct end, using time: " << time.get_microsecond_duration() << std::endl;

		cout << "translate H2H start " << endl;
		time.updateStart();

		H2H_dis_hash_D = H2H_dis_hash;
		H2H_dis_D = H2H_dis;

		TreeBFS_ID_D = TreeBFS_ID;
		TreeBFS_adj_D = TreeBFS_adj;
		TreeBFS_pos_D = TreeBFS_pos;
		//TreeBFS_changeTime_D = TreeBFS_changeTime;
		time.updateEnd();
		H2HTranslateTime = time.get_microsecond_duration();
		cout << "\t translate H2H end, using time: " << time.get_microsecond_duration() << endl;

		cout << "H2H construct in GPU start " << endl;
		H2HUsingTime_GPU = makeH2HLabel_noCommunication_noHub_3(tempHeight11111, cudaThreadNum, H2H_dis_hash_D, H2H_dis_D, 
			TreeBFS_ID_D, TreeBFS_adj_D, TreeBFS_pos_D, //TreeBFS_changeTime_D,
			father_D, TreeBFS_Hash);
		//H2HUsingTime_GPU = time.get_microsecond_duration();
		cout << "\t construct H2H end, using time: " << H2HUsingTime_GPU << endl;

		TreeBFS_ID.clear();
		TreeBFS_ID.shrink_to_fit();
		TreeBFS_ID_D.clear();
		TreeBFS_ID_D.shrink_to_fit();
		TreeBFS_adj.clear();
		TreeBFS_adj.shrink_to_fit();
		TreeBFS_adj_D.clear();
		TreeBFS_adj_D.shrink_to_fit();
		TreeBFS_pos.clear();
		TreeBFS_pos.shrink_to_fit();
		TreeBFS_pos_D.clear();
		TreeBFS_pos_D.shrink_to_fit();


		H2H_pos_ID_D = H2H_pos_ID;
		H2H_pos_POS_D = H2H_pos_POS;
		H2H_pos_hash_D = H2H_pos_hash;
	}

}