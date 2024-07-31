#pragma once
#include"Graph_D.cuh"

namespace Graph_D_H
{
	void Graph_D_H::Graph::displayAdjList()
	{
		for (int i = 0; i < adjList.size(); i++){
			cout << "Node ID : " << i << " ";
			for (int j = 0; j < adjList[i].size(); j++)
				cout << "<" << adjList[i][j].first << "," << adjList[i][j].second << "> ";
			cout << endl;
		}
	}

	void Graph_D_H::Graph::displayGraph(std::map<int, std::map<int, int>>& standard_graph_in_partition)
	{
		cout << "\tsubGraph" << endl;
		for (auto& i : standard_graph_in_partition) {
			cout << "\tvertex id: " << i.first<<" and adjscent: ";
			for (auto& j : standard_graph_in_partition.at(i.first)) {
				cout << j.first << ",";
			}
			cout << endl;
		}

	}

	void Graph_D_H::Graph::displayPartition()
	{
		cout << "______________________PARTITION________________________" << endl;
		int tempHeight = 1;
		while (tempHeight <= mainTreeHeight)
		{
			cout << "height:" << tempHeight << ";  ";
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 2;
			tempHeight++;
			for (int i = lowestIndexStart; i <= lowestIndexEnd; i++)
			{
				int left = partition_Tree[i].first;
				int right = partition_Tree[i].second;
				cout << "(";
				for (int j = left; j < right - 1; j++)
				{
					cout << NE_P[j].NodeID << ",";
				}
				cout << NE_P[right - 1].NodeID << ")";
			}
			cout << endl;
		}
		cout << "______________________ID_HASH________________________" << endl;
		for (int i = 0; i < NodeNumber; i++)
		{
			cout << "ID:" << i << " hash:" << ID_hash[i] << "; ";
		}
		cout << endl;
	}

	void Graph_D_H::Graph::displayOrderAndCHTree()
	{
		cout << "______________________ORDER_CHTREE________________________" << endl;
		int tempHeight = TreeHeight;
		while (tempHeight > 0)
		{
			cout << "height:" << tempHeight << ";  "<<endl;
			int lowestIndexStart = std::pow(2, tempHeight - 1) - 1;
			int lowestIndexEnd = std::pow(2, tempHeight) - 2;
			tempHeight--;
			for (int i = lowestIndexStart; i <= lowestIndexEnd; i++)
			{
				int nonAdjcentNode_index = NANHash[i].first;
				int nonAdjcentNode_index_end = nonAdjcentNode_index + NANHash[i].NodeID;
				cout << "start index: " << nonAdjcentNode_index << "; partition Node Size: " << nonAdjcentNode_index_end <<
					"; partition node max size: " << NANHash[i].second<< "; nonadjcent Node: ";
				for (int j = nonAdjcentNode_index; j < nonAdjcentNode_index_end; j++)
				{
					cout <<"<"<< nonAdjcentNode[j].first << "," << nonAdjcentNode[j].second << "> ";
				}
				cout << endl;
			}
			cout << endl;
		}
		cout << "______________________Rank________________________" << endl;
		for (int i = 0; i < nonAdjcentNode.size(); i++)
		{
			cout << nonAdjcentNode[i].second << " ";
		}
		cout << endl;

		cout << "______________________ID_HASH________________________" << endl;
		for (int i = 0; i < NodeNumber; i++)
		{
			cout<<"ID:"<<i<<" hash:" << ID_hash[i] << "; ";
		}
		cout << endl;

		cout << "______________________Tree_OLD________________________" << endl;
		for (int i = 0; i < CHTreeHash.size(); i++)
		{
			cout << "NodeID: " << NE_P[ID_hash[i]].NodeID << ": actualSize: "<<CHTreeHash[ID_hash[i]].first<<
				" maxSize: "<< CHTreeHash[ID_hash[i]].second<<endl;
			for (int j = CHTreeHash[ID_hash[i]].NodeID; j < CHTreeHash[ID_hash[i]].NodeID + CHTreeHash[ID_hash[i]].second; j++)
			{
				if (CHTree[j].first == INT_MAX)
					continue;
				cout << "(" << CHTree[j].first << "," << CHTree[j].second << "," << CHTree[j].NodeID << ") ";
			}
			cout << endl;
		}


	}

	void Graph_D_H::Graph::displayCHTree()
	{
		cout << "______________________Tree_AFTER_CH________________________" << endl;
		for (int i = 0; i < CHTreeHash.size(); i++)
		{
			cout << "NodeID: " << NE_P[ID_hash[i]].NodeID << ": actualSize: " << CHTreeHash[ID_hash[i]].first <<
				" maxSize: " << CHTreeHash[ID_hash[i]].second << endl;

			for (int j = CHTreeHash[ID_hash[i]].NodeID; j < CHTreeHash[ID_hash[i]].NodeID + CHTreeHash[ID_hash[i]].second; j++)
			{
				if (CHTree[j].first == INT_MAX)
					continue;
				cout << "(" << CHTree[j].first << "," << CHTree[j].second << "," << CHTree[j].NodeID << ") ";
			}
			cout << endl;
		}
	}

	void Graph_D_H::Graph::displayOULARMQ()
	{
		cout << "____________OULA Sequence_____________" << endl;
		for (auto i = 0; i < OULA_DFS_.size(); i++)
		{
			cout << " <" << OULA_DFS_[i].first << "," << OULA_DFS_[i].second << ">";
		}
		cout << endl;
		cout << "_____________FirstAppare______________" << endl;
		for (auto i = 0; i < firstAppeare.size(); i++)
		{
			cout << " <" << i << "," << firstAppeare[i] << ">";

		}

		cout << endl;
		cout << "______________RMQ_____________________" << endl;
		for (int j = 0; j < RMQ_Line_Size; j++)
		{
			cout << "j = " << j << " jump (1<<j)-1" << endl;
			for (int i = 0; i < RMQ_Size; i++)
			{
				if (RMQ_OneLine[i * RMQ_Line_Size + j].first != INT_MAX)
					cout << " <" << i << "," << RMQ_OneLine[i * RMQ_Line_Size + j].first << "," << RMQ_OneLine[i * RMQ_Line_Size + j].second << "> ";
			}
			cout << endl;
		}


	}

	void Graph_D_H::Graph::displayH2HLabel()
	{
		cout << "________________________H2H_CSR_LABEL__________________________" << endl;
		for (int i = 0; i < NodeNumber; i++)
		{
			cout << "NodeID: " << i << endl;
			int pos = H2H_startIndex[i].first;
			int size = H2H_startIndex[i].second;
			cout << "Node: ";
			for (int j = 0; j < size; j++)
			{
				cout << H2H_label[pos + j].Node << " ";
			}
			cout << endl;
			cout << "isTreeNode: ";
			for (int j = 0; j < size; j++)
			{
				cout << H2H_label[pos + j].isTreeNode << " ";
			}
			cout << endl;
			cout << "Hub: ";
			for (int j = 0; j < size; j++)
			{
				cout << H2H_label[pos + j].Hub << " ";
			}
			cout << endl;
			cout << "dis: ";
			for (int j = 0; j < size; j++)
			{
				cout << H2H_label[pos + j].dis << " ";
			}
			cout << endl;
			cout << "-------------------------------" << endl;
		}
	}

	void Graph_D_H::Graph::displaySingalLabel(int i)
	{
		cout << "NodeID: " << i << endl;
		int pos = H2H_startIndex[i].first;
		int size = H2H_startIndex[i].second;
		cout << "Node: ";
		for (int j = 0; j < size; j++)
		{
			cout << H2H_label[pos + j].Node << " ";
		}
		cout << endl;
		cout << "isTreeNode: ";
		for (int j = 0; j < size; j++)
		{
			cout << H2H_label[pos + j].isTreeNode << " ";
		}
		cout << endl;
		cout << "Hub: ";
		for (int j = 0; j < size; j++)
		{
			cout << H2H_label[pos + j].Hub << " ";
		}
		cout << endl;
		cout << "dis: ";
		for (int j = 0; j < size; j++)
		{
			cout << H2H_label[pos + j].dis << " ";
		}
		cout << endl;
		cout << "-------------------------------" << endl;
	}


	void Graph::displayCHT()
	{
		for(int j = 1;j < 10 ;j++)
		{
			int i = nonAdjcentNode[nonAdjcentNode.size() - j].second;
			cout << "NodeID: " << NE_P[ID_hash[i]].NodeID << ": actualSize: " << CHTreeHash[ID_hash[i]].first <<
				" maxSize: " << CHTreeHash[ID_hash[i]].second << endl;

			for (int j = CHTreeHash[ID_hash[i]].NodeID; j < CHTreeHash[ID_hash[i]].NodeID + CHTreeHash[ID_hash[i]].first; j++)
			{
				if (CHTree[j].first == INT_MAX)
					continue;
				cout << "(" << CHTree[j].first << "," << CHTree[j].second << "," << CHTree[j].NodeID << ") ";
			}
			cout << endl;
		}
	}

	void Graph::displayCHTree_mapVersion() {
		for (int j = 1; j <= 13; j++)
		{
			int i = nonAdjcentNode[nonAdjcentNode.size() - j].second;
			cout << "NodeID: " << NE_P[ID_hash[i]].NodeID << endl;

			for (auto& it : CHAdjlist[i]) {
				cout << "(" << it.first << "," << it.second.first << "," << it.second.second << ") ";
			}
			cout << endl;
		}
	}

	void Graph_D_H::Graph::displayCHAdjList()
	{
		cout << "_________________CHTree_____________________" << endl;
		for (int i = 0; i < CHAdjlist.size(); i++) {
			cout << "\tNodeID: " << i << " : ";
			for (auto& it : CHAdjlist[i]) {
				cout << "<" << it.first << "," << it.second.first << "," << it.second.second << ">";
			}
			cout << endl;
		}
	}

	void Graph_D_H::Graph::displayH2H_noHub()
	{
		cout << "________________________H2H label no Hub_____________" << endl;

		cout << "pos hash (id, index): ";
		for (int i = 0; i < H2H_pos_hash.size(); i++) {
			cout << "(" << i << "," << H2H_pos_hash[i] << ")";
		}
		cout << endl;

		cout << "pos ID and POS (index, id, pos): ";
		for (int i = 0; i < H2H_pos_ID.size(); i++) {
			cout << "(" << i << "," << H2H_pos_ID[i] << "," << H2H_pos_POS[i]<<")";
		}
		cout << endl;

		cout << "dis hash (id, index): ";
		for (int i = 0; i < H2H_dis_hash.size(); i++) {
			cout << "(" << i << "," << H2H_dis_hash[i] << ")";
		}
		cout << endl;

		cout << "dis (index, dis): ";
		for (int i = 0; i < H2H_dis.size(); i++) {
			cout << "(" << i << "," << H2H_dis[i] << ")";
		}
		cout << endl;


		cout << endl;
		cout << "____________H2H label no hub nodeID version______" << endl;
		for (int i = 0; i < NodeNumber; i++) {
			cout << "nodeID: " << i << endl;
			cout << "pos: "<<endl;
			cout << "\tID: ";
			for (int64_t j = H2H_pos_hash[i]; j < (int64_t)H2H_pos_hash[i + 1]; j++) {
				cout << H2H_pos_ID[j] << " ";
			}
			cout << endl;
			cout << "\tPOS: ";
			for (int64_t j = H2H_pos_hash[i]; j < (int64_t)H2H_pos_hash[i + 1]; j++) {
				cout << H2H_pos_POS[j] << " ";
			}
			cout << endl;
			cout << "dis: " << endl;
			cout << "\t";
			for (int64_t j = H2H_dis_hash[i]; j < (int64_t)H2H_dis_hash[i + 1]; j++) {
				cout << H2H_dis[j] << " ";
			}
			cout << endl;
			cout << endl;
		}
	}
}