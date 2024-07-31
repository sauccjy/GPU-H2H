#pragma once
#include"Graph_D.cuh"

namespace Graph_D_H
{
	void Graph_D_H::Graph::beforeH2H()
	{
		getrank();
		getTreeStructure();
		if (useRMQ)
			makeOULAAndRMQ();
		else
			onlyOULA();

		cout << "\t CHTreeHeight: " << CHTreeHeight << endl;
		TDTREEHeight = CHTreeHeight;
		cout << "\t RMQ size : " << ((double)(RMQ_Line_Size * RMQ_Size * sizeof(Graph_D_H::pairs))) / (1024 * 1024)
			<< "MB, H2H Label size: " << ((double)H2H_Size * sizeof(Graph_D_H::H2HLabel)) / (1024 * 1024) << "MB" << endl;
		RMQsize = ((double)(RMQ_Line_Size * RMQ_Size * sizeof(Graph_D_H::pairs))) / (1024 * 1024);
		H2HLabelSize = ((double)H2H_Size * sizeof(Graph_D_H::H2HLabel)) / (1024 * 1024);
		father_D = father;
		cleanBeforeConstruct();
	}

	void Graph_D_H::Graph::cleanBeforeConstruct()
	{
		CHInfo.clear();
		CHInfo.shrink_to_fit();
		//CHChild.clear();
		//CHChild.shrink_to_fit();
		ID_hash.clear();
		ID_hash.shrink_to_fit();
		ID_hash_D.clear();
		ID_hash_D.shrink_to_fit();
		visited.clear();
		visited.shrink_to_fit();
		visited_D.clear();
		visited_D.shrink_to_fit();
	}

	void Graph_D_H::Graph::beforeH2H_noHub()
	{
		getrank();
		getTreeStructure();
		if (useRMQ)
			makeOULAAndRMQ();
		else
			onlyOULA();

		father_D = father;
		cout << "H2H_Size: " << (double)(H2H_Size * sizeof(int)) / (1024 * 1024) << endl;
		cleanBeforeConstruct();
	}

	void Graph_D_H::Graph::getrank()
	{
		ranks.assign(NodeNumber,-1);
		for (int i = 0; i < nonAdjcentNode.size(); i++)
		{
			ranks[nonAdjcentNode[i].second] = i;
		}
		head = nonAdjcentNode[nonAdjcentNode.size() - 1].second;

		nonAdjcentNode.clear();
		nonAdjcentNode_D.clear();
		NANHash.clear();
		NANHash_D.clear();
		nonAdjcentNode.shrink_to_fit();
		nonAdjcentNode_D.shrink_to_fit();
		NANHash.shrink_to_fit();
		NANHash_D.shrink_to_fit();
	}

	void Graph_D_H::Graph::getTreeStructure() {
		father.assign(NodeNumber, -1);

		ChildHash.assign(father.size(), pairs(0, -1));
		thrust::host_vector<int> ChildHash_actual_pos(father.size(), 0);
		Childs.assign(father.size() - 1, -1);
		for (int i = 0; i < NodeNumber; i++) {
			int fatherID = -1;
			int tempRank = INT_MAX;
			for (auto& it : CHAdjlist[i]) {
				int adjID = it.first;
				if (tempRank > ranks[adjID])
				{
					tempRank = ranks[adjID];
					fatherID = adjID;
				}
			}

			father[i] = fatherID;
		}

		ranks.clear();
		ranks.shrink_to_fit();
		for (int i = 0; i < NodeNumber; i++) {
			if (father[i] == -1)
				continue;
			ChildHash[father[i]].first++;
		}
		int temp = 0;
		for (int i = 0; i < NodeNumber; i++)
		{
			ChildHash[i].second = temp;
			temp += ChildHash[i].first;
		}
		for (int i = 0; i < NodeNumber; i++)
		{
			if (father[i] == -1)
				continue;
			//int fatherID = father[i];
			Childs[ChildHash[father[i]].second + ChildHash_actual_pos[father[i]]] = i;
			ChildHash_actual_pos[father[i]]++;
		}
	}
	
	void Graph::onlyOULA()
	{
		oulaOnly = true;
		OULA_Only.assign(NodeNumber, 0);
		thrust::host_vector<bool> visit(NodeNumber, false);
		if (head == -1)
			return;
		onlyOULADFS(head, 1, visit);
	}

	void Graph::onlyOULADFS(const int headNow, const int depth, thrust::host_vector<bool>& visit)
	{
		OULA_Only[headNow] = depth;
		H2H_Size += depth;
		visit[headNow] = true;
		CHTreeHeight = max(CHTreeHeight, depth);
		for (int i = ChildHash[headNow].second; i < ChildHash[headNow].second + ChildHash[headNow].first; i++)
		{
			if (!visit[Childs[i]])
			{
				onlyOULADFS(Childs[i], depth + 1, visit);
			}
		}
	}

	void Graph::makeOULAAndRMQ()
	{
		thrust::host_vector<bool> visit(NodeNumber, false);
		firstAppeare.assign(NodeNumber, -1);
		H2H_TreeHeight_Hash.assign(NodeNumber, -1);
		if (head == -1)
			return;

		OULADFS(head, 1, visit);

		//construct RMQ by DP
		int jSize = (int)(log2((double)OULA_DFS_.size()));

		RMQ_Size = OULA_DFS_.size();
		RMQ_Line_Size = jSize + 1;

		//RMQ_OneLine = new pairs[RMQ_Line_Size * RMQ_Size];
		RMQ_OneLine.assign(RMQ_Line_Size * RMQ_Size, pairs());
		for (auto i = 0; i < RMQ_Size; i++)
		{
			RMQ_OneLine[i * RMQ_Line_Size].first = OULA_DFS_[i].second;
			RMQ_OneLine[i * RMQ_Line_Size].second = OULA_DFS_[i].first;
		}
		for (int j = 1; j <= jSize; j++)
		{
			auto davicate = (1 << (j - 1));
			for (auto i = 0; i + davicate < RMQ_Size; i++)
			{
				if (RMQ_OneLine[i * RMQ_Line_Size + j - 1].first <= RMQ_OneLine[(i + davicate) * RMQ_Line_Size + j - 1].first)
					RMQ_OneLine[i * RMQ_Line_Size + j].pairsCopy(RMQ_OneLine[i * RMQ_Line_Size + j - 1]);
				else
					RMQ_OneLine[i * RMQ_Line_Size + j].pairsCopy(RMQ_OneLine[(i + davicate) * RMQ_Line_Size + j - 1]);
			}
		}

		RMQHash.assign(OULA_DFS_.size(), -1);
		for (int64_t i = 0; i < RMQ_Size; i++) {
			RMQHash[i] = RMQ_ID.size();
			for (int j = 0; j < RMQ_Line_Size; j++) {
				if (RMQ_OneLine[i * RMQ_Line_Size + j].first == -1) {
					break;
				}
				//RMQ_Height.push_back(RMQ_OneLine[i * RMQ_Line_Size + j].first);
				RMQ_ID.push_back(RMQ_OneLine[i * RMQ_Line_Size + j].second);
			}
		}
		RMQ_OneLine.clear();
		RMQ_OneLine.shrink_to_fit();


	}

	void Graph::OULADFS(const int headNow, const int depth, thrust::host_vector<bool>& visit)
	{
		OULA_DFS_.push_back(pairs(headNow, depth));
		firstAppeare[headNow] = OULA_DFS_.size() - 1;
		H2H_TreeHeight_Hash[headNow] = depth;
		H2H_Size += depth;
		visit[headNow] = true;
		CHTreeHeight = max(CHTreeHeight, depth);
		//int child_hash_index = ID_hash[headNow];
		int CHildStart = ChildHash[headNow].second;
		int CHildEnd = CHildStart + ChildHash[headNow].first;
		for (int i = CHildStart; i < CHildEnd; i++)
		{
			int childID = Childs[i];
			if (!visit[childID])
			{
				OULADFS(childID, depth + 1, visit);
				OULA_DFS_.push_back(pairs(headNow, depth));
			}
		}
	}
}