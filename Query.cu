#pragma once
#include"Graph_D.cuh"

namespace Graph_D_H
{
	void Graph_D_H::Graph::generateQuery()
	{
		x.assign(queryMaxSize,0);
		y.assign(queryMaxSize,0);
		result.assign(queryMaxSize,INT_MAX);

		std::default_random_engine generator(std::time(0));
		std::uniform_int_distribution<int> distribution(0, NodeNumber - 1);
		for (int i = 0; i < queryMaxSize; ++i) {
			x[i] = distribution(generator);
			y[i] = distribution(generator);
		}
	}


	void Graph_D_H::Graph::uploadRealQuery()
	{
		fstream queryFile1("./Query/NY/normal.txt", ios::in | ios::out);
		x.clear();
		y.clear();
		for (int i = 0; i < 100000; i++) {
			int start, des, other;
			queryFile1 >> start >> des >> other;
			x.push_back(start);
			y.push_back(des);
		}
		queryFile1.close();

		fstream queryFile2("./Query/NY/sparse.txt", ios::in | ios::out);
		for (int i = 0; i < 100000; i++) {
			int start, des, other;
			queryFile2 >> start >> des >> other;
			x.push_back(start);
			y.push_back(des);
		}
		queryFile2.close();
		queryMaxSize = x.size();
		result.assign(x.size(), INT_MAX);
	}


	long long int Graph::H2Hquery_bunch(int size)
	{
		Graph_D_H::time_Mine time;
		cout << "H2H_sequence start query,querySize=" <<size<< endl;
		time.updateStart();

		for (int i = 0; i < size; i++) {
			result[i] = H2HDistancQuery_UsingLCA(x[i], y[i]);
		}

		time.updateEnd();
		cout << "using time: " << time.get_microsecond_duration() << "us" << endl;

		result.assign(queryMaxSize, INT_MAX);

		return time.get_microsecond_duration();

	}

	long long int Graph::H2Hquery_bunch_MultiThread(int size) {
		const unsigned int numThreads = threadNumber;
		const size_t chunkSize = size / numThreads + 1;
		std::vector<std::thread> threads;

		Graph_D_H::time_Mine time;

		cout << "H2H_multithread start query, thread size : "<< numThreads <<" ,querySize=" << size << endl;
		time.updateStart();

		for (int i = 0; i < numThreads; ++i) {
			size_t start = i * chunkSize;
			size_t end = (i == numThreads - 1) ? size : (i + 1) * chunkSize;
			end = min((size_t)size, (i + 1) * chunkSize);
			threads.emplace_back([&, start, end] {
				for (size_t j = start; j < end; ++j) {
					result[j] = H2HDistancQuery_UsingLCA(x[j], y[j]);
				}
				});
		}
		for (auto& thread : threads) {
			thread.join();
		}

		time.updateEnd();
		cout << "using time: " << time.get_microsecond_duration() << "us" << endl;

		result.assign(queryMaxSize, INT_MAX);
		return time.get_microsecond_duration();
	}


	void Graph_D_H::Graph::cleanResult()
	{

	}





	//noHub 
	int Graph_D_H::Graph::H2HDistancQuery_UsingLCA_noHub(int x, int y)
	{
		if (x == y)
			return 0;
		int LCA = findLCA_noHub(x, y);
		int resu = INT_MAX;
		for (int i = H2H_pos_hash[LCA]; i < H2H_pos_hash[LCA + 1]; i++) {

			int templength = H2H_dis[H2H_dis_hash[x] + H2H_pos_POS[i]] + H2H_dis[H2H_dis_hash[y] + H2H_pos_POS[i]];
			if (templength < resu) {
				resu = templength;
			}
		}
		return resu;
	}

	int Graph_D_H::Graph::H2HDistancQuery_noLCA_noHub(int x, int y)
	{
		if (x == y)
			return 0;
		int FirstA = firstAppeare[x];//position  in  OULA_DFS
		int FirstB = firstAppeare[y];
		int firstA = (FirstA < FirstB) * FirstA + (FirstB <= FirstA) * FirstB;
		int firstB = (FirstA < FirstB) * FirstB + (FirstB <= FirstA) * FirstA;
		int k = (int)(log2((double)(firstB - firstA + 1)));//k means log2(gap);

		int ID_1 = RMQ_ID[RMQHash[firstA] + (int64_t)k];
		int ID_2 = RMQ_ID[RMQHash[firstB - (1 << k) + 1] + (int64_t)k];
		int LCASize = H2H_TreeHeight_Hash[ID_2];
		if (H2H_TreeHeight_Hash[ID_1] <= H2H_TreeHeight_Hash[ID_2]) {
			LCASize = H2H_TreeHeight_Hash[ID_1];
		}
		int resu = INT_MAX;


		for (int i = 0; i < LCASize; i++) {

			int templength = H2H_dis[H2H_dis_hash[x] + (int64_t)i] + H2H_dis[H2H_dis_hash[y] + (int64_t)i];
			if (templength < resu) {
				resu = templength;
			}
		}
		return resu;
	}

	long long int Graph_D_H::Graph::H2Hquery_bunch_noHub(int size)
	{
		Graph_D_H::time_Mine time;
		cout << "H2H_sequence start query,querySize=" << size << endl;
		time.updateStart();

		for (int i = 0; i < size; i++) {
			result[i] = H2HDistancQuery_UsingLCA_noHub(x[i], y[i]);
		}

		time.updateEnd();
		cout << "using time: " << time.get_microsecond_duration() << "us" << endl;

		result.assign(queryMaxSize, INT_MAX);

		return time.get_microsecond_duration();
	}


	long long int Graph_D_H::Graph::H2Hquery_bunch_MultiThread_noHub(int size)
	{
		const unsigned int numThreads = threadNumber;
		const size_t chunkSize = size / numThreads + 1;
		std::vector<std::thread> threads;

		Graph_D_H::time_Mine time;

		cout << "H2H_multithread start query, thread size : " << numThreads << " ,querySize=" << size << endl;
		time.updateStart();

		for (int i = 0; i < numThreads; ++i) {
			size_t start = i * chunkSize;
			size_t end = (i == numThreads - 1) ? size : (i + 1) * chunkSize;
			end = min((size_t)size, (i + 1) * chunkSize);
			threads.emplace_back([&, start, end] {
				for (size_t j = start; j < end; ++j) {
					result[j] = H2HDistancQuery_UsingLCA_noHub(x[j], y[j]);
				}
				});
		}
		for (auto& thread : threads) {
			thread.join();
		}

		time.updateEnd();
		cout << "using time: " << time.get_microsecond_duration() << "us" << endl;

		result.assign(queryMaxSize, INT_MAX);
		return time.get_microsecond_duration();
	}

	void Graph_D_H::Graph::answerRandomQuery()
	{
		generateQuery();
		translateQuery();
		int endTime = 500000;
		bool seqContinue = true;
		bool MTContinue_4 = true;
		bool MTContinue_6 = true;
		bool MTContinue_8 = true;
		bool MTContinue_10 = true;
		bool MTContinue_12 = true;
		bool MTContinue_14 = true;
		bool MTContinue_16 = true;
		for (int i = 10000; i <= queryMaxSize; i += queryStep) {

			if (seqContinue) {
				long long int seqtime = H2Hquery_bunch(i);
				QueryTime_seq.push_back(seqtime);
				seqContinue = (seqtime < endTime);
			}
			else {
				QueryTime_seq.push_back(endTime);
			}

			if (MTContinue_4) {
				threadNumber = 4;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_4.push_back(seqtime);
				MTContinue_4 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_4.push_back(endTime);
			}

			if (MTContinue_6) {
				threadNumber = 6;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_6.push_back(seqtime);
				MTContinue_6 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_6.push_back(endTime);
			}

			if (MTContinue_8) {
				threadNumber = 8;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_8.push_back(seqtime);
				MTContinue_8 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_8.push_back(endTime);
			}

			if (MTContinue_10) {
				threadNumber = 10;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_10.push_back(seqtime);
				MTContinue_10 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_10.push_back(endTime);
			}

			if (MTContinue_12) {
				threadNumber = 12;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_12.push_back(seqtime);
				MTContinue_12 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_12.push_back(endTime);
			}

			if (MTContinue_14) {
				threadNumber = 14;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_14.push_back(seqtime);
				MTContinue_14 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_14.push_back(endTime);
			}

			if (MTContinue_16) {
				threadNumber = 16;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_16.push_back(seqtime);
				MTContinue_16 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_16.push_back(endTime);
			}
			//QueryTime_seq.push_back(H2Hquery_bunch(i));
			//QueryTime_MT.push_back(H2Hquery_bunch_MultiThread(i));
			cudaThreadNum = 64;
			QueryTime_noLCA_64.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_64.push_back(H2Hquery_bunch_UsingLCA_D(i));

			cudaThreadNum = 128;
			QueryTime_noLCA_128.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_128.push_back(H2Hquery_bunch_UsingLCA_D(i));

			cudaThreadNum = 256;
			QueryTime_noLCA_256.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_256.push_back(H2Hquery_bunch_UsingLCA_D(i));

			cudaThreadNum = 512;
			QueryTime_noLCA_512.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_512.push_back(H2Hquery_bunch_UsingLCA_D(i));

			cudaThreadNum = 1024;
			QueryTime_noLCA_1024.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_1024.push_back(H2Hquery_bunch_UsingLCA_D(i));
		}
		string resultFile = "./QueryResult/" + graphName + "/PTHeight-" + to_string(TreeHeight) + "-ChangeHeight-" + to_string(changeHeight) + "-random.csv";
		fstream queryresult(resultFile, ios::app | ios::in | ios::out);
		if (queryresult.is_open()) {
			//queryresult << "CPUThreadNum,CUDAThreadNum\n";
			//queryresult <<threadNumber<<","<< cudaThreadNum << "\n";
			//if(queryresult.peek() == std::ifstream::traits_type::eof())
			queryresult << "queryNumber,sequence,multiThread-4,multiThread-6,multiThread-8,multiThread-10,multiThread-12,multiThread-14,multiThread-16,GPUnoLCA-64,GPUwithLCA-64,GPUnoLCA-128,GPUwithLCA-128,GPUnoLCA-256,GPUwithLCA-256,GPUnoLCA-512,GPUwithLCA-512,GPUnoLCA-1024,GPUwithLCA-1024\n";
			for (int i = 0; i < QueryTime_seq.size(); i++) {
				queryresult << 10000 + i * queryStep << "," << QueryTime_seq[i] << "," << QueryTime_MT_4[i] << ","
					<< QueryTime_MT_6[i] << "," << QueryTime_MT_8[i] << "," << QueryTime_MT_10[i] << "," << QueryTime_MT_12[i] << ","
					<< QueryTime_MT_14[i] << "," << QueryTime_MT_16[i] << ","
					<< QueryTime_noLCA_64[i] << "," << QueryTime_LCA_64[i] << ","
					<< QueryTime_noLCA_128[i] << "," << QueryTime_LCA_128[i] << ","
					<< QueryTime_noLCA_256[i] << "," << QueryTime_LCA_256[i] << ","
					<< QueryTime_noLCA_512[i] << "," << QueryTime_LCA_512[i] << ","
					<< QueryTime_noLCA_1024[i] << "," << QueryTime_LCA_1024[i] << "\n";
			}
			queryresult.close();
		}

	}

	void Graph_D_H::Graph::answerRealNYQuery()
	{
		uploadRealQuery();
		translateQuery();
		int endTime = 500000;
		bool seqContinue = true;
		bool MTContinue_4 = true;
		bool MTContinue_6 = true;
		bool MTContinue_8 = true;
		bool MTContinue_10 = true;
		bool MTContinue_12 = true;
		bool MTContinue_14 = true;
		bool MTContinue_16 = true;
		for (int i = 10000; i <= queryMaxSize; i += queryStep) {

			if (seqContinue) {
				long long int seqtime = H2Hquery_bunch(i);
				QueryTime_seq.push_back(seqtime);
				seqContinue = (seqtime < endTime);
			}
			else {
				QueryTime_seq.push_back(endTime);
			}

			if (MTContinue_4) {
				threadNumber = 4;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_4.push_back(seqtime);
				MTContinue_4 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_4.push_back(endTime);
			}

			if (MTContinue_6) {
				threadNumber = 6;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_6.push_back(seqtime);
				MTContinue_6 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_6.push_back(endTime);
			}

			if (MTContinue_8) {
				threadNumber = 8;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_8.push_back(seqtime);
				MTContinue_8 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_8.push_back(endTime);
			}

			if (MTContinue_10) {
				threadNumber = 10;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_10.push_back(seqtime);
				MTContinue_10 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_10.push_back(endTime);
			}

			if (MTContinue_12) {
				threadNumber = 12;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_12.push_back(seqtime);
				MTContinue_12 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_12.push_back(endTime);
			}

			if (MTContinue_14) {
				threadNumber = 14;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_14.push_back(seqtime);
				MTContinue_14 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_14.push_back(endTime);
			}

			if (MTContinue_16) {
				threadNumber = 16;
				long long int seqtime = H2Hquery_bunch_MultiThread(i);
				QueryTime_MT_16.push_back(seqtime);
				MTContinue_16 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_16.push_back(endTime);
			}
			//QueryTime_seq.push_back(H2Hquery_bunch(i));
			//QueryTime_MT.push_back(H2Hquery_bunch_MultiThread(i));
			cudaThreadNum = 64;
			QueryTime_noLCA_64.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_64.push_back(H2Hquery_bunch_UsingLCA_D(i));

			cudaThreadNum = 128;
			QueryTime_noLCA_128.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_128.push_back(H2Hquery_bunch_UsingLCA_D(i));

			cudaThreadNum = 256;
			QueryTime_noLCA_256.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_256.push_back(H2Hquery_bunch_UsingLCA_D(i));

			cudaThreadNum = 512;
			QueryTime_noLCA_512.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_512.push_back(H2Hquery_bunch_UsingLCA_D(i));

			cudaThreadNum = 1024;
			QueryTime_noLCA_1024.push_back(H2Hquery_bunch_noLCA_D(i));
			QueryTime_LCA_1024.push_back(H2Hquery_bunch_UsingLCA_D(i));
		}
		string resultFile = "./QueryResult/" + graphName + "/PTHeight-" + to_string(TreeHeight) + "-ChangeHeight-" + to_string(changeHeight) + "-real.csv";
		fstream queryresult(resultFile, ios::app | ios::in | ios::out);
		if (queryresult.is_open()) {
			//queryresult << "CPUThreadNum,CUDAThreadNum\n";
			//queryresult <<threadNumber<<","<< cudaThreadNum << "\n";
			//if(queryresult.peek() == std::ifstream::traits_type::eof())
			queryresult << "queryNumber,sequence,multiThread-4,multiThread-6,multiThread-8,multiThread-10,multiThread-12,multiThread-14,multiThread-16,GPUnoLCA-64,GPUwithLCA-64,GPUnoLCA-128,GPUwithLCA-128,GPUnoLCA-256,GPUwithLCA-256,GPUnoLCA-512,GPUwithLCA-512,GPUnoLCA-1024,GPUwithLCA-1024\n";
			for (int i = 0; i < QueryTime_seq.size(); i++) {
				queryresult << 10000 + i * queryStep << "," << QueryTime_seq[i] << "," << QueryTime_MT_4[i] << ","
					<< QueryTime_MT_6[i] << "," << QueryTime_MT_8[i] << "," << QueryTime_MT_10[i] << "," << QueryTime_MT_12[i] << ","
					<< QueryTime_MT_14[i] << "," << QueryTime_MT_16[i] << ","
					<< QueryTime_noLCA_64[i] << "," << QueryTime_LCA_64[i] << ","
					<< QueryTime_noLCA_128[i] << "," << QueryTime_LCA_128[i] << ","
					<< QueryTime_noLCA_256[i] << "," << QueryTime_LCA_256[i] << ","
					<< QueryTime_noLCA_512[i] << "," << QueryTime_LCA_512[i] << ","
					<< QueryTime_noLCA_1024[i] << "," << QueryTime_LCA_1024[i] << "\n";
			}
			queryresult.close();
		}
	}



	void Graph_D_H::Graph::answerRandomQuery_noHub()
	{
		generateQuery();
		translateQuery();

		

		int endTime = 500000;
		bool seqContinue = true;
		bool MTContinue_4 = true;
		bool MTContinue_6 = true;
		bool MTContinue_8 = true;
		bool MTContinue_10 = true;
		bool MTContinue_12 = true;
		bool MTContinue_14 = true;
		bool MTContinue_16 = true;
		for (int i = 10000; i <= queryMaxSize; i += queryStep) {

			if (seqContinue) {
				long long int seqtime = H2Hquery_bunch_noHub(i);
				QueryTime_seq.push_back(seqtime);
				seqContinue = (seqtime < endTime);
			}
			else {
				QueryTime_seq.push_back(endTime);
			}

			if (MTContinue_4) {
				threadNumber = 4;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_4.push_back(seqtime);
				MTContinue_4 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_4.push_back(endTime);
			}

			if (MTContinue_6) {
				threadNumber = 6;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_6.push_back(seqtime);
				MTContinue_6 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_6.push_back(endTime);
			}

			if (MTContinue_8) {
				threadNumber = 8;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_8.push_back(seqtime);
				MTContinue_8 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_8.push_back(endTime);
			}

			if (MTContinue_10) {
				threadNumber = 10;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_10.push_back(seqtime);
				MTContinue_10 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_10.push_back(endTime);
			}

			if (MTContinue_12) {
				threadNumber = 12;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_12.push_back(seqtime);
				MTContinue_12 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_12.push_back(endTime);
			}

			if (MTContinue_14) {
				threadNumber = 14;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_14.push_back(seqtime);
				MTContinue_14 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_14.push_back(endTime);
			}

			if (MTContinue_16) {
				threadNumber = 16;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_16.push_back(seqtime);
				MTContinue_16 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_16.push_back(endTime);
			}
			//QueryTime_seq.push_back(H2Hquery_bunch(i));
			//QueryTime_MT.push_back(H2Hquery_bunch_MultiThread(i));
			cudaThreadNum = 64;
			QueryTime_noLCA_64.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_64.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));

			cudaThreadNum = 128;
			QueryTime_noLCA_128.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_128.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));

			cudaThreadNum = 256;
			QueryTime_noLCA_256.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_256.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));

			cudaThreadNum = 512;
			QueryTime_noLCA_512.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_512.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));

			cudaThreadNum = 1024;
			QueryTime_noLCA_1024.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_1024.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));
		}
		string resultFile = "./QueryResult/" + graphName + "/PTHeight-" + to_string(TreeHeight) + "-ChangeHeight-" + to_string(changeHeight) + "-random.csv";
		fstream queryresult(resultFile, ios::app | ios::in | ios::out);
		if (queryresult.is_open()) {
			//queryresult << "CPUThreadNum,CUDAThreadNum\n";
			//queryresult <<threadNumber<<","<< cudaThreadNum << "\n";
			//if(queryresult.peek() == std::ifstream::traits_type::eof())
			queryresult << "queryNumber,sequence,multiThread-4,multiThread-6,multiThread-8,multiThread-10,multiThread-12,multiThread-14,multiThread-16,GPUnoLCA-64,GPUwithLCA-64,GPUnoLCA-128,GPUwithLCA-128,GPUnoLCA-256,GPUwithLCA-256,GPUnoLCA-512,GPUwithLCA-512,GPUnoLCA-1024,GPUwithLCA-1024\n";
			for (int i = 0; i < QueryTime_seq.size(); i++) {
				queryresult << 10000 + i * queryStep << "," << QueryTime_seq[i] << "," << QueryTime_MT_4[i] << ","
					<< QueryTime_MT_6[i] << "," << QueryTime_MT_8[i] << "," << QueryTime_MT_10[i] << "," << QueryTime_MT_12[i] << ","
					<< QueryTime_MT_14[i] << "," << QueryTime_MT_16[i] << ","
					<< QueryTime_noLCA_64[i] << "," << QueryTime_LCA_64[i] << ","
					<< QueryTime_noLCA_128[i] << "," << QueryTime_LCA_128[i] << ","
					<< QueryTime_noLCA_256[i] << "," << QueryTime_LCA_256[i] << ","
					<< QueryTime_noLCA_512[i] << "," << QueryTime_LCA_512[i] << ","
					<< QueryTime_noLCA_1024[i] << "," << QueryTime_LCA_1024[i] << "\n";
			}
			queryresult.close();
		}

	}

	void Graph_D_H::Graph::answerRealNYQuery_noHub()
	{
		uploadRealQuery();
		translateQuery();
		int endTime = 500000;
		bool seqContinue = true;
		bool MTContinue_4 = true;
		bool MTContinue_6 = true;
		bool MTContinue_8 = true;
		bool MTContinue_10 = true;
		bool MTContinue_12 = true;
		bool MTContinue_14 = true;
		bool MTContinue_16 = true;
		for (int i = 10000; i <= queryMaxSize; i += queryStep) {

			if (seqContinue) {
				long long int seqtime = H2Hquery_bunch_noHub(i);
				QueryTime_seq.push_back(seqtime);
				seqContinue = (seqtime < endTime);
			}
			else {
				QueryTime_seq.push_back(endTime);
			}

			if (MTContinue_4) {
				threadNumber = 4;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_4.push_back(seqtime);
				MTContinue_4 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_4.push_back(endTime);
			}

			if (MTContinue_6) {
				threadNumber = 6;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_6.push_back(seqtime);
				MTContinue_6 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_6.push_back(endTime);
			}

			if (MTContinue_8) {
				threadNumber = 8;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_8.push_back(seqtime);
				MTContinue_8 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_8.push_back(endTime);
			}

			if (MTContinue_10) {
				threadNumber = 10;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_10.push_back(seqtime);
				MTContinue_10 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_10.push_back(endTime);
			}

			if (MTContinue_12) {
				threadNumber = 12;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_12.push_back(seqtime);
				MTContinue_12 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_12.push_back(endTime);
			}

			if (MTContinue_14) {
				threadNumber = 14;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_14.push_back(seqtime);
				MTContinue_14 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_14.push_back(endTime);
			}

			if (MTContinue_16) {
				threadNumber = 16;
				long long int seqtime = H2Hquery_bunch_MultiThread_noHub(i);
				QueryTime_MT_16.push_back(seqtime);
				MTContinue_16 = (seqtime < endTime);
			}
			else {
				QueryTime_MT_16.push_back(endTime);
			}
			//QueryTime_seq.push_back(H2Hquery_bunch(i));
			//QueryTime_MT.push_back(H2Hquery_bunch_MultiThread(i));
			cudaThreadNum = 64;
			QueryTime_noLCA_64.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_64.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));

			cudaThreadNum = 128;
			QueryTime_noLCA_128.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_128.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));

			cudaThreadNum = 256;
			QueryTime_noLCA_256.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_256.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));

			cudaThreadNum = 512;
			QueryTime_noLCA_512.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_512.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));

			cudaThreadNum = 1024;
			QueryTime_noLCA_1024.push_back(H2Hquery_bunch_noLCA_D_noHub(i));
			QueryTime_LCA_1024.push_back(H2Hquery_bunch_UsingLCA_D_noHub(i));
		}
		string resultFile = "./QueryResult/" + graphName + "/PTHeight-" + to_string(TreeHeight) + "-ChangeHeight-" + to_string(changeHeight) + "-real.csv";
		fstream queryresult(resultFile, ios::app | ios::in | ios::out);
		if (queryresult.is_open()) {
			//queryresult << "CPUThreadNum,CUDAThreadNum\n";
			//queryresult <<threadNumber<<","<< cudaThreadNum << "\n";
			//if(queryresult.peek() == std::ifstream::traits_type::eof())
			queryresult << "queryNumber,sequence,multiThread-4,multiThread-6,multiThread-8,multiThread-10,multiThread-12,multiThread-14,multiThread-16,GPUnoLCA-64,GPUwithLCA-64,GPUnoLCA-128,GPUwithLCA-128,GPUnoLCA-256,GPUwithLCA-256,GPUnoLCA-512,GPUwithLCA-512,GPUnoLCA-1024,GPUwithLCA-1024\n";
			for (int i = 0; i < QueryTime_seq.size(); i++) {
				queryresult << 10000 + i * queryStep << "," << QueryTime_seq[i] << "," << QueryTime_MT_4[i] << ","
					<< QueryTime_MT_6[i] << "," << QueryTime_MT_8[i] << "," << QueryTime_MT_10[i] << "," << QueryTime_MT_12[i] << ","
					<< QueryTime_MT_14[i] << "," << QueryTime_MT_16[i] << ","
					<< QueryTime_noLCA_64[i] << "," << QueryTime_LCA_64[i] << ","
					<< QueryTime_noLCA_128[i] << "," << QueryTime_LCA_128[i] << ","
					<< QueryTime_noLCA_256[i] << "," << QueryTime_LCA_256[i] << ","
					<< QueryTime_noLCA_512[i] << "," << QueryTime_LCA_512[i] << ","
					<< QueryTime_noLCA_1024[i] << "," << QueryTime_LCA_1024[i] << "\n";
			}
			queryresult.close();
		}
	}
}