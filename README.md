# GPU-H2H
The code for H2H GPU-version.

## Environment

We run experiments on a server with i5-13600k 20 threads CPU, a NVIDIA RTX 4090 24GB GPU, 32GB memory, and running Ubuntu 22.04-LTS system. 
We use the nvcc compiler to compile our code. Our g++ version is 11.2, nvcc version is 12.0, and GPU driver version is 12.0.525. GPU compute capacity is 6.0.

Although the Thrust library only requires a compute capacity of 2.0 or above, in order to exploit the GPU performance, please use a NVIDIA GPU with at least 5.0 or above compute capacity.

Some devices may need gcc_10 or lower compiler version: '-ccbin /usr/bin/gcc-10', and choose & link standard liberay: '-std=c++14' and '-lstdc++', '-lm'.

## Compile

Please run 
```./Zcompile.sh```
to compile our code.


## Data preparation

Graph files for testing like NY, FLA, E can be found at http://www.diag.uniroma1.it/~challenge9/download.shtml.
Directory ./Graph/USA covers a test network 13. 

The format of Edge file Edge.txt is:
```
13 44               //Node Number, Edge Number;
a 1 2 4             //Edge (1, 2) with weight = 4;
......
```

The format of Node file Node.txt is:
```
13                      //Node Number;
v 1 2.0 10.1
v 2 4.0 10.2            //vertex 2 and it's (Longitude, Latitude) = (4.0, 10.2)
......
```

We also use partitioning to generate orders. Directory ./Partition/13 covers two partition result. 
Such partitions are the results of https://github.com/sauccjy/PartitionForGPUH2H.git. 
The results in our paper are obtained by using latitude-longitude partitioning, but not minimum cut.

The partition results can be arbitrary as long as it meets the following format requirements.

We store there partition result as an Array. ''xxxID.txt'' records such array, and ''xxxCut.txt'' records the partition tree.

The format of ''xxxID.txt' (latitudeID.txt) is:
```
13                  //Node Number;
0 1.0 6.0 2         //index in array, Longitude, Latitude, vertex ID - 1;
1 6.0 1.1 6
......
```

The format of ''xxxCut.txt'' (latitudeCut.txt) is:
```
4                       //Partition Tree height
1 0 13                  //Sub graph at top-1 layer, [start index, end index) = [0, 13)
2 0 6                   //Sub graph at top-2 layer, [start index, end index) = [0, 6)
2 6 13
3 0 3
......
```

## Run
After compiling, we run next code to construct Label.
```
./G2H GraphName PartitionHeight ChangeToGlobalHeight threadNum partitionMethod useGPUContract parallelFine useGPUConstruct
```

'GraphName' represents the graph;

'Partition Height' represents the partition Tree Height you used;

'ChangeToGlobalHeight' represents the layer switch to CPU-Global Contraction;

'threadNum' represents max thread number on CPU;

'partitionMethod' represents partitioning way (only read files in ./Partition/'GraphName'/): 
1) if == 1, use XYCoord partitioning;
2) if == 2, use HC2L-like Partitioning;
3) if == 3, use METIS; 
4) if == 4, use SCOTCH;
   *Notice: there is no partition code in this link. Please copy the partitioning results to ./Partition/'GraphName'/;

'useGPUContract' is a bool type, if useGPUContract == 0, then contract on CPU, and == 1 on GPU-hybrid;

'parallelFine' is a int type, 
1) if == 0, construct labels Sequentially;
2) if == 1, PF method;
3) if == 2, LL method; 

'useGPUConstruct' is a bool type too, if useGPUConstruct == 0, then construct on CPU, and == 1 on GPU;


For example, 
```
./G2H 13 4 1 20 1 1 2 1
```
will contract 4 - 1 layers on GPU but rest on CPU,  and construct labels on GPU with LL frontier method. 

After the end of construction, construct information will be recorded at './ConstructInfo/GraphName/'.
If answered query, query time will be recorded at './QueryResult/graphName/'.


## Difference between L3-Label and L4-Label

We did not comment out the L4 code.

If readers need to run L3 label, please comment out L4 code in next functions:

Kernel_functions.cu: 

(line 738-742) makeH2HLabel_noCommunication_D_noHub_3();

(line 828-837) makeH2HLabel_noCommunication_D_noHub();
                     
H2HConstruction.cu: 

(line 723-729) makeH2HLabel_noHub_serial();

(line 791-797) makeH2HLabel_noHub_multiThred();

The query functions for answering L3 and L4 are defined in
```
Query.cu::H2HDistancQuery_UsingLCA_noHub()  //L4
Query.cu::H2HDistancQuery_noLCA_noHub()     //L3
```
