# GPU-H2H
The code for H2H GPU-version.

##Environment
We run experiments on a server with i5-13600k 20 threads CPU, a NVIDIA RTX 4090 24GB GPU, 32GB memory, and running Ubuntu 22.04-LTS system. 
We use the nvcc compiler to compile our code. Our g++ version is 11.2, nvcc version is 12.0, and GPU driver version is 12.0.525. GPU compute capacity is 6.0.

Although the Thrust library only requires a compute capacity of 2.0 or above, in order to exploit the GPU performance, please use a compute capacity of at least 5.0 or above.


##Compile
Please run 
```./Zcompile.sh```
or
```nvcc -O3 -arch=sm_60 -o G2H ./*.cu```
to compile our code.


##Data preparation

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
We store there partition result as an Array. ''xxxID.txt'' records such array, and ''xxxCut.txt'' records the partition tree.

The format of ''xxxID.txt' is:
```
13                  //Node Number;
0 1.0 6.0 2         //index in array, Longitude, Latitude, vertex ID - 1;
1 6.0 1.1 6
......
```

The format of ''xxxCut.txt'' is:
```
4                       //Partition Tree height
1 0 13                  //Sub graph at top-1 layer, [start index, end index) = [0, 13)
2 0 6                   //Sub graph at top-2 layer, [start index, end index) = [0, 6)
2 6 13
3 0 3
......
```

##Run
After compiling, we run
```
./G2H 
```


After the end of construct or query, datas will be recorded at ./ConstructInfo/13/ and 
