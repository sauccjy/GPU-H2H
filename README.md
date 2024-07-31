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

