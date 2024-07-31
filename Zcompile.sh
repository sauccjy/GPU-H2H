#! /bin/bash

nvcc -O3 -diag-suppress 20012 -arch=sm_60 -o G2H ./*.cu