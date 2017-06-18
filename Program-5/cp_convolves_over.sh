#! /bin/bash

cd /Users/glenncontreras/workspace/ubuntu/repos/ecpe893_hpc_spring2016/Program-5/

scp convolve_L0.cu gcontrer@ecs-cluster.serv.pacific.edu:Program-5/convolve_L0.cu
scp convolve_L1.cu gcontrer@ecs-cluster.serv.pacific.edu:Program-5/convolve_L1.cu
scp convolve_L2.cu gcontrer@ecs-cluster.serv.pacific.edu:Program-5/convolve_L2.cu

#ssh gcontrer@ecs-cluster.serv.pacific.edu 'sh runconvolve.sh'
