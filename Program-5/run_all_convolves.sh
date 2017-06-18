#1 /bin/bash

cd /Users/glenncontreras/workspace/ubuntu/repos/ecpe893_hpc_spring2016/Program-5/
scp convolve_L2.cu gcontrer@ecs-cluster.serv.pacific.edu:Program-5/convolve_L2.cu
ssh gcontrer@ecs-cluster.serv.pacific.edu 'sh runconvolve1.sh'
scp gcontrer@ecs-cluster.serv.pacific.edu:/exports/home/gcontrer/Program-5/horizontal_gradient.pgm .
scp gcontrer@ecs-cluster.serv.pacific.edu:/exports/home/gcontrer/Program-5/vertical_gradient.pgm .
open horizontal_gradient.pgm
open vertical_gradient.pgm
