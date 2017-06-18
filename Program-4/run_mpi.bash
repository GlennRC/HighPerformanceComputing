#! /bin/bash/

cd /root/rcglenn/repos/ecpe893_hpc_spring2016/Program-4
make clean
make
echo "times" > parallel_time.csv
for i in {1..1}
do
	for j in 2 4 6 8
	do
		for k in 256 512 1024 2048
		do
			mpiexec -n $j ./exec /root/rcglenn/repos/ecpe893_hpc_spring2016/Lenna_Images/Lenna_org_$k.pgm 1.0 >> parallel_time.csv
		done
	done
done	
