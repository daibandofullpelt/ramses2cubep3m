#!/bin/bash
#@ wall_clock_limit = 48:00:00
#@ job_type = MPICH
#@ job_name = runme_all_test
#@ class = general
#@ node = 128
#@ tasks_per_node = 16
#@ node_usage = not_shared
#@ network.MPI = sn_all,not_shared,us
#@ output = $(job_name).$(schedd_host).$(jobid).out
#@ error = $(job_name).$(schedd_host).$(jobid).err
#@ queue
#@ notification=always
#@ notify_user=aurel.schneider@sussex.ac.uk

. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.intel/4.1.0_gcc

#automatic export
set -a
#Uncomment for Debugging
set -xv

cd /home/hpc/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2

echo 15.132>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_15.132/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_15.132/AHF-v1.0-056 ./AHF_chunk/z_15.132/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_15.132/*
echo 14.699>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_14.699/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_14.699/AHF-v1.0-056 ./AHF_chunk/z_14.699/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_14.699/*
echo 14.294>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_14.294/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_14.294/AHF-v1.0-056 ./AHF_chunk/z_14.294/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_14.294/*
echo 13.914>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_13.914/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_13.914/AHF-v1.0-056 ./AHF_chunk/z_13.914/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_13.914/*
echo 13.557>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_13.557/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_13.557/AHF-v1.0-056 ./AHF_chunk/z_13.557/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_13.557/*
echo 13.221>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_13.221/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_13.221/AHF-v1.0-056 ./AHF_chunk/z_13.221/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_13.221/*
echo 12.903>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_12.903/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_12.903/AHF-v1.0-056 ./AHF_chunk/z_12.903/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_12.903/*
echo 12.603>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_12.603/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_12.603/AHF-v1.0-056 ./AHF_chunk/z_12.603/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_12.603/*
echo 12.318>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_12.318/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_12.318/AHF-v1.0-056 ./AHF_chunk/z_12.318/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_12.318/*
echo 12.048>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_12.048/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_12.048/AHF-v1.0-056 ./AHF_chunk/z_12.048/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_12.048/*
echo 11.791>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_11.791/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_11.791/AHF-v1.0-056 ./AHF_chunk/z_11.791/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_11.791/*
echo 11.546>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_11.546/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_11.546/AHF-v1.0-056 ./AHF_chunk/z_11.546/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_11.546/*
echo 11.313>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_11.313/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_11.313/AHF-v1.0-056 ./AHF_chunk/z_11.313/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_11.313/*
echo 11.090>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_11.090/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_11.090/AHF-v1.0-056 ./AHF_chunk/z_11.090/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_11.090/*
echo 10.877>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_10.877/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_10.877/AHF-v1.0-056 ./AHF_chunk/z_10.877/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_10.877/*
echo 10.673>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_10.673/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_10.673/AHF-v1.0-056 ./AHF_chunk/z_10.673/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_10.673/*
echo 10.478>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_10.478/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_10.478/AHF-v1.0-056 ./AHF_chunk/z_10.478/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_10.478/*
echo 10.290>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_10.290/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_10.290/AHF-v1.0-056 ./AHF_chunk/z_10.290/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_10.290/*
echo 10.110>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_10.110/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_10.110/AHF-v1.0-056 ./AHF_chunk/z_10.110/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_10.110/*
echo 9.938>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_9.938/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_9.938/AHF-v1.0-056 ./AHF_chunk/z_9.938/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_9.938/*
echo 9.771>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_9.771/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_9.771/AHF-v1.0-056 ./AHF_chunk/z_9.771/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_9.771/*
echo 9.611>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_9.611/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_9.611/AHF-v1.0-056 ./AHF_chunk/z_9.611/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_9.611/*
echo 9.457>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_9.457/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_9.457/AHF-v1.0-056 ./AHF_chunk/z_9.457/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_9.457/*
echo 9.308>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_9.308/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_9.308/AHF-v1.0-056 ./AHF_chunk/z_9.308/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_9.308/*
echo 9.164>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_9.164/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_9.164/AHF-v1.0-056 ./AHF_chunk/z_9.164/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_9.164/*
echo 9.026>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_9.026/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_9.026/AHF-v1.0-056 ./AHF_chunk/z_9.026/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_9.026/*
echo 8.892>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_8.892/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_8.892/AHF-v1.0-056 ./AHF_chunk/z_8.892/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_8.892/*
echo 8.762>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_8.762/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_8.762/AHF-v1.0-056 ./AHF_chunk/z_8.762/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_8.762/*
echo 8.636>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_8.636/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_8.636/AHF-v1.0-056 ./AHF_chunk/z_8.636/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_8.636/*
echo 8.515>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_8.515/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_8.515/AHF-v1.0-056 ./AHF_chunk/z_8.515/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_8.515/*
echo 8.397>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_8.397/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_8.397/AHF-v1.0-056 ./AHF_chunk/z_8.397/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_8.397/*
echo 8.283>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_8.283/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_8.283/AHF-v1.0-056 ./AHF_chunk/z_8.283/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_8.283/*
echo 8.172>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_8.172/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_8.172/AHF-v1.0-056 ./AHF_chunk/z_8.172/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_8.172/*
echo 8.064>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_8.064/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_8.064/AHF-v1.0-056 ./AHF_chunk/z_8.064/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_8.064/*
echo 7.960>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.960/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.960/AHF-v1.0-056 ./AHF_chunk/z_7.960/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.960/*
echo 7.859>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.859/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.859/AHF-v1.0-056 ./AHF_chunk/z_7.859/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.859/*
echo 7.760>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.760/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.760/AHF-v1.0-056 ./AHF_chunk/z_7.760/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.760/*
echo 7.664>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.664/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.664/AHF-v1.0-056 ./AHF_chunk/z_7.664/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.664/*
echo 7.570>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.570/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.570/AHF-v1.0-056 ./AHF_chunk/z_7.570/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.570/*
echo 7.480>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.480/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.480/AHF-v1.0-056 ./AHF_chunk/z_7.480/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.480/*
echo 7.391>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.391/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.391/AHF-v1.0-056 ./AHF_chunk/z_7.391/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.391/*
echo 7.305>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.305/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.305/AHF-v1.0-056 ./AHF_chunk/z_7.305/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.305/*
echo 7.221>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.221/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.221/AHF-v1.0-056 ./AHF_chunk/z_7.221/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.221/*
echo 7.139>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.139/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.139/AHF-v1.0-056 ./AHF_chunk/z_7.139/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.139/*
echo 7.059>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_7.059/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_7.059/AHF-v1.0-056 ./AHF_chunk/z_7.059/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_7.059/*
echo 6.981>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_6.981/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_6.981/AHF-v1.0-056 ./AHF_chunk/z_6.981/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_6.981/*
echo 6.905>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_6.905/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_6.905/AHF-v1.0-056 ./AHF_chunk/z_6.905/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_6.905/*
echo 6.830>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_6.830/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_6.830/AHF-v1.0-056 ./AHF_chunk/z_6.830/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_6.830/*
echo 6.757>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_6.757/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_6.757/AHF-v1.0-056 ./AHF_chunk/z_6.757/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_6.757/*
echo 6.686>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_6.686/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_6.686/AHF-v1.0-056 ./AHF_chunk/z_6.686/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_6.686/*
echo 6.617>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_6.617/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_6.617/AHF-v1.0-056 ./AHF_chunk/z_6.617/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_6.617/*
echo 6.549>redshift_checkpoint
mpiexec -n 512 ./chunk_cubep3m/chunk chunking_paramters.dat

rm  mfile.*

export SUBJOB
export NUMBER_OF_SUBJOBS=64
export SIZE_OF_SUBJOB=32
export EXE=./AHF_chunk/z_6.549/AHF-v1.0-056

for X in `seq 0 $(($NUMBER_OF_SUBJOBS-1))`

do
   SUBJOB=$X
   L1=$(( $X*SIZE_OF_SUBJOB + 1 ))
   L2=$(( $X*SIZE_OF_SUBJOB + $SIZE_OF_SUBJOB  ))
   #cut appropriate lines from Hostfile
   cp $LOADL_HOSTFILE mfile
   sed -n -e "$L1,${L2}p" $LOADL_HOSTFILE >mfile.$X
   #start subjobs with different hostfilms, pinning is done in subjob.intel
   I_MPI_PIN=yes
   I_MPI_PIN_CELL=core
   I_MPI_HOSTFILE=`pwd`/mfile.$X
   KMP_AFFINITY="granularity=fine,compact,1,0"
   mpiexec -n $SIZE_OF_SUBJOB -f mfile.$X ./AHF_chunk/z_6.549/AHF-v1.0-056 ./AHF_chunk/z_6.549/chunk_"$SUBJOB".input  &
done
wait
rm -r /gpfs/scratch/pr86be/pr3df064/AHF_analysis/cubepm_130708_8_2048_64Mpc_ext2/chunked_output/z_6.549/*

