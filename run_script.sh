#!/bin/sh

NPROCS=6
EXE=gpe
DATA=parameters.f90
RUN_DIR=`pwd`
PROC_DIR=proc

for i in `seq 0 $(( $NPROCS-1 ))`
do
  if [ $i -lt 10 ]; then
    mkdir ${PROC_DIR}0$i
  else
    mkdir $PROC_DIR$i
  fi
done

echo Making directories and copying files...
for node in `cat ${HOME}/mpd.hosts`
do
  ssh $node "mkdir -p $RUN_DIR"
  if [ $? -eq 0 ]; then
    scp -r $PROC_DIR* $EXE $DATA $node:$RUN_DIR 2> /dev/null
  else
    echo "ERROR: Make directory failed on node $node"
    exit 1
  fi
done
echo Done

echo Starting job...
time mpiexec -1 -l -n $NPROCS $EXE
echo Done

echo Tarring up directories and copying back to master...
for node in `cat ${HOME}/mpd.hosts`
do
  ssh $node "cd $RUN_DIR && tar cvf dirs.tar $PROC_DIR* &> /dev/null \
             && rm -r $PROC_DIR*"
  if [ $? -eq 0 ]; then
    scp $node:$RUN_DIR/* $RUN_DIR &> /dev/null && \
    tar xvf dirs.tar &> /dev/null && \
    rm dirs.tar && \
    ssh $node "rm -r $RUN_DIR"
  else
    echo "ERROR: tar failed on node $node"
    exit 1
  fi
done
echo Done
