#!/bin/sh
#PBS -S /bin/sh
#PBS -j oe -o run.log

TALISMAN=master
NODED=noded
GIGA=giga
HOST=`hostname -s` && HOST=${HOST/[0-9]*}

if [ $HOST == $GIGA ]; then
  cd $PBS_O_WORKDIR
  #uniq $PBS_NODEFILE
fi

#****************************************************************************
#PARAMETERS
#****************************************************************************
NPROCS=16
EXE=ulimit_hack.sh #gpe
DATA='parameters.f90 gpe'
DATA_DIR=`pwd`
PROC_DIR=proc
case $HOST in
  $GIGA)
    RUN_DIR=/work/ay_$PBS_JOBID
    SSH=/usr/local/bin/ssh
    SCP=/usr/local/bin/scp
    HOSTS=`uniq $PBS_NODEFILE`
    HOSTFILE=mpd.hosts
  ;;
  $TALISMAN)
    RUN_DIR=`pwd`
    EXE=ulimit_hack.sh
    DATA='parameters.f90 gpe'
    HOSTFILE=mpd.hosts
  ;;
  $NODED)
    RUN_DIR=`pwd`
    EXE=ulimit_hack.sh
    DATA='parameters.f90 gpe'
    HOSTFILE=mpd.hosts.64
  ;;
esac

#****************************************************************************
#NO CHANGES NECESSARY BELOW HERE
#****************************************************************************
case $HOST in
  $GIGA)
    mkdir $RUN_DIR
    for i in `seq 0 $(( $NPROCS-1 ))`
    do
      if [ $i -lt 10 ]; then
        mkdir $RUN_DIR/${PROC_DIR}0$i
      else
        mkdir $RUN_DIR/$PROC_DIR$i
      fi
    done
    cp $EXE $RUN_DIR
    cp $DATA $RUN_DIR
    cp -rf $PROC_DIR* $RUN_DIR
    cd $RUN_DIR
    
    echo No. processors = $NPROCS
    for NODE in $HOSTS
    do
      if [ `hostname` != $NODE ]; then
        $SSH $NODE "mkdir $RUN_DIR"
        $SCP -r $PROC_DIR* $EXE $DATA $NODE:$RUN_DIR 2> /dev/null
      fi
    done
    uniq $PBS_NODEFILE > $HOSTFILE
    NUMHOSTS=`cat $HOSTFILE | wc -l`
    #mpdboot -n $NUMHOSTS --ncpus=2 --rsh=$SSH
    #mpdtrace
    time mpiexec -l -n $NPROCS $EXE
    #mpdallexit
    rm $EXE $HOSTFILE
    
    TARFILE=data.tar
    cd $RUN_DIR
    tar cvf $TARFILE * &> /dev/null
    cp $TARFILE $DATA_DIR
    cd $DATA_DIR && tar xvf $TARFILE &> /dev/null
    rm $TARFILE
    for NODE in $HOSTS
    do
      if [ `hostname` != $NODE ]; then
        $SSH $NODE "cd $RUN_DIR && \
                    tar cvf $TARFILE * &> /dev/null && \
                    cp $TARFILE $DATA_DIR && \
                    cd $DATA_DIR && tar xvf $TARFILE &> /dev/null && \
                    rm $TARFILE"
        $SSH $NODE "rm -rf $RUN_DIR"
      fi
    done
    rm -rf $RUN_DIR
  ;;
  $TALISMAN|$NODED)
    for i in `seq 0 $(( $NPROCS-1 ))`
    do
      if [ $i -lt 10 ]; then
        mkdir ${PROC_DIR}0$i
      else
        mkdir $PROC_DIR$i
      fi
    done
    
    echo Making directories and copying files...
    for node in `cat ${HOME}/$HOSTFILE`
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
    if [ $HOST = $TALISMAN ]; then
      time mpiexec -1 -l -n $NPROCS $EXE
    else
      time mpiexec -l -n $NPROCS $EXE
    fi
    echo Done
    
    echo Tarring up directories and copying back to master...
    TARFILE=dirs.tar
    for node in `cat ${HOME}/$HOSTFILE | grep -v noded1`
    do
      ssh $node "cd $RUN_DIR && tar cvf $TARFILE $PROC_DIR* &> /dev/null \
                 && rm -r $PROC_DIR*"
      if [ $? -eq 0 ]; then
        scp $node:$RUN_DIR/* $RUN_DIR &> /dev/null && \
        tar xvf $TARFILE &> /dev/null && \
        rm $TARFILE && \
        ssh $node "rm -r $RUN_DIR"
      else
        echo "ERROR: tar failed on node $node"
        exit 1
      fi
    done
    echo Done
  ;;
esac
