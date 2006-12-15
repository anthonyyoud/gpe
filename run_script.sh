#!/bin/sh
#PBS -S /bin/sh
#PBS -j oe -o run.log

FILTER=0
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
if [ $FILTER -eq 1 ]; then
  DATA='parameters.f90 gpe p_saved.dat'
else
  DATA='parameters.f90 gpe'
fi
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
    if [ ! -e ${PROC_DIR}00 ]; then
      for i in `seq 0 $(( $NPROCS-1 ))`
      do
        if [ $i -lt 10 ]; then
          mkdir ${PROC_DIR}0$i
        else
          mkdir $PROC_DIR$i
        fi
      done
    fi
    TARFILE=data.tar
    cp $EXE $RUN_DIR
    cp $DATA $RUN_DIR
    tar hcf $TARFILE $PROC_DIR* 2> /dev/null
    cp $TARFILE $RUN_DIR
    rm -rf $PROC_DIR*
    cd $RUN_DIR
    tar xf $TARFILE
    rm $TARFILE
    cd $DATA_DIR
    
    echo No. processors = $NPROCS
    for NODE in $HOSTS
    do
      if [ `hostname` != $NODE ]; then
        $SSH $NODE "mkdir $RUN_DIR"
        $SSH $NODE "cd $DATA_DIR && cp $TARFILE $EXE $DATA $RUN_DIR"
        $SSH $NODE "cd $RUN_DIR && tar xf $TARFILE && rm $TARFILE"
      fi
    done
    rm $TARFILE
    cd $RUN_DIR
    uniq $PBS_NODEFILE > $HOSTFILE
    NUMHOSTS=`cat $HOSTFILE | wc -l`
    #mpdboot -n $NUMHOSTS --ncpus=2 --rsh=$SSH
    #mpdtrace
    time mpiexec -l -n $NPROCS $EXE
    #mpdallexit
    rm $EXE $HOSTFILE
    
    for dir in $PROC_DIR*
    do
      if [ $FILTER -eq 1 ]; then
        rm $dir/dens*
      fi
      num_files=$(ls -1 $dir | wc -l)
      if [ $num_files -eq 1 ]; then
        rm -r $dir
      fi
    done
    tar cf $TARFILE * &> /dev/null
    cp $TARFILE $DATA_DIR
    cd $DATA_DIR && tar xf $TARFILE &> /dev/null && rm $TARFILE
    for NODE in $HOSTS
    do
      if [ `hostname` != $NODE ]; then
        $SSH $NODE "cd $RUN_DIR && \
                    for dir in $PROC_DIR*
                    do
                      if [ $FILTER -eq 1 ]; then
                        rm \$dir/dens*
                      fi
                      num_files=\$(ls \$dir | wc -l)
                      if [ \$num_files -eq 1 ]; then
                        rm -r \$dir
                      fi
                    done && \
                    tar cf $TARFILE * &> /dev/null && \
                    cp $TARFILE $DATA_DIR && \
                    cd $DATA_DIR && \
                    tar xf $TARFILE &> /dev/null && \
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
      ssh $node "cd $RUN_DIR && tar cf $TARFILE $PROC_DIR* &> /dev/null \
                 && rm -r $PROC_DIR*"
      if [ $? -eq 0 ]; then
        scp $node:$RUN_DIR/* $RUN_DIR &> /dev/null && \
        tar xf $TARFILE &> /dev/null && \
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
