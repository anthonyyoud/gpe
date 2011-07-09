#!/bin/sh
#Copyright 2011 Anthony Youd/Newcastle University
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


FILTER=0
TARFILE=data_`hostname`.tar
END_STATE=end_state.tar

#****************************************************************************
#PARAMETERS
#****************************************************************************
NPROCS=$_CONDOR_NPROCS
echo `hostname` $NPROCS
#EXE=hello_mpi
#EXE=fast2.5D_MPI
EXE=gpe
PROC_DIR=proc

#****************************************************************************
#NO CHANGES NECESSARY BELOW HERE
#****************************************************************************
if [ -e $END_STATE ]; then
  tar xf $END_STATE
  rm $END_STATE
fi

for i in `seq 0 $(( $NPROCS-1 ))`
do
  if [ $i -lt 10 ]; then
    mkdir ${PROC_DIR}0$i
  else
    mkdir $PROC_DIR$i
  fi
done

touch modified
start_time=`ls -cl --time-style=+%s -d modified | cut -d" " -f6`
rm modified

echo No. processors = $NPROCS
time $EXE

rm $EXE
rm -f end_state??.dat

for dir in ${PROC_DIR}*
do
  if [ $FILTER -eq 1 ]; then
    rm $dir/dens*
  fi
  mod_time=`ls -cl --time-style=+%s -d $dir | cut -d" " -f6`
  if [ $mod_time -le $start_time ]; then
    rm -r $dir
  fi
done
#tar cf $TARFILE --exclude="chirp.config" --exclude="condor_exec.exe" --exclude="errfile.*" --exclude="outfile.*" --exclude="$0" * &> /dev/null
tar cf $TARFILE $PROC_DIR* &> /dev/null
hostname
echo "----------"
pwd
ls -l --time-style=+%s
echo NPROCS $_CONDOR_NPROCS
echo "----------"
echo "----------"
