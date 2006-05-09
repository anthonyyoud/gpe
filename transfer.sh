#!/bin/sh

EXE=ulimit_hack.sh
DATA="gpe parameters.f90 run_script.sh"
CUR_DIR=`pwd`
RUN_DIR=`echo $CUR_DIR | cut -c29-`
TAL_DIR=/local/ajy25/$RUN_DIR

echo Directory on talisman will be $TAL_DIR

ssh talisman "mkdir -p $TAL_DIR"
if [ $? -eq 0 ]; then
  scp $EXE $DATA talisman:$TAL_DIR
else
  echo ERROR:  Failed to make directory on talisman
  exit 1
fi
