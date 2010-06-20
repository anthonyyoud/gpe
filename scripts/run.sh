#!/bin/bash
# $Id$
#----------------------------------------------------------------------------

# Simple script to run a job.

PATH=.:/bin:/usr/bin

usage() {
cat <<EOF
run.sh
---------

DESCRIPTION:

Run the GPE code on <NPROCS> processes either from the start or by doing a
restart from an existing run.

USAGE:

        run.sh [OPTIONS] <NPROCS>

OPTIONS:

        -r <directory>|--restart=<directory>
                    Do a restart of an existing run.  The directory where the
                    previous run exists must be given as an argument.  The
                    end_state.dat files are then copied to the new run
                    directory.

        -c|--clean
                    If the script detects the existence of process directories
                    in the current directory, new process directories are not
                    created and the script aborts.  Use -c to clean the current
                    directory, removing all files and directories except those
                    needed for a fresh run.

        -p|--post-process
                    Should be used when post-processing isosurfaces.

        -h|--help
                    Show usage information.
EOF
}

run() {
  echo Running job...
  nohup /usr/bin/time -p mpiexec -n $NPROCS $EXE &> $LOGFILE &
}

# Get command line options
options=`getopt -o r:cph -l restart:,clean,post-process,help -n run.sh -- "$@"`

# If no options, show the help
if [ $# == 0 ]; then
  usage
  exit 1
fi

eval set -- "$options"

CLEAN=0
POST=0
# Set parameters depending on options
while true
do
  case "$1" in
    -r|--restart) PREDIR=$2; shift 2;;
    -c|--clean) CLEAN=1; shift;;
    -p|--post-process) POST=1; shift;;
    -h|--help) usage; exit 1;;
    --) shift ; break ;;
    *) echo "Invalid flag" ; exit 1 ;;
  esac
done

if [ $# -eq 0 ];
then
  usage
  exit 1
fi

NPROCS=$1
DIGITS=2
EXE=gpe
LOGFILE=log.txt

NYPROCS=`grep nyprocs parameters.in | cut -d= -f 2`
NZPROCS=`grep nzprocs parameters.in | cut -d= -f 2`

if [ ! $((NYPROCS*NZPROCS)) -eq $NPROCS ]; then
  echo "ERROR: Number of processes ($NPROCS) does not match nyprocs*nzprocs in parameters.in."
  exit 1
fi

echo Going to run on $NPROCS processes

if [ $POST -eq 1 ]; then
  # If post-processing just do the run.
  echo Post-processing isosurfaces...
  run
  exit 0
fi

for i in `seq -f "%0${DIGITS}g" 0 $(($NPROCS-1))`
do
  PROCDIR=proc$i
  if [ -d $PROCDIR ]; then
    if [ $CLEAN -eq 0 ]; then
      echo WARNING: At least one process directory already exists - aborting.
      echo Use -c to clean current directory and begin a fresh run.
      exit 1
    else
      rm -rf proc* *.dat *.txt idf pdf spectrum vcf RUNNING ERROR links
    fi
  fi
  mkdir $PROCDIR
done

if [ ! -z $PREDIR ]; then
  echo Going to do a restart from $PREDIR
  for i in `seq -f "%0${DIGITS}g" 0 $(($NPROCS-1))`
  do
    PROCDIR=proc$i
    cp $PREDIR/$PROCDIR/end_state.dat end_state$i.dat
  done
fi

# Make directories in which to hold the PDFs, VCFs, spectra, and IDF.
mkdir pdf vcf spectrum idf &> /dev/null

# Run the job.
run
