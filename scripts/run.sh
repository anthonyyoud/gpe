#!/bin/bash
#$ -cwd -V
#$ -e error_log.$JOB_ID
#$ -o output_log.$JOB_ID
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

DOMAINNAME=$(dnsdomainname)
# Node list
POLARIS="polaris.leeds.ac.uk"

case "$DOMAINNAME" in
    $POLARIS);;
    *) PATH=.:/bin:/usr/bin;;
esac

usage() {
case "$DOMAINNAME" in
    $POLARIS)
    DESCRIPTION="
        Submit a GPE run to the polaris queue either from the start or by doing
        a restart from an existing run."
    USAGE="
        qsub [QSUB_OPTIONS] run-polaris.sh [OPTIONS]"
    QSUB_OPTIONS="
QSUB_OPTIONS:

        These are the options to pass to the qsub command.  See man qsub
        for more information.
    
        To request one hour of run time on four cores, for example, specify
        -l h_rt=1:00:00 and -pe ib 4."
;;
    *)
    DESCRIPTION="
        Run the GPE code on <NPROCS> processes either from the start or by
        doing a restart from an existing run."
    USAGE="
        run.sh [OPTIONS] <NPROCS>"
        
;;
esac

cat <<EOF
run.sh
---------

DESCRIPTION:
$DESCRIPTION

USAGE:
$USAGE

$QSUB_OPTIONS

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
case "$DOMAINNAME" in
    $POLARIS)
        mpiexec $EXE
    ;;
    *)
        nohup /usr/bin/time -p mpiexec -n $NPROCS $EXE &> $LOGFILE &
    ;;
esac
}

# Get command line options
options=`getopt -o r:cph -l restart:,clean,post-process,help -n run.sh -- "$@"`

case "$DOMAINNAME" in
    $POLARIS)
        if [ "BATCH" != "$ENVIRONMENT" ]; then
            usage
            exit 1
        fi
    ;;
    *)
        # If no options, show the help
        if [ $# == 0 ]; then
          usage
          exit 1
        fi
    ;;
esac

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

case "$DOMAINNAME" in
    $POLARIS)
        NPROCS=$NSLOTS
    ;;
    *)
        NPROCS=$1
    ;;
esac

DIGITS=4
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
