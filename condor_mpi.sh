universe = parallel
#executable = openmpi_script
executable = /quantum/bin/condor-parallel.sh
#executable = run_script_quantum.sh
arguments = run_script_quantum.sh
#arguments = fast2.5D_MPI
machine_count = 8
log = logfile
output = outfile.$(NODE)
error = errfile.$(NODE)
notification = Never
notify_user = a.j.youd@ncl.ac.uk
should_transfer_files = yes
when_to_transfer_output = on_exit
#transfer_input_files = run_script_quantum.sh, hello_mpi
#transfer_input_files = run_script_quantum.sh, fast2.5D_MPI
transfer_input_files = run_script_quantum.sh, gpe
#transfer_input_files = fast2.5D_MPI
#transfer_output_files = quantum
environment="JOBID=$(Cluster).$(Process)"
queue
