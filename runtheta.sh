#!/bin/bash -evx
####COBALT -q flat-quad
#COBALT -t 60
#COBALT -A Performance ###EarlyPerf_theta 

EXE=./testcase

echo "Running Cobalt Job $COBALT_JOBID on $COBALT_PARTNAME."

#io tracing
export DXT_ENABLE_IO_TRACE=4

#mpich_abort_on_rw_wrror 
export ATP_ENABLED=1

export MPICH_MPIIO_HINTS_DISPLAY=1
export MPICH_MPIIO_AGGREGATOR_PLACEMENT_DISPLAY=1 
export MPICH_MPIIO_STATS=1
export MPICH_MPIIO_XSTATS=1
export MPICH_MPIIO_TIMERS=1
export MPICH_MPIIO_ABORT_ON_RW_ERROR=enable

export MPICH_CPUMASK_DISPLAY=1
export MPICH_RANK_REORDER_DISPLAY=1
export MPICH_RANK_REORDER_DISPLAY=true
export MPICH_CPUMASK_DISPLAY=true
export KMP_AFFINITY=verbose
export MPICH_VERSION_DISPLAY=1
export MPICH_NEMESIS_ASYNC_PROGRESS=1
export MPICH_ENV_DISPLAY=1
export MPICH_MPIIO_HINTS="*:cb_nodes=2"

locfile=loc_${COBALT_JOBSIZE}_${COBALT_JOBID}.txt
echo ${COBALT_PARTNAME} > $locfile 
jobmapfile=jobmap_${COBALT_JOBSIZE}_${COBALT_JOBID}.txt
#python parsejobnodes.py theta.computenodes $locfile > $jobmapfile

aprun -n 1 -N 1 -d 1 -j 1 -r 1 ./location.x
#aprun -n 1 -N 1 -d 1 -j 1 -r 1 ./info

THREADS=1
#export PAT_RT_SUMMARY=0

startiter=$1
enditer=$(($startiter+1))

ENVVARS="" #"--env MPICH_MPIIO_HINTS "*:cb_nodes=4""
#--env MPICH_ENV_DISPLAY=1 --env MPICH_VERSION_DISPLAY=1 --env MPICH_SMP_SINGLE_COPY_SIZE=1024 --env MPICH_NEMESIS_ASYNC_PROGRESS=1 --env MPICH_SHARED_MEM_COLL_OPT=1 --env MPICH_GNI_ASYNC_PROGRESS_STATS=enabled -r 1"

for iter in `seq $startiter $enditer` 
do
 for ppn in 64
 do
  RANKS=$((${COBALT_PARTSIZE}*$ppn))
	echo 
  APRUNPARAMS=" -n ${RANKS} -N ${ppn} -d 1 -j 1 -r 1 " #--attrs mcdram=cache:numa=quad "
  for STRIPECNT in 8 # 16 24 32 
  do 
   for STRIPESZ in 8M #2M  
   do
    for size in 32 1024 4096
    do
      OUTPUT=output_${COBALT_PARTSIZE}_${RANKS}_R${ppn}_${STRIPECNT}_${STRIPESZ}_${size}_${iter}_${COBALT_JOBID}

      ARG=" $size 0 0 1 0"
	    echo "Starting $OUTPUT with $ARG"
      #mkdir pat_${OUTPUT}
      #export PAT_RT_EXPFILE_DIR=pat_${OUTPUT}

      FNAME=TestFile-${RANKS}
      rm -f $FNAME 2>/dev/null
      echo "Testing: echo $FNAME"
      lfs setstripe -c ${STRIPECNT} -S ${STRIPESZ} $FNAME
      echo "Testing done: echo $FNAME"
      lfs getstripe $FNAME

	    #srun $ENVVARS -n ${RANKS} -N ${COBALT_PARTSIZE} --cpu_bind=verbose,cores -c ${num_logical_cores} $EXE ${ARG} > $OUTPUT
		  xtnodestat > xtnodestat.start.${OUTPUT}
		  qstat -f > qstat.start.${OUTPUT}
			aprun ${ENVVARS} ${APRUNPARAMS} ${EXE} ${ARG} > ${OUTPUT}
		  qstat -f > qstat.end.${OUTPUT}
 			xtnodestat > xtnodestat.end.${OUTPUT}
    done
   done
  done

	echo 
	echo "* * * * *"
	echo 
 done
done

EXIT_STATUS=$?
echo "Job $COBALT_JOBID completed."
exit $EXIT_STATUS
