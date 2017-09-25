#!/bin/bash
#SBATCH -p regular
#SBATCH -A m2885
#SBATCH -C knl,quad,flat
#SBATCH -t 06:00:00
#SBATCH -J my_job
#SBATCH -o my_job.o%j
#SBATCH -S 2

#module swap craype-haswell craype-mic-knl
#module load darshan
module load hwloc

hwloc-info
hwloc-gather-topology

lstopo 

EXE=./testcase

#io tracing
export DXT_ENABLE_IO_TRACE=4

#mpich_abort_on_rw_wrror 
export ATP_ENABLED=1

export MPICH_MPIIO_HINTS_DISPLAY=1
export MPICH_MPIIO_AGGREGATOR_PLACEMENT_DISPLAY=1 
export MPICH_MPIIO_STATS=1
export MPICH_MPIIO_TIMERS=1
export MPICH_MPIIO_ABORT_ON_RW_ERROR=enable

export MPICH_RANK_REORDER_DISPLAY=true
export MPICH_RANK_REORDER_DISPLAY=1
export MPICH_CPUMASK_DISPLAY=true
export MPICH_CPUMASK_DISPLAY=1
export KMP_AFFINITY=verbose
export MPICH_VERSION_DISPLAY=1
export MPICH_NEMESIS_ASYNC_PROGRESS=1
export MPICH_ENV_DISPLAY=1

ENVVARS="" #"--export=MPICH_ENV_DISPLAY=1 --export=MPICH_VERSION_DISPLAY=1 --export=OMP_NUM_THREADS=2"

echo "Running Job $SLURM_JOB_ID on $SLURM_JOB_NUM_NODES $SLURM_JOB_NODELIST"
scontrol show hostnames $SLURM_JOB_NODELIST

maxranks=64
totalranks=$((${maxranks}*${SLURM_JOB_NUM_NODES}))

srun -n ${SLURM_JOB_NUM_NODES} -N ${SLURM_JOB_NUM_NODES} ./location.x 
srun -n ${SLURM_JOB_NUM_NODES} -N ${SLURM_JOB_NUM_NODES} ./status.knl 

startiter=1
startiter=$1
enditer=$(($startiter))

squeue | grep -v " 1 nid" > squeue.start.$SLURM_JOB_ID

for iter in `seq $startiter $enditer` 
do
for ppn in 32 #64
do
  RANKS=$(($SLURM_JOB_NUM_NODES * $ppn)) 
  num_logical_cores=$(($maxranks / $ppn))
	echo 
  for STRIPECNT in 16 
  do 
   for STRIPESZ in 16M 
   do
    for size in 512
    do 

      OUTPUT=output_${SLURM_JOB_NUM_NODES}_${RANKS}_R${ppn}_${STRIPECNT}_${STRIPESZ}_${size}_${iter}_${SLURM_JOB_ID}

	    echo "Starting $OUTPUT with $ARG"
      echo `date`
      #mkdir pat_${OUTPUT}
      #export PAT_RT_EXPFILE_DIR=pat_${OUTPUT}
      FNAME=TestFile-${RANKS}
      rm $FNAME 2>/dev/null
      lfs setstripe -c ${STRIPECNT} -S ${STRIPESZ} $FNAME
      lfs getstripe $FNAME
	    srun --profile=All $ENVVARS -n ${RANKS} -N ${SLURM_JOB_NUM_NODES} --cpu_bind=verbose,cores -c ${num_logical_cores} $EXE ${ARG} > $OUTPUT
	#srun $ENVVARS -n ${RANKS} -N ${SLURM_JOB_NUM_NODES} --cpu_bind=verbose,rank -c 2 $EXE ${ARG} > $OUTPUT
	#srun $ENVVARS -n ${RANKS} -N ${SLURM_JOB_NUM_NODES} -c 2 $EXE ${ARG} > $OUTPUT

      echo `date`
    done
   done
  done

	echo 
	echo "* * * * *"
	echo 
done
done

squeue | grep -v " 1 nid" > squeue.end.$SLURM_JOB_ID

EXIT_STATUS=$?
echo "Job $SLURM_JOB_ID completed."
exit $EXIT_STATUS

