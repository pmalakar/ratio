#!/bin/bash -x
#COBALT --disable_preboot
 
#export L1P_POLICY=std
#export BG_THREADLAYOUT=1   # 1 - default next core first; 2 - my core first

#Free bootable blocks
boot-block --reboot
 
NODES=$1

for PROG in testcase
do
for iter in 1 # 2 
do
for THRD in 1 #2 4 8 16
do
for ppn in 1 #16 #8 #4 #2 1
do
for MSG in 16 #128 #256 512 1024 2048 # 4096 8192
do 
 for coalesced in 0 #1 
 do
 for streams in 0 #1 #2 #0 #1 2
 do
 for blocking in 0 #1
 do
 rm -f dummy*
 for type in 1 #0 1 2 
 do
	if [ $type -gt 0 ] && [ $coalesced -eq 1 ] 
	then
			continue;
	fi
	if [ $type -gt 0 ] && [ $streams -gt 0 ] 
	then
			continue;
	fi
  if [ $type -eq 0 ] && [ $streams -gt 0 ] && [ $coalesced -eq 0 ] 
	then
			continue;
	fi
	RANKS=`echo "$NODES*$ppn"|bc`
	OUTPUT=${PROG}_N${NODES}_R${ppn}_${MSG}_${coalesced}_${blocking}_${type}_${streams}
	rm -f ${OUTPUT}.cobaltlog ${OUTPUT}.output ${OUTPUT}.error
	echo 
	echo "* * * * *"
	echo 
	echo "Starting $OUTPUT with numthreads=$THRD ppn=$ppn args=${MSG} ${coalesced} ${blocking} ${type} ${streams}"

	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}

	continue;

	if [ $type -gt 0 ] || [ $coalesced -eq 0 ] 
	then
			continue;
	fi
	if [ $streams -eq 0 ] 
	then
			continue;
	fi

 	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 PAMID_RZV_LOCAL=4M --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2_rzvlocal_4M
  	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV_LOCAL=4M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_local_4M

 	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 PAMID_RZV_LOCAL=4M MUSPI_RECFIFOSIZE=2097152 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2_rec_2M_rzvlocal_4M




	continue;



 	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 PAMID_RZV_LOCAL=8M --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2_rzvlocal_8M
 	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 PAMID_RZV_LOCAL=16M --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2_rzvlocal_16M

	if [ $MSG -gt 128 ] && [ $streams -gt 0 ] 
	then
    rm -f dummy*
	  runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=8 MUSPI_NUMRECFIFOS=8 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_8
    rm -f dummy*
	  runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_RECFIFOSIZE=8388608" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_recfifo_8M
  	#rm -f dummy*
	#	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_INJFIFOSIZE=16777216" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injfifo_16M
  	#rm -f dummy*
	#	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_RECFIFOSIZE=16777216" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_recfifo_16M
	fi

  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2
  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=4 MUSPI_NUMRECFIFOS=4 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_4
  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs PAMID_RZV_LOCAL=2M MUSPI_RECFIFOSIZE=2097152 --envs PAMID_STATISTICS=1 --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_local_2M_recfifo_2M
  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs PAMID_RZV_LOCAL=4M MUSPI_RECFIFOSIZE=2097152 --envs PAMID_STATISTICS=1 --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_local_4M_recfifo_2M
 	rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 PAMID_RZV_LOCAL=4M MUSPI_RECFIFOSIZE=4194304 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2_rec_4M_rzvlocal_4M

	if [ $MSG -gt 64 ] 
	then
			continue;
	fi

  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_RECFIFOSIZE=2097152" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_recfifo_2M
  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_INJFIFOSIZE=2097152" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injfifo_2M

	#streams > 1
	if [ $streams -gt 1 ] 
	then
  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_RECFIFOSIZE=4194304" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_recfifo_4M
	fi

	if [ $MSG -gt 4 ] 
	then
			continue;
	fi

  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_RECFIFOSIZE=2097152 MUSPI_INJFIFOSIZE=2097152 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injrec_2M
  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV=4M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_remote_4M
  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_INJFIFOSIZE=4194304" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injfifo_4M
  rm -f dummy*
	runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV_LOCAL=4M" --envs "PAMID_RZV=4M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_4M

  #rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_THREAD_MULTIPLE=1" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_ptm
  #rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "MUSPI_INJFIFOSIZE=8388608" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injfifo_8M

	#if [ $MSG -gt 32 ] 
	#then
	#		continue;
	#fi
	
  #rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_RECFIFOSIZE=4194304 MUSPI_INJFIFOSIZE=4194304 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_injrec_4M
  #rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs PAMID_RZV_LOCAL=4M PAMID_RZV=4M MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_4M_numfifos2
  #rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs "PAMID_RZV=8M" --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_rzv_remote_8M
 	#rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 PAMID_RZV=4M PAMID_RZV_LOCAL=4M MUSPI_INJFIFOSIZE=4194304 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2_inj_4M_rzv_4M
 	#rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=2 MUSPI_NUMRECFIFOS=2 PAMID_RZV=8M PAMID_RZV_LOCAL=8M MUSPI_INJFIFOSIZE=4194304 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_2_inj_4M_rzv_8M

 	#rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=4 MUSPI_NUMRECFIFOS=4 PAMID_RZV=4M PAMID_RZV_LOCAL=4M MUSPI_INJFIFOSIZE=4194304 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_4_inj_4M_rzv_4M
 	#rm -f dummy*
	#runjob --np $RANKS -p $ppn --block $COBALT_PARTNAME --verbose=INFO --envs "OMP_MAX_NUM_THREADS=${THRD}" --envs MUSPI_NUMINJFIFOS=4 MUSPI_NUMRECFIFOS=4 PAMID_RZV=8M PAMID_RZV_LOCAL=8M MUSPI_INJFIFOSIZE=4194304 --envs "PAMID_STATISTICS=1" --envs "PAMID_VERBOSE=1" : ${PROG} ${MSG} ${coalesced} ${blocking} ${type} ${streams} > ${OUTPUT}_numfifos_4_inj_4M_rzv_8M


	echo 
	echo "* * * * *"
	echo
done 
done 
done 
done 
done 
done 
done 
done 
done

exit

