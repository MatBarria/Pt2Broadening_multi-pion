TARNAME="$@"

# set dirs
MAINDIR=/eos/user/${USER:0:1}/${USER}/Pt2Broadening_multi-pion              # dir of the program
JOBDIR=/eos/user/${USER:0:1}/${USER}/Pt2Broadening_multi-pion/jobs          # dir to store logs and job scripts
#MAINDIR=/home/matias/proyecto/Pt2Broadening_multi-pion              # dir of the program
#JOBDIR=/home/matias/proyecto/Pt2Broadening_multi-pion/jobs          # dir to store logs and job scripts

mkdir -p ${JOBDIR} # just in case

# setting jobname
jobname="ACC_MulPion_${TARNAME}"
jobfile="${JOBDIR}/${jobname}.sh"

echo ${jobname}

echo "#!/bin/bash"                                                 > ${jobfile}
echo "#SBATCH -J ${jobname}"                                      >> ${jobfile}
echo "#SBATCH -o ${JOBDIR}/${jobname}.out"                        >> ${jobfile}
echo "#SBATCH -e ${JOBDIR}/${jobname}.err"                        >> ${jobfile}
echo "#SBATCH --time=4:00:00"                                     >> ${jobfile} # 4hrs or 15min for test
echo "#SBATCH --mem=2GB"                                          >> ${jobfile}
echo ""                                                           >> ${jobfile}
echo "source ${HOME}/.bashrc"                                     >> ${jobfile}
echo "cp ${MAINDIR}/bin/Acc ${JOBDIR}"                            >> ${jobfile} # retrieve executable
echo "cd ${JOBDIR}"                                               >> ${jobfile}
echo "./Acc ${TARNAME}"                                           >> ${jobfile}
echo "rm Acc"                                                     >> ${jobfile} # remove binary from temp dir
echo "Submitting job: ${jobfile}"
sbatch ${jobfile} # submit job!
