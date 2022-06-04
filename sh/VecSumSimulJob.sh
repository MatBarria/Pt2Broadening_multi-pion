TARNAME="$@"

# set dirs
MAINDIR=/eos/user/${USER:0:1}/${USER}/Pt2Broadening_multi-pion              # dir of the program
JOBDIR=/eos/user/${USER:0:1}/${USER}/Pt2Broadening_multi-pion/bin          # dir to store logs and job scripts
#MAINDIR=/home/matias/proyecto/Pt2Broadening_multi-pion              # dir of the program
#JOBDIR=/home/matias/proyecto/Pt2Broadening_multi-pion/jobs          # dir to store logs and job scripts

mkdir -p ${JOBDIR} # just in case

# setting jobname
jobname="VecSumSim_${TARNAME}"
jobfile="${JOBDIR}/${jobname}.sh"

echo ${jobname}

echo "#!/bin/bash"                                                 > ${jobfile}
echo "#SBATCH -J ${jobname}"                                      >> ${jobfile}
echo "#SBATCH -o ${JOBDIR}/${jobname}.out"                        >> ${jobfile}
echo "#SBATCH -e ${JOBDIR}/${jobname}.err"                        >> ${jobfile}
echo "#SBATCH --time=3:00:00"                                     >> ${jobfile}
echo "#SBATCH --mem=1GB"                                          >> ${jobfile}
echo ""                                                           >> ${jobfile}
echo "source ${HOME}/.bashrc"                                     >> ${jobfile}
echo "cd ${JOBDIR}"                                               >> ${jobfile}
echo "./VecSumSimul ${TARNAME}"                                   >> ${jobfile}
echo "Submitting job: ${jobfile}"
sbatch ${jobfile} # submit job!
