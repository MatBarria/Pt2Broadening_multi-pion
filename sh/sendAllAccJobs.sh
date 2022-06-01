MAINDIR=/eos/user/${USER:0:1}/${USER}/Pt2Broadening_multi-pion
BINDIR=${MAINDIR}/bin

mkdir -p ${BINDIR}

g++ -Wall -fPIC -I${MAINDIR}/include `root-config --cflags` ${MAINDIR}/Acc.cpp -o ${BINDIR}/Acc `root-config --glibs` ${MAINDIR}/include/Acc.h

bash Acc_job.sh DC
bash Acc_job.sh DFe
bash Acc_job.sh DPb

bash Acc_job.sh C
bash Acc_job.sh Fe
bash Acc_job.sh Pb

echo "All jobs sumitted"
