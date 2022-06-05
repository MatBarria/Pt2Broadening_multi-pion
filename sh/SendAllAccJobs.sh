MAINDIR=/eos/user/${USER:0:1}/${USER}/Pt2Broadening_multi-pion
BINDIR=${MAINDIR}/bin

mkdir -p ${BINDIR}

g++ -Wall -fPIC -I${MAINDIR}/include `root-config --cflags` ${MAINDIR}/AccCorrection.cpp -o ${BINDIR}/AccCorrection  `root-config --glibs` ${MAINDIR}/include/Acc.h

bash AccJob.sh DC
bash AccJob.sh DFe
bash AccJob.sh DPb

bash AccJob.sh C
bash AccJob.sh Fe
bash AccJob.sh Pb

echo "All jobs sumitted"
