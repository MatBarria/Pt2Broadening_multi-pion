cd ..

mkdir -p bin

g++ -Wall -fPIC -I./include `root-config --cflags` Integration.cpp -o ./bin/Integration  `root-config --glibs` ./include/Broad.h
g++ -Wall -fPIC -I./include `root-config --cflags` AccCorrection.cpp -o ./bin/AccCorrection  `root-config --glibs` ./include/Acc.h
g++ -Wall -fPIC -I./include `root-config --cflags` Broadening.cpp -o ./bin/Broadening  `root-config --glibs` ./include/Broad.h

cd bin

./AccCorrection DC
./AccCorrection DFe
./AccCorrection DPb
./AccCorrection C
./AccCorrection Pb
./AccCorrection Fe

cd ..
cd Data

hadd -f corr_data_Phi.root corr_data_Phi_C.root corr_data_Phi_Fe.root corr_data_Phi_Pb.root corr_data_Phi_DC.root corr_data_Phi_DFe.root corr_data_Phi_DPb.root

cd ..
cd bin
./Integration
./Broadening
