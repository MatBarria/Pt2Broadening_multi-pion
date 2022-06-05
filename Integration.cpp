#include "Broad.h"

int main(){

  TStopwatch t;

  TString inputDirectory  = "/home/matias/proyecto/Pt2Broadening_multi-pion/Data/";
  TString outputDirectory = "/home/matias/proyecto/Pt2Broadening_multi-pion/Data/";

  std::cout << "PhiPQ integration" << std::endl;
  CallPhiIntegration(   inputDirectory , outputDirectory);
  std::cout << "Nu and Zh integration" << std::endl;
  CallNuZhIntegration(  inputDirectory , outputDirectory);
  std::cout << "Q2 and Nu integration" << std::endl;
  CallQ2NuIntegration(  inputDirectory , outputDirectory);
  std::cout << "Q2 and Zh integration" << std::endl;
  CallQ2ZhIntegration(  inputDirectory , outputDirectory);
  std::cout << "Q2, Nu and Zh integration" << std::endl;
  CallQ2NuZhIntegration(inputDirectory , outputDirectory);

  // PtBroadeningQ2Integrated(  inputDirectory, plotDirectory);
  // PtBroadeningNuIntegrated(  inputDirectory, plotDirectory);
  // PtBroadeningZhIntegrated(  inputDirectory, plotDirectory);
  // PtBroadeningFullIntegrated(inputDirectory, plotDirectory);

  t.Print();

}
