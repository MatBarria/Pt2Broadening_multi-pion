#include "Broad.h"

int main(){

  TStopwatch t;

  TString inputDirectory  = "/home/matias/proyecto/Pt2Broadening_multi-pion/Data/";
  TString outputDirectory = "/home/matias/proyecto/Pt2Broadening_multi-pion/Data/";
  TString plotDirectory   = "/home/matias/proyecto/Pt2Broadening_multi-pion/Plots/";

  PtBroadeningQ2Integrated(  inputDirectory, plotDirectory);
  PtBroadeningNuIntegrated(  inputDirectory, plotDirectory);
  PtBroadeningZhIntegrated(  inputDirectory, plotDirectory);
  PtBroadeningFullIntegrated(inputDirectory, plotDirectory);

  t.Print();

}
