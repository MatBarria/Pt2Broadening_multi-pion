#include "Broad.h"

int main(){

  TStopwatch t;

  TString inputDirectory  = "~/proyecto/Piones/Data/VecSum/test/";
  TString outputDirectory = "~/proyecto/Piones/Data/VecSum/test/";
  TString plotDirectory   = "~/proyecto/Piones/Plots/VecSum/test/";

  PtBroadeningQ2Integrated(  inputDirectory, plotDirectory);
  PtBroadeningNuIntegrated(  inputDirectory, plotDirectory);
  PtBroadeningZhIntegrated(  inputDirectory, plotDirectory);
  PtBroadeningFullIntegrated(inputDirectory, plotDirectory);

  t.Print();

}
