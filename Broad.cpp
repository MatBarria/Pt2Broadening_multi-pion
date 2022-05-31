#include "Broad.h"

int main(){

  TStopwatch t;

  TString inputDirectory  = "~/proyecto/Piones/Data/VecSum/test/";
  TString outputDirectory = "~/proyecto/Piones/Data/VecSum/test/";
  TString plotDirectory   = "~/proyecto/Piones/Plots/VecSum/test/";

  CallPhiIntegration(inputDirectory , outputDirectory);
  CallNuZhIntegration(inputDirectory , outputDirectory);
  CallQ2NuIntegration(inputDirectory , outputDirectory);
  CallQ2ZhIntegration(inputDirectory , outputDirectory);
  CallQ2NuZhIntegration(inputDirectory , outputDirectory);

  PtBroadeningQ2Integrated(inputDirectory, plotDirectory);
  PtBroadeningNuIntegrated(inputDirectory, plotDirectory);
  PtBroadeningZhIntegrated(inputDirectory, plotDirectory);
  PtBroadeningFullIntegrated(inputDirectory, plotDirectory);

  t.Print();

}
