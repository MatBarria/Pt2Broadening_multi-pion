#include "Acc.h"

int main(int argc, char* argv[]) {

  if(argc != 2) {
    std::cout << "Insert (just) the target name as a parameter" << std::endl;
    return 0;
  }

  TStopwatch t;

  std::cout << "Start" << std::endl;

  // For the Target name use (C,Fe,Pb) for the solids targets and (DC,DFe,DPb) for the liquid target
  std::string target = argv[1];
  // Creating a array of chars instead of a string to use Form method
  int n = target.length();
  char targetArr[n + 1];
  strcpy(targetArr, target.c_str());

  // Select the solid target

  int m;
  TString fileDataName;
  if(targetArr[0] == 'D') {
    m = 2;
    char solidTarget[n];
    for(int i = 0; i < n; i++){
      solidTarget[i] = targetArr[i+1];
    }
    //fileDataName = Form("/eos/user/m/mbarrial/Data/VecSum_%s.root", solidTarget);
    fileDataName = Form("/home/matias/proyecto/Piones/Data/VecSum/VecSum_%s.root", solidTarget);
  } else{
    m = n;
    //fileDataName = Form("/eos/user/m/mbarrial/Data/VecSum_%s.root", targetArr);
    fileDataName = Form("/home/matias/proyecto/Piones/Data/VecSum/VecSum_%s.root", targetArr);
  }
  char simulTarget[m];
  if(targetArr[0] == 'D') {
    simulTarget[0] = 'D';
    simulTarget[1] = '\0';
  } else{
    for(int i = 0; i < n; i++){
      simulTarget[i] = targetArr[i];
    }
  }

  // Open the input and output files
  TFile* fileData   = new TFile(fileDataName,"READ");
  TFile *fileOutput= new TFile(Form("/home/matias/proyecto/Piones/Data/VecSum/test/corr_data_Phi_%s.root", targetArr), "RECREATE");
  gROOT->cd();

  // Create some variables to use inside the for loops
  TString tupleDataName;
  TCut Q2Cut, NuCut, ZhCut, Pt2Cut, YCCut, VCData, cutsData;

  // Create all the necessary histograms
  TH1F *histData        = new TH1F("Data",            "", N_Phi, -180, 180);

  // Store the sum of the weights A.K.A the erros (in the other histograms if save it by other methods)
  histData->Sumw2();
  // histThrown->Sumw2();
  // histTotDetected->Sumw2();
  // histDetected->Sumw2();

  std::cout << Form("Simul target %s, Target %s", simulTarget, targetArr) << std::endl;
  // Select liquid liquid or solid target
  if(targetArr[0] == 'D') { VCData  = "VC_TM == 1.";}
  else {VCData  = "VC_TM == 2.";}

  for(int gen = 1; gen <= N_PION ; gen++) { // Loops in every number of generated pions
    // Open the tuple of gen number of pions
    // TNtuple* ntuple_data = (TNtuple*) fileData->Get(Form("ntuple_%i_pion",gen));
    // Loops in every bin
    for(int Q2Counter = 0; Q2Counter < N_Q2; Q2Counter++) { // Loops in every Q2 bin
      for(int NuCounter = 0; NuCounter < N_Nu; NuCounter++) { // Loops in every Nu bin
        for(int ZhCounter = 0; ZhCounter < N_Zh; ZhCounter++) { // Loops in every Zh bin

          std::cout << "Bin selected: " << gen << Q2Counter << NuCounter << ZhCounter << std::endl;
          Q2Cut   = Form("Q2>%f&&Q2<%f", Q2_BINS[Q2Counter],   Q2_BINS[Q2Counter+1]);
          NuCut   = Form("Nu>%f&&Nu<%f", Nu_BINS[NuCounter],   Nu_BINS[NuCounter+1]);
          ZhCut   = Form("Zh>%f&&Zh<%f", Zh_BINS[ZhCounter],   Zh_BINS[ZhCounter+1]);
          YCCut   = "TMath::Abs(YC)<1.4";

          cutsData = Q2Cut&&NuCut&&ZhCut&&YCCut&&VCData;

          TNtuple* ntupleData = (TNtuple*) fileData->Get(Form("ntuple_%i_pion", gen));
          ntupleData->Draw(">>listData", cutsData);
          TEventList* evntData = (TEventList*) gDirectory->Get("listData");
          ntupleData->SetEventList(evntData);

          for(int Pt2Counter = 0; Pt2Counter < N_Pt2; Pt2Counter++) { // Loops in every Pt2 bin

            // Generate the cuts depeding on the selected bin
            // Q2_cut   = Form("Q2>%f&&Q2<%f", Q2_BINS[Q2Counter],   Q2_BINS[Q2Counter+1]);
            // Nu_cut   = Form("Nu>%f&&Nu<%f", Nu_BINS[NuCounter],   Nu_BINS[NuCounter+1]);
            // Zh_cut   = Form("Zh>%f&&Zh<%f", Zh_BINS[ZhCounter],   Zh_BINS[ZhCounter+1]);
            // Pt2_cut  = Form("Pt2>%f&&Pt2<%f", Pt2_BINS[Pt2Counter], Pt2_BINS[Pt2Counter+1]);
            // YC_cut   = "TMath::Abs(YC)<1.4";


            //std::cout << "LOOP CUT = " << cuts_data << std::endl;
            Pt2Cut  = Form("Pt2>%f&&Pt2<%f", Pt2_BINS[Pt2Counter], Pt2_BINS[Pt2Counter+1]);
            ntupleData->Project("Data","PhiPQ", Pt2Cut);
            if(EmptyHist(histData) == 1){continue;}


            // Save the histograms in the output file
            fileOutput->cd();

            histData->Write(Form("DataCorr2_%s_%i%i%i%i_%i",             targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));


            gROOT->cd();

            // Set the histograms values to 0
            histData->Reset();


          }
          delete ntupleData;
          delete evntData;
        }
      }
    }
    //delete ntuple_data;
  }
  fileOutput->Close();
  fileData->Close();
  t.Print();

}
