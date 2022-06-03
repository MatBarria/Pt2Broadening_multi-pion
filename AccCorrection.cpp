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
  TFile* fileGenDec = new TFile(Form("/home/matias/proyecto/Piones/Data/VecSum/Acc/Gen_Dec_%s.root", simulTarget), "READ");
  //TFile *fileOutput = new TFile(Form("/eos/user/m/mbarrial/Data/Acc/corr_data_Phi_%s.root", targetArr), "RECREATE");
  TFile *fileOutput= new TFile(Form("/home/matias/proyecto/Piones/Data/VecSum/Acc/corr_data_Phi_%s1.root", targetArr), "RECREATE");
  gROOT->cd();

  // Create some variables to use inside the for loops
  TString tupleDataName;
  TCut Q2Cut, NuCut, ZhCut, Pt2Cut, YCCut, VCData, cutsData;

  // Create all the necessary histograms
  TH1F *histDetected    = new TH1F("histDetected",    "", N_Phi, -180, 180);
  TH1F *histTotDetected = new TH1F("histTotDetected", "", N_Phi, -180, 180);
  TH1F *histThrown      = new TH1F("histThrown",      "", N_Phi, -180, 180);
  TH1F *histData        = new TH1F("Data",            "", N_Phi, -180, 180);
  TH1F* histFalPos      = new TH1F("FalPosFactor",    "", N_Phi, -180, 180);
  TH1F* histAccFactors  = new TH1F("AcceFactor",      "", N_Phi, -180, 180);
  TH1F* histDataCorr    = new TH1F("DataCorr",        "", N_Phi, -180, 180);
  TH1F* histDataCorr2   = new TH1F("DataCorr2",       "", N_Phi, -180, 180);

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
            // Generate histograms of the all dectected pion, all generated pion, and the pions that was correct dectected
            histDetected    = (TH1F*) fileGenDec->Get(Form("histDetected_%s_%i%i%i%i_%i",    simulTarget, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));
            histThrown      = (TH1F*) fileGenDec->Get(Form("histThrown_%s_%i%i%i%i_%i",      simulTarget, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));
            histTotDetected = (TH1F*) fileGenDec->Get(Form("histTotDetected_%s_%i%i%i%i_%i", simulTarget, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));

            // Calculate the Acceptance factor
            histAccFactors->Divide(histDetected, histThrown, 1, 1, "B");
            // Calculate a factor that reprensent how many of the detected as N pion events are truly N pions events
            histFalPos->Divide(histDetected, histTotDetected, 1, 1, "B");
            // Create the data histogram depending on of the cuts
            //ntuple_data->Project("Data","PhiPQ", cuts_data);
            AccHist1(histAccFactors);
            AccHist1(histFalPos);
            // Apply the correction factors
            histDataCorr->Divide(histData, histAccFactors, 1, 1);
            histDataCorr2->Multiply(histDataCorr, histFalPos, 1, 1);


            // Save the histograms in the output file
            fileOutput->cd();

            //histData->Write(Form("Data_%s_%i%i%i%i_%i",             targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));
            histDataCorr2->Write(Form("DataCorr2_%s_%i%i%i%i_%i",   targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));
            histDataCorr->Write(Form("DataCorr_%s_%i%i%i%i_%i",     targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));
            //histFalPos->Write(Form("FalPosFactor_%s_%i%i%i%i_%i",   targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));
            //histAccFactors->Write(Form("AcceFactor_%s_%i%i%i%i_%i", targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));

            gROOT->cd();

            // Set the histograms values to 0
            histData->Reset();
            histDataCorr2->Reset();
            histDataCorr->Reset();
            histFalPos->Reset();
            histAccFactors->Reset();
            histThrown->Reset();
            histDetected->Reset();
            histTotDetected->Reset();

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
  fileGenDec->Close();
  t.Print();

}
