void pionSeparator(string strMat) {

  int n = strMat.length();
  char material[n + 1];
  strcpy(material, strMat.c_str());

  TString inputName;
  TString tupleName;
  TString outputName;

  inputName = Form("/work/emolina/%s_data_Npion.root", material);
  tupleName = "ntuple_data";
  outputName = Form("/eos/user/m/mbarrial/Data/pSeparated_%s.root", material);
  std::cout << "data " << std::endl;

  TFile* file = new TFile(inputName ,"READ");
  TNtuple* tuple = (TNtuple*)file->Get(tupleName);

  const char* VarList = "Q2:Nu:PhiPQ:Zh:Pt2:Xf:Xb:YC:PID:Mass2:VC_TM";
  Float_t *vars         = new Float_t[12];
  tuple->SetBranchAddress("Q2",&vars[0]);
  tuple->SetBranchAddress("Nu",&vars[1]);
  tuple->SetBranchAddress("PhiPQ",&vars[2]);
  tuple->SetBranchAddress("Zh",&vars[3]);
  tuple->SetBranchAddress("Pt2",&vars[4]);
  tuple->SetBranchAddress("Xf",&vars[5]);
  tuple->SetBranchAddress("Xb",&vars[6]);
  tuple->SetBranchAddress("YC",&vars[7]);
  tuple->SetBranchAddress("PID",&vars[8]);
  tuple->SetBranchAddress("Mass2",&vars[9]);;
  tuple->SetBranchAddress("VC_TM",&vars[10]);
  tuple->SetBranchAddress("NmbPion",&nmbPion);

  gROOT->cd();

  TNtuple *ntuple_1pion = new TNtuple("ntuple_1_pion","1 Pion pluses",VarList);
  TNtuple *ntuple_2pion = new TNtuple("ntuple_2_pion","2 Pion pluses",VarList);
  TNtuple *ntuple_3pion = new TNtuple("ntuple_3_pion","3 Pion pluses",VarList);

  int tmpCounter = 0;
  float tmpEvnt;
  float tmpZh[4];
  for(int i = 0; i < tuple->GetEntries() - 2 ; i++) {
    tuple->GetEntry(i);
    tmpEvnt = vars[0];
    tuple->GetEntry(i + 1);
    while(tmpEvnt == vars[0]) {
      tmpCounter++;
      tmpEvnt = vars[0];
      tuple->GetEntry(i + 1 + tmpCounter);
    }
    if(tmpCounter == 0) {
      tuple->GetEntry(i);
      vars[32] = 1;
      ntuple_1pion->Fill(vars);
    }else {
      vars[5] = 0;m
      TVector2* vec = new TVector2(0,0);
      for(int k = 0; k <= tmpCounter; k++) {
        // Calculate de tranvers momentum vector
        TVector2 *tmpVec = new TVector2(tmpPt[k]*TMath::Cos((tmpPhi[k] + 180)*TMath::DegToRad()), tmpPt[k]*TMath::Sin((tmpPhi[k] + 180)*TMath::DegToRad()));
        // Sum the vector and save the sum of Zh
        vars[5] += tmpZh[k];
        *vec += *tmpVec;
        //vecTemp->Print();
        delete tmpVec;
      }
      // Save the Pt2 of the sum vector
      vars[6] = std::pow(vec->Mod(),2);
      // Save the PhiPQ of the sum vector
      vars[7] = vec->Phi()*TMath::RadToDeg()-180;
      delete vec;

    }
    if(tmpCounter == 1) {
      for(int j = 0; j <= tmpCounter; j++) {
        tuple->GetEntry(i + j);
        vars[32] = 1 + j;
        ntuple_2pion->Fill(vars);
      }
    }
    if(tmpCounter == 2){
      for(int j = 0; j <= tmpCounter; j++) {
        tuple->GetEntry(i + j);
        vars[32] = 1 + j;
        ntuple_3pion->Fill(vars);
      }
    }
    // Jump to the next event
    i += tmpCounter;
    tmpCounter = 0;
  }

  TFile* fOutput = new TFile(outputName,"RECREATE");

  fOutput->cd();
  ntuple_1pion->Write();
  ntuple_2pion->Write();
  ntuple_3pion->Write();
  gROOT->cd();
  fOutput->Close();
  cout << "Done." << endl;
  delete ntuple_1pion;
  delete ntuple_2pion;
  delete ntuple_3pion;
  file->Close();

}
