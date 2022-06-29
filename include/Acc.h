#include <iostream>
#include <string>
#include "Binning.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TVector2.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TH1.h"
#include "TCut.h"
#include "TEventList.h"

// If the histogram if empty return 1 if not return 0
int EmptyHist(TH1F* h) {

  int empty = 0;
  for(int i = 1 ; i <= h->GetNbinsX() ; i++) {
    if(h->GetBinContent(i) == 0){ empty++; }
  }
  if(empty == h->GetNbinsX()) { return 1; }
  else { return 0; }

}

// Check if the number of the correctly detected events is bigger than one
void AccCondition(TH1F* hist) {

  for(Int_t bin = 1; bin <= hist->GetNbinsX(); bin++) {
    if(hist->GetBinContent(bin) == 1) {
      hist->SetBinContent(bin, 0);
      hist->SetBinError(bin, 0);
    }
  }

}

// If Acceptance Factor > 1 set it to 0
void AccHist1(TH1F* hist) {

  for(Int_t bin = 1; bin <= hist->GetNbinsX(); bin++) {
    if(hist->GetBinContent(bin) > 1) {
      hist->SetBinContent(bin, 0);
      hist->SetBinError(bin, 0);
    }
  }

}
