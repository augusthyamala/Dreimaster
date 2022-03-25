#include "TH1.h"
#include "TFile.h"
#include <iostream>

void biascorrection( Int_t run )
{

  TFile *f = new TFile( Form("drei-r%i.root",run));

  if( f->IsOpen() )
    cout << "success!" << endl;
  else
    cout << "failed to open file" << endl; 

  TH1I *xBcog = new TH1I();
  TH1I *xBias = new TH1I();

  xBcog = (TH1I*) f->Get("xBmodptch");
  xBias = (TH1I*) xBcog->GetCumulative();
  xBias->SetName("xBias");
  xBias->Scale( 25/xBcog->Integral() );  //Scale( pitch/.... )
  xBias->Draw("hist");
 
  TFile *file = new TFile( Form("bias-r%i.root", run ), "RECREATE" );
  xBias -> Write();

}
