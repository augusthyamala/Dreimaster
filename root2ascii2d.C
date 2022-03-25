#include "TDirectory.h"
#include "TObject.h"
#include "TH2.h"
#include <iostream>
#include <fstream>

using namespace std;

void root2ascii2d( string h, Int_t run )
{
  TObject * obj = gDirectory->Get( h.c_str() );

  if( obj == NULL ) {
    cout << h << " does not exist\n";
    return;
  }

  if( !obj->InheritsFrom( "TH2" ) ) {
    cout << h << " is a not a 2-D historgram\n";
    return;
  }

  TH2 * hs = (TH2*)obj;
  hs -> Draw();

  Int_t l = hs -> GetNbinsX();
  Int_t m = hs -> GetNbinsY();

  ofstream file;
  file.open( Form("drei-2D-%i.txt",run) );

  file << "run " << run << endl;
  file << "plot : " << h << endl;
  file << "no. of bins along x-axis : " << l << endl;
  file << "no. of bins along y-axis : " << m << endl;
  file << " " << endl;
  file << "   clsz   " << "   x   " << "   value   " << endl;

  Int_t i, j; 

  for( i = 0 ; i <= l + 1 ; i++ ){
    for(j = 0 ; j <= m + 1 ; j++ ){ 
      Double_t ax = hs->GetXaxis()->GetBinLowEdge(i);
      Double_t ay = hs->GetYaxis()->GetBinLowEdge(j) + ( hs->GetYaxis()->GetBinWidth(j) )/2;
      Double_t b = hs->GetBinContent(i,j);
      file << ax << "        " << ay << "        " << b << endl;
    }
  }

  cout << l << " " << m << endl;
}
