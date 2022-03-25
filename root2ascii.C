#include "TDirectory.h"
#include "TObject.h"
#include "TH1.h"
#include <iostream>
#include <fstream>

using namespace std;

void root2ascii( string h, Int_t run )
{
  TObject * obj = gDirectory->Get( h.c_str() );

  if( obj == NULL ) {
    cout << h << " does not exist\n";
    return;
  }

  if( !obj->InheritsFrom( "TH1" ) ) {
    cout << h << " is a not a 1-D historgram\n";
    return;
  }

  TH1 * hs = (TH1*)obj;
  hs -> Draw();
 
    Int_t n = hs ->GetNbinsX();
 
    ofstream file;
    file.open( Form("drei-%i.txt",run) );

    file << "run " << run << endl;
    file << "plot : " << h << endl;
    file << "no. of bins : " << n << endl;
    file << " " << endl;
    file << "   x   " << "   value   " << endl;

    Int_t i;

    for( i = 0; i <= n + 1 ; i++ )
	file << hs->GetBinLowEdge(i) + hs->GetBinWidth(i)/2 << "        " << hs->GetBinContent(i) << endl;
      
    cout << "drei-"<< run << ".txt file has been created" << endl;
 
}
