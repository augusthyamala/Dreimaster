
// Daniel Pitzl (DESY) Sep 2017
// read region-of-interest data: event display

// make edroi

// ./edroi B/roi001020.txt

#include <stdlib.h> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <cmath>
#include <unistd.h> // usleep
#include <sys/ioctl.h>

#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2.h>
#include <TLine.h>

using namespace std;

struct pixel {
  int col;   // column 0..154
  int row;   // 0..159
  double ph; // [ADC] pulse height
  double q;  // [ke] charge
};

//------------------------------------------------------------------------------
bool kbhit()
{
  usleep(1000); // [us]
  int byteswaiting;
  ioctl( 0, FIONREAD, &byteswaiting );
  return byteswaiting > 0;
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give file name" << endl;
    return 1;
  }

  // file name = last argument:

  string evFileName( argv[argc-1] );
  cout << "try to open  " << evFileName;

  ifstream evFile( argv[argc-1] );

  if( !evFile ) {
    cout << " : failed " << endl;
    return 2;
  }

  cout << " : succeed " << endl;

  // further arguments:

  bool fifty = 0;
  bool sixteen = 0;

  int eoi = -1; // event of interest

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-f" ) )
      fifty = 1;

    if( !strcmp( argv[i], "-s" ) )
      sixteen = 1;

    if( !strcmp( argv[i], "-e" ) )
      eoi = atoi( argv[++i] );

  } // argc

  // set styles:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleFont( 62, "XYZ" );

  bool quer = 0; // 1 for edge-on
  gStyle->SetTitleOffset( 1.3, "x" );
  if( quer ) {
    gStyle->SetTitleOffset( 0.8, "y" );
    gStyle->SetTitleOffset( 0.9, "z" );
  }
  else {
    gStyle->SetTitleOffset( 1.5, "y" );
    gStyle->SetTitleOffset( 1.7, "z" );
  }

  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetLabelOffset( 0.022, "x" );
  if( quer ) {
    gStyle->SetLabelOffset( 0.012, "y" ); // quer
    gStyle->SetLabelOffset( 0.012, "z" );
  }
  else {
    gStyle->SetLabelOffset( 0.022, "y" );
    gStyle->SetLabelOffset( 0.022, "z" );
  }

  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside
  if( quer )
    gStyle->SetTickLength( -0.01, "yz" ); // tick marks outside

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(10);
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  gStyle->SetStatBorderSize(1); // no 'shadow'
  gStyle->SetStatX(0.82);
  gStyle->SetStatY(0.92);

  //gStyle->SetPalette(55); // sunset
  gStyle->SetPalette(56); // white to brown
  //gStyle->SetPalette(73); // blue
  //gStyle->SetPalette(90); // green to magenta
  //gStyle->SetPalette(109); // sea blue to magenta
  gStyle->SetNumberContours(32); // -20..300

  TApplication app( "app", 0, 0 );
  TCanvas * c1;
  if( quer ) {
    c1 = new TCanvas( "c1", "pixel event display", 1800, 800 );
    c1->SetLeftMargin( 0.06 );
    c1->SetRightMargin( 0.12 );
  }
  else {
    c1 = new TCanvas( "c1", "pixel event display", 900, 800 ); // square
    c1->SetRightMargin( 0.18 );
  }

  c1->SetTopMargin( 0.12 );

  unsigned nbz =  78;
  int nbx = 320;
  if( fifty ) {
    nbz = 155;
    nbx = 160;
  }
  else if( sixteen ) {
    nbz = 52;
    nbx = 480;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Read file by lines:

  string START {"START"};
  string hd;

  while( hd != START ) {
    getline( evFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  string F {"F"}; // filled flag

  string evseed;
  getline( evFile, evseed ); // read one line into string

  bool more = 1;
  string q{"q"};

  while( evFile.good() && ! evFile.eof() && more ) { // events

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;
    cout << "ev " << iev;

    string filled;
    iss >> filled;

    if( iev < eoi ) { // not wanted

      if( filled == F ) {
	string roi;
	getline( evFile, roi ); // read and continue
      }

    }

    else if( filled == F ) {

      string roi;
      getline( evFile, roi );

      istringstream css( roi ); // tokenize string

      int npx = 0;
      vector <pixel> vpx;
      vpx.reserve(35);

      cout << " px";

      while( ! css.eof() ) { // pixels

	int col;
	int row;
	double ph;
	css >> col;
	css >> row;
	css >> ph;

	++npx;

	pixel px { col, row, ph, ph };
	vpx.push_back(px);

      } // roi px

      if( vpx.size() > 1999 ) {
	cout << " " << iev << ":" << vpx.size() << flush;
	vpx.clear();
      }

      TH2D hpxmap( "pxmap",
		   Form( "pixel map %i;col;row;PH [ADC]", iev ),
		   nbz, -0.5, nbz-0.5, nbx, -0.5, nbx-0.5 );
      //hpxmap.SetMinimum(-20);
      //hpxmap.SetMaximum(300);

      int nhit{0}; // counter

      bool draw{0}; // flag

      if( eoi > 0 ) draw = 1; // event of interest

      // columns-wise common mode correction:

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {

	int col4 = vpx[ipx].col;
	int row4 = vpx[ipx].row;
	double ph4 = vpx[ipx].ph;

	int row1 = row4;
	int row7 = row4;
	double ph1 = ph4;
	double ph7 = ph4;

	for( unsigned jpx = 0; jpx < vpx.size(); ++jpx ) {

	  if( jpx == ipx ) continue;
	  if( vpx[jpx].col != col4 ) continue; // want same column

	  int jrow = vpx[jpx].row;

	  if( jrow < row1 ) {
	    row1 = jrow;
	    ph1 = vpx[jpx].ph;
	  }

	  if( jrow > row7 ) {
	    row7 = jrow;
	    ph7 = vpx[jpx].ph;
	  }

	} // jpx

	if( row4 == row1 ) continue; // Randpixel
	if( row4 == row7 ) continue;

	double dph;
	if( row4 - row1 < row7 - row4 )
	  dph = ph4 - ph1;
	else
	  dph = ph4 - ph7;

	//cout << " " << col << " " << row << " " << ph;

	pixel px;

	if( fifty ) {
	  px.col = col4;
	  px.row = row4;
	}

	else if( sixteen ) {

	  px.col = (col4+2)/3; // 150 um

	  if(       (col4+1)%6 == 0 )
	    px.row = 3*row4 + 0;
	  else if(  (col4+1)%6 == 1 )
	    px.row = 3*row4 + 1;
	  else if(  (col4+1)%6 == 2 )
	    px.row = 3*row4 + 1;
	  else if(  (col4+1)%6 == 3 )
	    px.row = 3*row4 + 2;
	  else if(  (col4+1)%6 == 4 )
	    px.row = 3*row4 + 0;
	  else
	    px.row = 3*row4 + 2;

	}

	else {

	  px.col = (col4+1)/2; // 100 um

	  if( col4%2 )
	    px.row = 2*row4 + 0;
	  else
	    px.row = 2*row4 + 1;

	}

	hpxmap.Fill( px.col, px.row, dph ); // ROI, dph can go negative

	if( dph > 12 ) { // threshold

	  ++nhit;

	  //hpxmap.Fill( px.col, px.row, dph ); // hits

	} // dph

      } // roi px

      cout << ", roi size " << npx
	   << ", hits " << nhit
	   << endl;

      //if( nhit > 1 ) draw = 1;
      if( nhit > 11 ) draw = 1; // delta rays
      //if( nhit > 33 ) draw = 1;

      if( draw ) {

	hpxmap.Draw( "colz" );

	c1->Update();

	cout << "enter any key, q to stop" << endl;

	while( !kbhit() ) gSystem->ProcessEvents(); // ROOT

	string any;
	cin >> any;
	if( any == q )
	  more = 0;

      } // show

    } // filled

    else
      cout << "  empty" << endl;

    getline( evFile, evseed ); // read ahead

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << "done " << evFileName
       << endl;

  return 0;

}
