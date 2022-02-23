
// Daniel Pitzl (DESY) Sep 2017
// read region-of-interest data

// make rdroi

// rdroi B/roi001020.txt

#include <stdlib.h> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <cmath>
#include <sys/time.h> // gettimeofday, timeval

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>

using namespace std;

struct pixel {
  int col;   // column 0..154
  int row;   // 0..159
  double ph; // [ADC] pulse height
  double q;  // [ke] charge
};

struct cluster {
  vector <pixel> vpix;
  int sum; // [ADC]
  double q; // [ke] charge
  double col, row; // [pixel] center-of-gravity
};

//------------------------------------------------------------------------------
vector<cluster> getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with local coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> v;
  if( pb.size() == 0 ) return v;

  int * gone = new int[pb.size()] {0};

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut)
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list of clusters:

    c.sum = 0;
    c.q = 0;
    c.col = 0;
    c.row = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      double ph = p->ph;
      c.sum += ph;
      double q = p->q;
      c.q += q;
      c.col += (*p).col*q;
      c.row += (*p).row*q;
    }

    c.col /= c.q;
    c.row /= c.q;

    v.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  // nothing left,  return clusters

  delete [] gone;

  return v;

} // getClus

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

  // B/roi001020.txt

  int run = stoi( evFileName.substr( 5, 6 ) );
  cout << evFileName << "  " << evFileName.substr( 5, 6 ) << "  " << run << endl;

  // further arguments:

  int Nev = 9000*1000;
  bool fifty = 0;

  for( int i = 1; i < argc; i++ ) {

    if( !strcmp( argv[i], "-n" ) )
      Nev = atoi( argv[++i] );

    if( !strcmp( argv[i], "-f" ) )
      fifty = 1;

  } // argc

  // gain:

  double p0[155][160]; // Fermi
  double p1[155][160];
  double p2[155][160];
  double p3[155][160];

  string gain{ "B/r146-scancal-tb21-1208.dat" };

  double ke = 0.0354; // Landau peak at 11 ke

  ifstream gainFile( gain );

  bool haveGain = 0;

  if( gainFile ) {

    haveGain = 1;

    while( ! gainFile.eof() ) {

      int icol;
      int irow;
      gainFile >> icol;
      gainFile >> irow;
      gainFile >> p0[icol][irow];
      gainFile >> p1[icol][irow];
      gainFile >> p2[icol][irow];
      gainFile >> p3[icol][irow];

    } // while

  } // gain
  else
    cout << "missing " << gain << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile * histoFile = new TFile( Form( "rdroi-%i.root", run ), "RECREATE" );

  // book histos:

  TH1I hsz( "sz", "event size;event size [ROI pixels];events", 500, 0, 500 );
  TH1I hszl( "szl", "event size;log_{10}(event size [ROI pixels]);events", 80, 1, 5 );

  TH1I hph( "ph", "PH;ADC-PED [ADC];pixels", 1000, -100, 900 );
  TH1I hdph( "dph", "dPH;#DeltaPH [ADC];pixel", 1000, -100, 900 );

  TH1I hg200( "g200", "200mV/PH;pixel gain [mv/ADC];pixels in clusters", 100, 0, 2 );

  int nbx =  78;
  int nby = 320;
  if( fifty ) {
    nbx = 155;
    nby = 160;
  }
  TH2I * hpxmap = new TH2I( "pxmap", "pixel map, dph > cut;col;row;pixels above cut",
			    nbx, -0.5, nbx-0.5, nby, -0.5, nby-0.5 );
  TH1I hnpx( "npx", "PH pixels per event;PH pixels;events", 80, 0.5, 80.5 );

  TH1D hncl( "ncl", "cluster per event;cluster;events", 20, 0.5,  20.5 );

  TH2I * hclmap = new TH2I( "clmap", "cluster map;col;row;clusters above noise",
			    nbx, -0.5, nbx-0.5, nby, -0.5, nby-0.5 );
  TH1I hclcol( "clcol", "cluster column;column;clusters", nbx, -0.5, nbx-0.5 );
  TH1I hclrow( "clrow", "cluster row;row;clusters", nby, -0.5, nby-0.5 );
  TH1I hclsz( "clsz", "cluster size;cluster size [pixels];clusters", 80, 0.5, 80.5 );

  TH1I hclph( "clph", "cluster PH;cluster ph [ADC];clusters", 200, 0, 1000 );
  TH1I hclq( "clq", "cluster charge;cluster charge [ke];cluster", 100, 0, 20 );
  TH1I hclq1( "clq1", "cluster charge;cluster charge [ke];cluster", 200, 0, 100 );

  TProfile clqvscol( "clqvscol", "Q vs column;column;<cluster charge> [ke]", nbx, -0.5, nbx-0.5 );

  TH1I hpxph( "pxph", "pixel PH;pixel ph [ADC];pixels in clusters", 200, 0, 400 );
  TH1I hpxq( "pxq", "pixel charge;pixel charge [ke];pixels in clusters", 200, 0, 20 );

  TH1I hcolph( "colph", "column PH;column ph [ADC];columns", 100, 0, 500 );

  TH1I hncol( "ncol", "cluster cols;cluster size [cols];clusters", 90, 0.5, 90.5 );
  TH1I hnrow( "nrow", "cluster rows;cluster size [rows];clusters", 90, 0.5, 90.5 );

  const double log10 = log(10);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Read file by lines:

  string START {"START"};
  string hd;

  while( hd != START ) {
    getline( evFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  timeval tv;
  gettimeofday( &tv, NULL );
  long s0 = tv.tv_sec; // seconds since 1.1.1970
  long u0 = tv.tv_usec; // microseconds

  string F {"F"}; // filled flag
  string BLANK{" "};

  string evseed;
  getline( evFile, evseed ); // read one line into string
  int nev = 0;

  while( evFile.good() && ! evFile.eof() && nev < Nev ) {

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;

    if( iev%1000 == 0 ) cout << " " << iev << flush;

    vector <pixel> pb; // for clustering

    string filled;
    iss >> filled;

    if( filled == F ) {

      int npx = 0;
      vector <pixel> vpx;
      vpx.reserve(35);

      string roi;
      getline( evFile, roi );
      //cout << roi << " (" << roi.size() << ")" << endl;

      /* istringstream css( roi ); // tokenize string
	 while( ! css.eof() ) { // slower
	 int col;
	 int row;
	 double ph;
	 css >> col;
	 css >> row;
	 css >> ph;
	 pixel px { col, row, ph, ph };
	 vpx.push_back(px);
	 ++npx;
	 } // roi stream
      */

      // faster:

      size_t start = 0;
      size_t gap = 0;
      while( gap < roi.size()-1 ) { // data have trailing blank

        gap = roi.find( BLANK, start );
	string s1( roi.substr( start, gap - start ) );
	int col = stoi(s1);
	start = gap + BLANK.size();

        gap = roi.find( BLANK, start );
	string s2( roi.substr( start, gap - start ) );
	int row = stoi(s2);
	start = gap + BLANK.size();

        gap = roi.find( BLANK, start );
	string s3( roi.substr( start, gap - start ) );
	double ph = stod(s3);
	start = gap + BLANK.size();

	pixel px { col, row, ph, ph };
	vpx.push_back(px); // comment out = no clustering
	++npx;

      }

      hsz.Fill( vpx.size() + 0.5 );
      if( vpx.size() )
	hszl.Fill( log(vpx.size())/log10 );

      if( vpx.size() > 1999 ) // huge event
	vpx.clear(); // save time

      // column-wise common mode correction: (pedestal subtraction)

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {

	int col4 = vpx[ipx].col; // pixel in roi
	int row4 = vpx[ipx].row;
	double ph4 = vpx[ipx].ph; // [ADC]

	int row1 = row4; // 1 will become bottom row in roi
	int row7 = row4; // 7 will vecome top row in roi
	double ph1 = ph4;
	double ph7 = ph4;

	// find bottom and top row in this column:

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

	if( row4 == row1 ) continue; // bottom of roi
	if( row4 == row7 ) continue; // top of roi

	double dph; // difference pulse height
	if( row4 - row1 < row7 - row4 )
	  dph = ph4 - ph1; // 4 is closer to 1
	else
	  dph = ph4 - ph7; // 4 is closer to 7

	hph.Fill( ph4 );
	hdph.Fill( dph );

	if( dph > 12 ) { // offline threshold [ADC]

	  pixel px;

	  if( fifty ) {
	    px.col = col4;
	    px.row = row4;
	  }
	  else{
	    px.col = (col4+1)/2; // 100 um
	    if( col4%2 )
	      px.row = 2*row4 + 0;
	    else
	      px.row = 2*row4 + 1;
	  }

	  px.ph = dph;

	  if( haveGain ) { // from r4scal.C

	    double U = ( dph - p3[col4][row4] ) / p2[col4][row4];

	    if( U >= 1 )
	      U = 0.9999999; // avoid overflow

	    double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	    // subtract Vcal offset:

	    double U0 = -p3[col4][row4] / p2[col4][row4]; // dph = 0
	    double v0 = p0[col4][row4] - p1[col4][row4] * log( (1-U0)/U0 ); // inverse Fermi

	    double q = ke*(vcal-v0);

	    px.q = q;

	  }
	  else
	    px.q = dph;

	  pb.push_back(px);

	  hpxmap->Fill( px.col, px.row );

	} // dph

      } // ipx

      //cout << "ev " << iev << " roi " << npx << ", hits " << pb.size() << endl;

    } // filled

    //else cout << "  empty" << endl;

    // clustering:

    hnpx.Fill( pb.size() ); // pixels above threshold

    vector <cluster> vcl = getClus(pb);

    hncl.Fill( vcl.size() );
    //if( vcl.size() ) cout << "  clusters " << vcl.size();

    for( unsigned icl = 0; icl < vcl.size(); ++ icl ) {

      hclmap->Fill( vcl[icl].col, vcl[icl].row );
      hclcol.Fill( vcl[icl].col );
      hclrow.Fill( vcl[icl].row );
      hclsz.Fill( vcl[icl].vpix.size() );

      hclph.Fill( vcl[icl].sum );
      hclq.Fill( vcl[icl].q ); // peak
      hclq1.Fill( vcl[icl].q ); // tail

      clqvscol.Fill( vcl[icl].col, vcl[icl].q );

      // pixels in clusters:

      int colmin = 999;
      int colmax = 0;
      int rowmin = 999;
      int rowmax = 0;
      double colph[155] { 0 };

      for ( unsigned ipx = 0; ipx < vcl[icl].vpix.size(); ++ipx ) {

	int col = vcl[icl].vpix[ipx].col;
	if( col < colmin ) colmin = col;
	if( col > colmax ) colmax = col;

	int row = vcl[icl].vpix[ipx].row;
	if( row < rowmin ) rowmin = row;
	if( row > rowmax ) rowmax = row;

	colph[col] += vcl[icl].vpix[ipx].ph;

	hpxph.Fill( vcl[icl].vpix[ipx].ph );
	hpxq.Fill( vcl[icl].vpix[ipx].q );

      } // px

      int ncol = colmax - colmin + 1;
      hncol.Fill( ncol );

      int nrow = rowmax - rowmin + 1;
      hnrow.Fill( nrow );

    } // cl

    ++nev;

    getline( evFile, evseed ); // read ahead

  } // while events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  gettimeofday( &tv, NULL );
  long s9 = tv.tv_sec; // seconds since 1.1.1970
  long u9 = tv.tv_usec; // microseconds

  cout << endl << "done " << evFileName
       << endl << nev << " events"
       << " in " << s9 - s0 + ( u9 - u0 ) * 1e-6 << " s"
       << endl;

  cout << endl << histoFile->GetName()
       << endl << endl;

  histoFile->Write();
  histoFile->Close();

  return 0;
}
