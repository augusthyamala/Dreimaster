
// Daniel Pitzl (DESY) Sep 2017
// Jan 2018: version for openMP (faster)
// Jun 2021: sensor-level cross talk
// 3 x R4S, 17x150, 25x100 or 50x50 on rot90 PCB

// drei -s 4685   # 17x150
// drei    4702   # 25x100

#include <cstdlib> // atoi
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <list>
#include <cmath>
#include <time.h> // clock_gettime
#include <sched.h> // getcpu
#include <sys/resource.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

using namespace std;

struct evInfo {
  int iev;
  uint64_t evtime;
  bool skip;
  string filled;
};

struct pixel {
  int col;
  int row;
  double q;
};

struct cluster {
  vector <pixel> vpix;
  int size; // [px]
  double q; // [ke]
  double col, row; // [px]
  bool iso;
};

const int A{0};
const int B{1};
const int C{2};
string PN[]{"A","B","C"};

double p0[3][155][160]; // Fermi
double p1[3][155][160];
double p2[3][155][160];
double p3[3][155][160];
double ke[3];

TProfile phvsprev[3];
TProfile dphvsprev[3];
TH1I hph[3];
TH1I hdph[3];
TH1I hnpx[3];
TProfile npxvsev[3];
TH1I hnht[3];
TH2I * hpxmap[3];
TH1I hncl[3];
TH2I * hclmap[3];
TH1I hclsz[3];
TH1I hclq[3];

list < evInfo > infoA;
list < evInfo > infoB;
list < evInfo > infoC;

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

    // start a new cluster:
    
    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ) { // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( ( dr >= -fCluCut ) && ( dr <= fCluCut ) &&
		( dc >= -fCluCut ) && ( dc <= fCluCut ) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // p, important!
            }
          } // p loop over vpix
        } // not gone
      } // i loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list of clusters:

    c.q = 0;
    c.size = 0;
    c.col = 0;
    c.row = 0;
    c.iso = 1;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      double q = p->q;
      c.q += q;
      c.col += p->col * q;
      c.row += p->row * q;
    }

    c.size = c.vpix.size();
    c.col /= c.q;
    c.row /= c.q;

    v.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  // nothing left,  return clusters

  delete [] gone;
  return v;
}

//------------------------------------------------------------------------------
list < vector < cluster > > oneplane( int plane, string runnum, unsigned Nev, int geo )
{
  int run = stoi( runnum );

  list < vector < cluster > > evlist;

  string Xfile;
  if( plane == A ) {
    Xfile = "A/roi000" + runnum + ".txt";
    if( run > 999 )
      Xfile = "A/roi00" + runnum + ".txt";
  }
  if( plane == B ) {
    Xfile = "B/roi000" + runnum + ".txt";
    if( run > 999 )
      Xfile = "B/roi00" + runnum + ".txt";
  }
  if( plane == C ) {
    Xfile = "C/roi000" + runnum + ".txt";
    if( run > 999 )
      Xfile = "C/roi00" + runnum + ".txt";
  }
  cout << "try to open  " << Xfile;
  ifstream Xstream( Xfile.c_str() );
  if( !Xstream ) {
    cout << " : failed " << endl;
    return evlist;
  }
  cout << " : succeed " << endl;

  // typical: Moyal peak 200 ADC = 11.2 ke => thr 12 = 670 e
  //double pxqcut = 0.4; // [ke] 25x100 4702 2.97
  //double pxqcut = 0.5; // [ke] 25x100 4702 2.96
  double pxqcut = 1.1; // [ke] 25x100 4702 2.96
  //double pxqcut = 0.7; // [ke] 25x100 4702 2.96
  //double pxqcut = 0.9; // [ke] 25x100 4702 2.965
  //double pxqcut = 1.2; // [ke] 25x100 4702 2.973 flatter eta
  //double pxqcut = 1.5; // [ke] 25x100 4702 2.984 flatter eta
  //double pxqcut = 1.8; // [ke] 25x100 4702 3.014 flatter eta
  //double pxqcut = 2.1; // [ke] 25x100 4702 3.061 flatter eta

  double cx{0}; // cross talk
  if( geo == 16 )
    //cx = 0.00; // 4685: 2.22
    //cx = 0.09; // 4685: 2.18
    //cx = 0.12; // 4685: 2.17
    //cx = 0.15; // 4685: 2.15
    //cx = 0.18; // 4685: 2.13
    cx = 0.19; // 4685: 2.13
    //cx = 0.21; // 4685: 2.13
    //cx = 0.24; // 4685: 2.16
  else if( geo == 25 )
    //cx = -0.03; // 4702 3.00
    cx = 0.00; // 4702 2.96
    //cx = 0.03; // 4702 3.04
    //cx = 0.06; // 4702 3.09

  string START {"START"};
  string hd; // header

  while( hd != START ) {
    getline( Xstream, hd ); // read one line into string
    cout << "  " << hd << endl;
  }

  string F {"F"}; // filled flag
  string E {"E"}; // empty flag
  string Add {"A"}; // added flag
  string BLANK{" "};

  bool ldb = 0;
  uint64_t prevtime = 0;

  if( plane == A && run == 2786 ) {
    string evseed;
    getline( Xstream, evseed ); // fix
    getline( Xstream, evseed ); // fix
  }
  if( plane == C && run == 2786 ) {
    string evseed;
    getline( Xstream, evseed ); // fix
  }

  while( Xstream.good() && ! Xstream.eof() &&
	 evlist.size() < Nev ) {

    string evseed;
    getline( Xstream, evseed ); // read one line into string

    istringstream iss( evseed ); // tokenize string
    int iev;
    iss >> iev;

    if( iev%1000 == 0 )
      cout << "  " << PN[plane] << " " << iev << flush;

    string filled;
    iss >> filled;

    int iblk; // event block number: 100, 200, 300...
    iss >> iblk;

    if( plane == C && run == 2786 && iev == 800 && filled == Add )
      getline( Xstream, evseed ); // fix

    uint64_t evtime = prevtime; // for added events

    if( filled == F )
      iss >> evtime; // from run 456

    else if( filled == E )
      iss >> evtime; // from run 456

    vector <pixel> pb; // for clustering
    vector <cluster> vcl;
    evInfo evinf;
    evinf.iev = iev;
    evinf.evtime = evtime;
    evinf.filled = filled;
    evinf.skip = 0;
    prevtime = evtime;

    if( filled == F ) {

      string roi;
      getline( Xstream, roi );

      size_t start = 0;
      size_t gap = 0;
      unsigned ng = 0; // 3 gaps per pix
      string BLANK{" "};
      while( gap < roi.size()-1 ) { // data have trailing blank
	gap = roi.find( BLANK, start );
	start = gap + BLANK.size();
	++ng;
      }
      hnpx[plane].Fill( ng/3 );
      npxvsev[plane].Fill( iev, ng/3 ); // ROI size

      map <int,double> mpx[155]; // mpx[col][row]=dph;

      if( ng/3 < 400 ) {

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
	  hph[plane].Fill( ph );
	  start = gap + BLANK.size();

	  mpx[col][row] = ph;

	}
      } // size
      else {
	evinf.skip = 1;
	//cout << " (" << iev << ": B ROI " << ng/3 << " skip)";
      }

      // R4S column-wise common mode correction:

      double phprev = 0;
      double dphprev = 0;

      map < int, map <int,double> > spx; // spx[col][row]=dph sensor pixel without common mode

      for( int col = 0; col < 155; ++col ) {

	if( mpx[col].size() < 3 ) continue; // sensor edge

	auto px = mpx[col].begin();
	//int row1 = px->first;
	double ph1 = px->second;

	auto px7 = mpx[col].end();
	--px7; // last
	//int row7 = px7->first;
	double ph7 = px7->second;

	++px;
	for( ; px != px7; ++px ) {

	  int row = px->first;
	  double ph = px->second;
	  phvsprev[plane].Fill( phprev, ph ); // correlation?

	  phprev = px->second; // original ph

	  double dph = ph - 0.5*(ph1+ph7); // subtract baseline and common mode

	  hdph[plane].Fill( dph ); // sig 2.7

	  dphvsprev[plane].Fill( dphprev, dph );
	  dphprev = dph;

	  // r4scal.C

	  double U = ( dph - p3[plane][col][row] ) / p2[plane][col][row];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = p0[plane][col][row] - p1[plane][col][row] * log( (1-U)/U );
	  // inverse Fermi

	  // subtract Vcal offset:

	  double U0 = -p3[plane][col][row] / p2[plane][col][row];
	  double v0 = p0[plane][col][row] - p1[plane][col][row] * log( (1-U0)/U0 );
	  // inverse Fermi

	  double pxq = ke[plane] * ( vcal - v0 );

	  // map r4s to sens:

	  int scol{col};
	  int srow{row}; // for 50x50

	  if( geo == 16 ) {

	    scol = (col+2)/3; // 150 um

	    if(       (col+1)%6 == 0 )
	      srow = 3*row + 0;
	    else if(  (col+1)%6 == 1 )
	      srow = 3*row + 1;
	    else if(  (col+1)%6 == 2 )
	      srow = 3*row + 1;
	    else if(  (col+1)%6 == 3 )
	      srow = 3*row + 2;
	    else if(  (col+1)%6 == 4 )
	      srow = 3*row + 0;
	    else
	      srow = 3*row + 2;

	  }
	  else if( geo == 25 ) {

	    scol = (col+1)/2; // 100 um

	    if( col%2 )
	      srow = 2*row + 0;
	    else
	      srow = 2*row + 1;
	  }

	  spx[scol][srow] = pxq; // sensor pixel

	} // px

      } // col

      // sensor cross talk correction:

      for( auto col = spx.begin(); col != spx.end(); ++col ) {

	if( col->second.size() < 2 ) continue; // sensor edge

	auto row1{col->second.begin()};
	double q1{row1->second};

	auto row2{row1};
	++row2;

	while( row2 != col->second.end() ) {

	  double q2{row2->second};
	  row1->second = (1-cx)*q1 + cx*q2;
	  row2->second = (1-cx)*q2 + cx*q1;
	  q1 = q2;
	  ++row1;
	  ++row2;

	} // row

      } // col

      for( auto col = spx.begin(); col != spx.end(); ++col ) {

	if( col->second.size() < 2 ) continue; // sensor edge

	for( auto row = col->second.begin();row != col->second.end(); ++row ) {

	  if( row->second < pxqcut ) continue;
	 
	  hpxmap[plane]->Fill( col->first, row->first );

	  pixel px;
	  px.col = col->first;
	  px.row = row->first;
	  px.q = row->second;

	  pb.push_back(px); // sensor pixels for clustering

	} // row

      } // col

    } // filled

    // clustering:

    hnht[plane].Fill( pb.size() );

    if( pb.size() > 50 ) pb.clear(); // speed

    vcl = getClus(pb);

    hncl[plane].Fill( vcl.size() );
    if( vcl.size() )
      if( ldb ) cout << "  clusters " << vcl.size() << endl;

    for( unsigned icl = 0; icl < vcl.size(); ++icl ) {

      hclmap[plane]->Fill( vcl[icl].col, vcl[icl].row );

      hclq[plane].Fill( vcl[icl].q );

      hclsz[plane].Fill( vcl[icl].size );

      // cluster isolation:

      for( unsigned jcl = icl+1; jcl < vcl.size(); ++jcl ) {

	bool done = 0;

	for( unsigned ipx = 0; ipx < vcl[icl].vpix.size(); ++ipx ) {

	  for( unsigned jpx = 0; jpx < vcl[jcl].vpix.size(); ++jpx )

	    if( fabs( vcl[icl].vpix[ipx].col - vcl[jcl].vpix[jpx].col ) < 3 &&
		fabs( vcl[icl].vpix[ipx].row - vcl[jcl].vpix[jpx].row ) < 3 ) {

	      if( vcl[icl].q < vcl[jcl].q ) // Thu 22.3.2018
		vcl[icl].iso = 0; // flag smaller cluster
	      else
		vcl[jcl].iso = 0;

	      done = 1;
	      break; // jpx

	    }

	  if( done ) break; // ipx

	} // ipx

      } // jcl

    } // icl

    evlist.push_back(vcl);
    if( plane == A )
      infoA.push_back(evinf);
    else if( plane == B )
      infoB.push_back(evinf);
    else if( plane == C )
      infoC.push_back(evinf);

  } // events

  cout << endl;
  cout << "done " << Xfile << ", read " << evlist.size() << " events" << endl;

  return evlist;

} // oneplane

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{

  //std::ios_base::sync_with_stdio( false ); // faster read ?

  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give run number" << endl;
    return 1;
  }

  // file name = last argument:

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  // further arguments:

  int Nev = 20*1000*1000;
  double pp = 5.6; // [GeV]
  int geo = 25;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      Nev = atoi( argv[++i] );

    if( !strcmp( argv[i], "-p" ) )
      pp = atof( argv[++i] );

    if( !strcmp( argv[i], "-f" ) )
      geo = 50; // rot90

    if( !strcmp( argv[i], "-s" ) )
      geo = 16;

  } // argc

  cout << "geo " << geo << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // gains:

  string gainA{ "A/r113-scancal-tb21-0923.dat"};
  if( run >=  423 ) gainA = "A/r113-scancal-tb21-0923.dat";
  if( run >=  430 ) gainA = "A/r112-scancal-tb21-0925.dat";
  if( run >=  439 ) gainA = "A/r113-scancal-tb21-0923.dat";
  if( run >=  866 ) gainA = "A/r109-scancal-tb21-1112.dat";
  if( run >=  998 ) gainA = "A/r146-scancal-tb21-1208.dat";
  if( run >= 1010 ) gainA = "A/r163-scancal-tb21-1209.dat";
  if( run >= 2731 ) gainA = "B/scm146-scancal1-drei-2018-06-16-ia92-hold33.dat";
  if( run >= 2769 ) gainA = "B/scm148-scancal1-drei-2018-06-18-ia99-hold33.dat";
  if( run >= 2822 ) gainA = "B/scm148-scancal1-drei-2018-06-21-ia82-hold33.dat";
  if( run >= 4626 ) gainA = "A/scancal-c2506-ia45-hold50.dat";
  if( run >= 4645 ) gainA = "A/scancal_c2506_ia44_hold50_vb200_20210513.dat";
  if( run >= 4699 ) gainA = "A/scancal_c2508_ia45_20210515.dat";
  if( run >= 4751 ) gainA = "A/scancal-c1113-ia45-hold50-20210520.dat";
  if( run >= 4761 ) gainA = "A/scancal_c2514_ia45_20210521.dat";

  ke[A] = 0.039; // Landau peak at 11 ke
  if( run >= 423 ) ke[A] = 0.0396; // Landau peak at 11 ke
  if( run >= 432 ) ke[A] = 0.039; // r112
  if( run >= 866 ) ke[A] = 0.038; // r109
  if( run >= 998 ) ke[A] = 0.0354; // r146
  if( run >= 1010 ) ke[A] = 0.0385; // r163 -V0
  if( run >= 2731 ) ke[A] = 0.0430; // 11 ke 146
  if( run >= 2769 ) ke[A] = 0.0450; // 11 ke 148
  if( run >= 2822 ) ke[A] = 0.0475; // 11.3 ke 148
  if( run >= 2833 ) ke[A] = 0.0472; // 11 ke 148
  if( run >= 4639 ) ke[A] = 0.0381; // 22 ke in 285 mu
  if( run >= 4699 ) ke[A] = 0.0435; // 22 ke in 285 mu
  if( run >= 4751 ) ke[A] = 0.0457; // 22 ke in 285 mu
  if( run >= 4761 ) ke[A] = 0.0430; // 22 ke in 285 mu

  ifstream gainFileA( gainA );

  if( ! gainFileA ) {
    cout << "gain file for A not found" << endl;
    return 1;
  }

  while( ! gainFileA.eof() ) {

    int icol;
    int irow;
    gainFileA >> icol;
    gainFileA >> irow;
    gainFileA >> p0[A][icol][irow];
    gainFileA >> p1[A][icol][irow];
    gainFileA >> p2[A][icol][irow];
    gainFileA >> p3[A][icol][irow];

  } // while

  // B:

  string gainB{ "B/r108-scancal-tb21-0921.dat" };
  if( run >=  423 ) gainB = "B/r108-scancal-tb21-0923-hold25.dat";
  if( run >=  430 ) gainB = "A/r117-scancal-tb21-1005.dat";
  if( run >=  432 ) gainB = "B/r110-scancal-tb21-0925-hold25.dat"; // 12.6 ke
  if( run >=  444 ) gainB = "B/r110-scancal-tb21-0928-hold24.dat";
  if( run >=  866 ) gainB = "B/r148-scancal-tb21-1112.dat";
  if( run >=  998 ) gainB = "B/r150-scancal-tb21-1208.dat";
  if( run >= 1010 ) gainB = "B/r146-scancal-tb21-1208.dat";
  if( run >= 1024 ) gainB = "B/r152-scancal-tb21-1210.dat";
  if( run >= 1037 ) gainB = "B/r160-scancal-tb21-2017-12-11.dat";
  if( run >= 2731 ) gainB = "B/scm148-scancal1-drei-2018-06-16-ia99-hold33.dat";
  if( run >= 2769 ) gainB = "B/scm146-scancal1-drei-2018-06-16-ia90-hold33.dat";
  if( run >= 2775 ) gainB = "B/scm150-scancal1-drei-2018-06-18-ia96-hold33.dat";
  if( run >= 2783 ) gainB = "B/scm120i-scancal1-drei-cold-2018-06-18-ia40-vb800-hold33.dat";
  if( run >= 2822 ) gainB = "B/scm163-scancal1-2018-06-21-ia83-hold33.dat";
  if( run >= 4626 ) gainB = "B/scancal-c1111-ia45-hld50.dat";
  if( run >= 4645 ) gainB = "B/scancal_c1111_ia44_hold50_vb200_20210513.dat";
  if( run >= 4699 ) gainB = "B/scancal_c1113_ia45_20210515.dat";
  if( run >= 4751 ) gainB = "B/scancal-c2509-ia45-hold50-20210520.dat";
  if( run >= 4761 ) gainB = "B/scancal-c2509-ia45-hold50-20210520.dat";

  ke[B] = 0.0276; // Landau peak at 11 ke
  if( run >= 423 ) ke[B] = 0.026;
  if( run >= 432 ) ke[B] = 0.0326; // r110
  if( run >= 443 ) ke[B] = 0.036; // cmspixel-daq
  if( run >= 866 ) ke[B] = 0.029; // c148
  if( run >= 998 ) ke[B] = 0.0264; // c150 clcut 2
  if( run >= 1010 ) ke[B] = 0.0354; // r146 -V0
  if( run >= 1024 ) ke[B] = 0.0235; // r152 thicker deep diff at 11 ke
  if( run >= 1037 ) ke[B] = 0.029; // r160
  if( run >= 2731 ) ke[B] = 0.0450; // 11 ke 148
  if( run >= 2769 ) ke[B] = 0.0430; // 11 ke 146
  if( run >= 2775 ) ke[B] = 0.0390; // 11 ke 150
  if( run >= 2783 ) ke[B] = 0.0367; // default for irrad
  if( run >= 2822 ) ke[B] = 0.0486; // 11.3 ke 163
  if( run >= 2833 ) ke[B] = 0.0493; // 11 ke 163
  if( run >= 4639 ) ke[B] = 0.0490; // 22 ke in 285 mu
  if( run >= 4699 ) ke[B] = 0.0394; // 22 ke in 285 mu
  if( run >= 4751 ) ke[B] = 0.0449; // 22 ke in 285 mu
  if( run >= 4761 ) ke[B] = 0.0435; // 22 ke in 285 mu

  ifstream gainFileB( gainB );

  if( ! gainFileB ) {
    cout << "gain file for B not found" << endl;
    return 1;
  }

  while( ! gainFileB.eof() ) {

    int icol;
    int irow;
    gainFileB >> icol;
    gainFileB >> irow;
    gainFileB >> p0[B][icol][irow];
    gainFileB >> p1[B][icol][irow];
    gainFileB >> p2[B][icol][irow];
    gainFileB >> p3[B][icol][irow];

  } // while

  // C:

  string gainC{ "C/r110-scancal-tb21-0921.dat" };
  if( run >= 423 ) gainC = "C/r110-scancal-tb21-0923-hold25.dat";
  if( run >= 432 ) gainC = "C/r114-scancal-tb21-0925-hold25.dat"; // 15.2
  if( run >= 866 ) gainC = "C/r110-scancal-tb21-1112.dat";
  if( run >= 998 ) gainC = "C/r148-scancal-tb21-1208.dat";
  if( run >= 1024 ) gainC = "C/r159-scancal-tb21-1210.dat";
  if( run >= 2731 ) gainC = "B/scm163-scancal1-drei-2018-06-16-ia91-hold33.dat";
  if( run >= 2783 ) gainC = "B/scm163-scancal1-drei-2018-06-18-ia89-hold33.dat";
  if( run >= 2822 ) gainC = "B/scm150-scancal1-drei-2018-06-21-ia97-hold33.dat";
  if( run >= 4626 ) gainC = "C/scancal_c2511_ia44_hold50_vb200_20210513.dat";
  if( run >= 4699 ) gainC = "C/scancal_c2513_ia45_20210515.dat";
  if( run >= 4751 ) gainC = "C/scancal_c2513_ia45_vb149_20210520.dat";
  if( run >= 4761 ) gainC = "C/scancal_c1114_ia45_20210521.dat";

  ke[C] = 0.0366; // Landau peak at 11 ke
  if( run >=  432 ) ke[C] = 0.028; // Landau peak at 11 ke
  if( run >=  866 ) ke[C] = 0.034; // c110
  if( run >=  998 ) ke[C] = 0.0423; // c148 -V0
  if( run >= 1024 ) ke[C] = 0.0381; // c159
  if( run >= 2731 ) ke[C] = 0.0489; // 11 ke 163
  if( run >= 2783 ) ke[C] = 0.0445; // 11 ke 163
  if( run >= 2822 ) ke[C] = 0.0434; // 11.3 ke 150
  if( run >= 2833 ) ke[C] = 0.0439; // 11 ke 150
  if( run >= 4639 ) ke[C] = 0.0414; // 22 ke in 285 mu
  if( run >= 4699 ) ke[C] = 0.0469; // 22 ke in 285 mu
  if( run >= 4751 ) ke[C] = 0.0480; // 22 ke in 285 mu
  if( run >= 4761 ) ke[C] = 0.0474; // 22 ke in 285 mu

  ifstream gainFileC( gainC );

  if( ! gainFileC ) {
    cout << "gain file for C not found" << endl;
    return 1;
  }

  while( ! gainFileC.eof() ) {

    int icol;
    int irow;
    gainFileC >> icol;
    gainFileC >> irow;
    gainFileC >> p0[C][icol][irow];
    gainFileC >> p1[C][icol][irow];
    gainFileC >> p2[C][icol][irow];
    gainFileC >> p3[C][icol][irow];

  } // while

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int aligniteration = 0;

  double alignxA = 0.0; // [mm] same sign as dx
  double alignyA = 0.0; // [mm] same sign as dy
  double alignfA = 0.0; // [rad] same sign dxvsy
  double aligntA = 0.0; // [rad] turn dxvsxAB

  double alignxC = 0.0; // [mm] same sign as dx
  double alignyC = 0.0; // [mm] same sign as dy
  double alignfC = 0.0; // [rad] same sign dxvsy
  double aligntC = 0.0; // [rad] turn dxvsxCB

  string alignFileName = "align_" + runnum + ".dat";

  ifstream alignFile( alignFileName );

  cout << endl;

  if( alignFile.bad() || ! alignFile.is_open() ) {
    cout << "no " << alignFileName << ", will bootstrap" << endl;
  }
  else {

    cout << "read alignment from " << alignFileName << endl;

    string HASH( "#" );
    string ITER( "iteration" );
    string ALXA( "alignxA" );
    string ALYA( "alignyA" );
    string ALFA( "alignfA" );
    string ALTA( "aligntA" );
    string ALXC( "alignxC" );
    string ALYC( "alignyC" );
    string ALFC( "alignfC" );
    string ALTC( "aligntC" );

    while( ! alignFile.eof() ) {

      string line;
      getline( alignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == HASH ) // comments start with #
	continue;

      if( tag == ITER )
	tokenizer >> aligniteration;

      double val;
      tokenizer >> val;
      if(      tag == ALXA )
	alignxA = val;
      else if( tag == ALYA )
	alignyA = val;
      else if( tag == ALFA )
	alignfA = val;
      else if( tag == ALTA )
	aligntA = val;
      else if( tag == ALXC )
	alignxC = val;
      else if( tag == ALYC )
	alignyC = val;
      else if( tag == ALFC )
	alignfC = val;
      else if( tag == ALTC )
	aligntC = val;

      // anything else on the line and in the file gets ignored

    } // while getline

    alignFile.close();

  } // alignFile

  double cfA = cos(alignfA);
  double sfA = sin(alignfA);
  double cfC = cos(alignfC);
  double sfC = sin(alignfC);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile * histoFile = new TFile( Form( "drei-r%i.root", run ), "RECREATE" );

   // book histos:

  int nbc =  80;
  int nbr = 320; // rows for 25x100
  if( abs(geo) == 50 ) {
    nbc = 160;
    nbr = 160; // 50x50
  }
  else if( geo == 16 ) {
    nbc = 52;
    nbr = 480; // 17x150
  }

  for( unsigned ipl = 0; ipl < 3; ++ipl ) {

    phvsprev[ipl] = TProfile( Form( "phvsprev%s", PN[ipl].c_str() ),
			      Form( "%s Tsunami;previous PH [ADC];%s <PH> [ADC]",
				    PN[ipl].c_str(), PN[ipl].c_str() ),
			      80, 0, 800, -999, 1999 );
    dphvsprev[ipl] = TProfile( Form( "dphvsprev%s", PN[ipl].c_str() ),
			      Form( "%s Tsunami;previous #DeltaPH [ADC];%s <#DeltaPH> [ADC]",
				    PN[ipl].c_str(), PN[ipl].c_str() ),
			       80, 0, 800, -999, 1999 );

    hph[ipl] = TH1I( Form( "ph%s", PN[ipl].c_str() ),
		     Form("%s PH;ADC-PED [ADC];%s pixels",
			  PN[ipl].c_str(), PN[ipl].c_str() ),
		     500, -100, 900 );
    hdph[ipl] = TH1I( Form( "dph%s", PN[ipl].c_str() ),
		      Form( "%s #DeltaPH;#DeltaPH [ADC];%s pixels",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      500, -100, 900 );
    hnpx[ipl] = TH1I( Form( "npx%s", PN[ipl].c_str() ),
		      Form( "%s ROI pixels per event;ROI pixels;%s events",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      200, 0, 1000 );
    npxvsev[ipl] = TProfile( Form( "npxvsev%s", PN[ipl].c_str() ),
			     Form( "%s pixels vs time;time [events];<%s ROI size> [pixels]",
				   PN[ipl].c_str(), PN[ipl].c_str() ),
			     7000, 0, 7000*1000, -1, 99999 );
    hnht[ipl] = TH1I( Form( "nht%s", PN[ipl].c_str() ),
		      Form( "%s pixel hits per event;pixel hits;%s events",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      50, 0.5, 50.5 );
    hpxmap[ipl] = new TH2I( Form( "pxmap%s", PN[ipl].c_str() ),
			    Form( "%s pixel map, PH > cut;col;row;%s PH pixels",
				  PN[ipl].c_str(), PN[ipl].c_str() ),
			    155, -0.5, 154.5, 160, -0.5, 159.5 );

    hncl[ipl] = TH1I( Form( "ncl%s", PN[ipl].c_str() ),
		      Form( "%s cluster per event;clusters;%s events",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      21, -0.5, 20.5 );
    hclmap[ipl] = new TH2I( Form( "clmap%s", PN[ipl].c_str() ),
			    Form( "%s cluster map;col;row;%s clusters",
				  PN[ipl].c_str(), PN[ipl].c_str() ),
			    nbc, 0, nbc, nbr, 0, nbr );
    hclsz[ipl] = TH1I( Form( "clsz%s", PN[ipl].c_str() ),
		       Form( "%s cluster size;cluster size [pixels];%s clusters",
			     PN[ipl].c_str(), PN[ipl].c_str() ),
		       20, 0.5, 20.5 );
    hclq[ipl] = TH1I( Form( "clq%s", PN[ipl].c_str() ),
		      Form( "%s cluster charge;cluster charge [ke];%s clusters",
			    PN[ipl].c_str(), PN[ipl].c_str() ),
		      100, 0, 50 );

  } // ipl

  TH1I hdt( "dt", "time between events;log_{10}(#Deltat [s]);events", 100, -4, 1 );
  TH1I hddtAB( "ddtAB", "dtA - dtB;A-B #delta#Deltat [clocks];events", 100, -1000, 1000 );
  TH1I hddtCB( "ddtCB", "dtC - dtB;C-B #delta#Deltat [clocks];events", 100, -1000, 1000 );
  TH1I hddtCA( "ddtCA", "dtC - dtA;C-A #delta#Deltat [clocks];events", 100, -1000, 1000 );
  TProfile ddtvsdtAB( "ddtvsdtAB",
		      "A-B time lag vs intervall;log_{10}(#Deltat [s]);<#delta#Deltat> [clocks]",
		      100, -4, 1,-1E99, 1E99 );
  TProfile ddtvsevAB( "ddtvsevAB", "dtA - dtB vs time;time [events];<#delta#Deltat AB>",
		      7000, 0, 7000*1000, -1E99, 1E99 );
  TProfile ddtvsevCB( "ddtvsevCB", "dtC - dtB vs time;time [events];<#delta#Deltat CB>",
		      7000, 0, 7000*1000, -1E99, 1E99 );
  TProfile ddtvsevCA( "ddtvsevCA", "dtC - dtA vs time;time [events];<#delta#Deltat CA>",
		      7000, 0, 7000*1000, -1E99, 1E99 );

  TH1I skipvsevA( "skipvsevA", "A skipped events vs time;time [events];A skipped events",
		      7000, 0, 7000*1000 );
  TH1I skipvsevB( "skipvsevB", "B skipped events vs time;time [events];B skipped events",
		      7000, 0, 7000*1000 );
  TH1I skipvsevC( "skipvsevC", "C skipped events vs time;time [events];C skipped events",
		      7000, 0, 7000*1000 );

  // correlations:

  TH1I hxA( "xA", "x A;x [mm];clusters A", 100, -5, 5 );
  TH1I hyA( "yA", "y A;y [mm];clusters A", 100, -5, 5 );
  TH1I hxAi( "xAi", "x A isolated;x [mm];isolated clusters A", 100, -5, 5 );
  TH1I hyAi( "yAi", "y A isolated;y [mm];isolated clusters A", 100, -5, 5 );
  TH1I hclqAi( "clqAi", "A isolated cluster charge;cluster charge [ke];A isolatewd clusters",
	       100, 0, 50 );
  
 TH2I * hxxAB = new
    TH2I( "xxAB", "x B vs A;x_{A} [mm];x_{B} [mm];clusters", nbr, -4, 4, nbr, -4, 4 );
  TH2I * hyyAB = new
    TH2I( "yyAB", "y B vs A;y_{A} [mm];y_{B} [mm];clusters", nbc, -4, 4, nbc, -4, 4 );

  double f = 0.5;
  //double f = 0.1; // aligned

  TH1I hdxAB( "dxAB", "Bx-Ax;x-x [mm];cluster pairs", 800, -2, 2 );
  TH1I hdyAB( "dyAB", "By-Ay;y-y [mm];cluster pairs", 400, -2, 2 );
  TProfile dxvsxAB( "dxvsxAB", "dx vs x A-B;x [mm];<dx> [mm]", nbr, -4, 4, -f, f );
  TProfile dxvsyAB( "dxvsyAB", "dx vs y A-B;y [mm];<dx> [mm]",  nbc, -4, 4, -f, f );

  TH2I * hdxvsev = new
    TH2I( "dxvsev", "Bx-Ax vs events;events;#Deltax [px];clusters",
	  100, 0, 10000, 100, -f, f );

  TH1I hdxABC( "dxABC", "AB vs C in x;xC - xAB [mm];AB tracks at C clusters", 800, -2, 2 );
  TH1I hdyABC( "dyABC", "AB vs C in y;yC - yAB [mm];AB tracks at C clusters", 400, -2, 2 );
  TProfile dxvsxABC( "dxvsxABC", "dxABC vs xC;xC [mm];<dxABC> [mm]", nbr, -4, 4, -1, 1 );
  TProfile dxvsyABC( "dxvsyABC", "dxABC vs yC;yC [mm];<dxABC> [mm]",  nbc, -4, 4, -1, 1 );

  TH1I hdminC( "dminC", "dmin AB vs C;mind dABC [mm];AB tracks at C clusters", 200, 0, 1 );

  TProfile effCvsx( "effCvsx", "C efficiency vs x;xAB at C [mm];<eff(C)>",  100, -5, 5, -1, 2 );
  TProfile effCvsy( "effCvsy", "C efficiency vs y;yAB at C [mm];<eff(C)>",  100, -5, 5, -1, 2 );
  TProfile2D * effCvsxy = new
    TProfile2D( "effCvsxy", "C efficiency vs x-y;xAB at C [mm];yAB at C [mm];<eff(C)>",
		100, -5, 5, 100, -5, 5, -1, 2 );
  TProfile effCvsev200( "effCvsev200", "C efficiency vs time;time [events];<eff(C)>",
			1000, 0, 200*1000, -1, 2 );
  TProfile effCvsev( "effCvsev", "C efficiency vs time;time [events];<eff(C)>",
		     7000, 0, 7000*1000, -1, 2 );
  TProfile effCvsdx( "effCvsdx", "C efficiency vs dx;dxAB [mm];<eff(C)>",
		     40, -0.2, 0.2, -1, 2 );
  TProfile effCvsdy( "effCvsdy", "C efficiency vs dy;dyAB [mm];<eff(C)>",
		     40, -0.2, 0.2, -1, 2 );
  TProfile effCvsncl( "effCvsncl", "C efficiency vs multiplicity;clusters A * clusters B;<eff(C)>",
		     40, 0.5, 40.5, -1, 2 );

  TProfile nmvsevAB( "nmvsevAB", "AB matches vs time;time [events];<AB matches>",
		     7000, 0, 7000*1000, -1, 99 );

  TH1I hxB( "xB", "x B;x [mm];clusters B", 100, -5, 5 );
  TH1I hyB( "yB", "y B;y [mm];clusters B", 100, -5, 5 );
  TH1I hxBi( "xBi", "x B isolated;x [mm];isolated clusters B", 100, -5, 5 );
  TH1I hyBi( "yBi", "y B isolated;y [mm];isolated clusters B", 100, -5, 5 );
  TH1I hclqBi( "clqBi", "B isolated cluster charge;cluster charge [ke];B isolatewd clusters",
	       100, 0, 50 );

  TH2I hxcvsxb("xcvsxb","xc vs xb; xB; xC", 100,0,1,100,0,1);

  TH2I * hxxCB = new TH2I( "xxCB", "C vs B;row B;row C;clusters", nbr, -4, 4, nbr, -4, 4 );
  TH2I * hyyCB = new TH2I( "yyCB", "C vs B;col B;col C;clusters",  nbc, -4, 4,  nbc, -4, 4 );

  TH1I hdxCB( "dxCB", "Cx-Bx;x-x [mm];cluster pairs", 800, -2, 2 );
  TH1I hdyCB( "dyCB", "Cy-By;y-y [mm];cluster pairs", 400, -2, 2 );
  TProfile dxvsxCB( "dxvsxCB", "dx vs x C-B;x [mm];<dx> [mm]", nbr, -4, 4, -f, f );
  TProfile dxvsyCB( "dxvsyCB", "dx vs y C-B;y [mm];<dx> [mm]",  nbc, -4, 4, -f, f );

  TH1I hdxCBA( "dxCBA", "CB vs A in x;xA - xCB [mm];CB tracks at A clusters", 800, -2, 2 );
  TH1I hdyCBA( "dyCBA", "CB vs A in y;yA - yCB [mm];CB tracks at A clusters", 400, -2, 2 );
  TProfile dxvsxCBA( "dxvsxCBA", "dxCBA vs xA;xA [mm];<dxCBA> [mm]", nbr, -4, 4, -1, 1 );
  TProfile dxvsyCBA( "dxvsyCBA", "dxCBA vs yA;yA [mm];<dxCBA> [mm]",  nbc, -4, 4, -1, 1 );

  TH1I hdminA( "dminA", "dmin CB vs A;mind dCBA [mm];CB tracks at A clusters", 200, 0, 1 );

  TProfile effAvsx( "effAvsx", "A efficiency vs x;xCB at A [mm];<eff(A)>",  100, -5, 5, -1, 2 );
  TProfile effAvsy( "effAvsy", "A efficiency vs y;yCB at A [mm];<eff(A)>",  100, -5, 5, -1, 2 );
  TProfile2D * effAvsxy = new
    TProfile2D( "effAvsxy", "A efficiency vs x-y;xCB at A [mm];yCB at A [mm];<eff(A)>",
		100, -5, 5, 100, -5, 5, -1, 2 );
  TProfile effAvsev200( "effAvsev200", "A efficiency vs time;time [events];<eff(A)>",
			1000, 0, 200*1000, -1, 2 );
  TProfile effAvsev( "effAvsev", "A efficiency vs time;time [events];<eff(A)>",
		     7000, 0, 7000*1000, -1, 2 );
  TProfile effAvsdx( "effAvsdx", "A efficiency vs dx;dxCB [mm];<eff(A)>",
		     40, -0.2, 0.2, -1, 2 );
  TProfile effAvsdy( "effAvsdy", "A efficiency vs dy;dyCB [mm];<eff(A)>",
		     40, -0.2, 0.2, -1, 2 );
  TProfile effAvsncl( "effAvsncl",
		      "A efficiency vs multiplicity;clusters A * clusters B;<eff(A)>",
		     40, 0.5, 40.5, -1, 2 );
  TProfile nmvsevCB( "nmvsevCB", "CB matches vs time;time [events];CB matches",
		     7000, 0, 7000*1000, -1, 99 );

  // triplets:

  TH1I hxC( "xC", "x C;x [mm];clusters C", 100, -5, 5 );
  TH1I hyC( "yC", "y C;y [mm];clusters C", 100, -5, 5 );
  TH1I hxCi( "xCi", "x C isolated;x [mm];isolated clusters C", 100, -5, 5 );
  TH1I hyCi( "yCi", "y C isolated;y [mm];isolated clusters C", 100, -5, 5 );
  TH1I hclqCi( "clqCi", "C isolated cluster charge;cluster charge [ke];C isolatewd clusters",
	       100, 0, 50 );

  TH2I * hxxCA = new TH2I( "xxCA", "C vs A;row A;row C;clusters", nbr, -4, 4, nbr, -4, 4 );
  TH2I * hyyCA = new TH2I( "yyCA", "C vs A;col A;col C;clusters",  nbc, -4, 4,  nbc, -4, 4 );

  TH1I hdxCA( "dxCA", "Cx-Ax;x-x [mm];cluster pairs", 400, -1, 1 );
  TH1I hdyCA( "dyCA", "Cy-Ay;y-y [mm];cluster pairs", 200, -1, 1 );

  TProfile dxvsxCA( "dxvsxCA", "dx vs x C-A;x [mm];<dx> [mm]", nbr, -4, 4, -f, f );
  TProfile dxvsyCA( "dxvsyCA", "dx vs y C-A;y [mm];<dx> [mm]",  nbc, -4, 4, -f, f );
  TH1I hdyCAc( "dyCAc", "Cy-Ay, cut dx;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dyvsyCA( "dyvsyCA", "dy vs y C-A;y [mm];<dy> [mm]",  nbc, -4, 4, -f, f );
  TProfile nmvsevCA( "nmvsevCA", "CA matches vs time;time [events];CA matches",
		     7000, 0, 7000*1000, -1, 99 );

  TH1I hdx3( "dx3", "triplet dx;dx [#mum];triplets", 500, -500, 500 );
  TH1I hdy3( "dy3", "triplet dy;dy [#mum];triplets", 200, -1000, 1000 );

  TH1I hdx3c( "dx3c", "triplet dx, cut dy;dx [#mum];triplets", 400, -100, 100 );
  TH1I hdu3c( "du3c", "triplet du, cut dy;du [#mum];triplets", 400, -100, 100 );
  TH1I htx3( "tx3", "track angle x;#theta_{x} [mrad];triplets", 400, -10, 10 );

  TH1I hdu3c2( "du3c2", "triplet du, B npx 2;du [#mum];triplets B npx 2", 400, -100, 100 );
  TH1I hdu3c222( "du3c222", "triplet du, npx 2;du [#mum];triplets npx 2", 400, -100, 100 );

  TH1I hdx3ci( "dx3ci",
	       "triplet dx, cut dy, isolated;dx [#mum];B isolated triplets",
	       400, -100, 100 );
  TH1I hdx3cii( "dx3cii",
		"triplet dx, cut dy, isolated;dx [#mum];AC isolated triplets",
		400, -100, 100 );
  TH1I hdx3ciii( "dx3ciii",
		 "triplet dx, cut dy, isolated;dx [#mum];ABC isolated triplets",
		 400, -100, 100 );

  TH1I hdx3c1( "dx3c1", "triplet dx, cut dy, npx 1;dx [#mum];triplets, B npx 1",
	       400, -100, 100 );
  TH1I hdx3c2( "dx3c2", "triplet dx, cut dy, npx 2;dx [#mum];triplets, B npx 2",
	       400, -100, 100 );
  TH1I hdx3c3( "dx3c3", "triplet dx, cut dy, npx 3;dx [#mum];triplets, B npx 3",
	       400, -100, 100 );
  TH1I hdx3c4( "dx3c4", "triplet dx, cut dy, npx 4;dx [#mum];triplets, B npx 4",
	       400, -100, 100 );
  TH1I hdx3c5( "dx3c5", "triplet dx, cut dy, npx 5;dx [#mum];triplets, B npx 5",
	       400, -100, 100 );
  TH1I hdx3c6( "dx3c6", "triplet dx, cut dy, npx 6;dx [#mum];triplets, B npx 6",
	       400, -100, 100 );
  TH1I hdx3c7( "dx3c7", "triplet dx, cut dy, npx > 6;dx [#mum];triplets, B npx > 6",
	       400, -100, 100 );

  TH1I hdx3m( "dx3m", "triplet dx, x < 0;dx [#mum];triplets", 400, -100, 100 );
  TH1I hdx3p( "dx3p", "triplet dx, x > 0;dx [#mum];triplets", 400, -100, 100 );
  TH1I hdx3ct( "dx3ct", "triplet dx, cut dy, tx;dx [#mum];triplets",
	       400, -100, 100 );
  TProfile madx3vsq( "madx3vsq", "MAD(dx) vs Q;B cluster charge [ke];MAD dx [#mum]",
		     100, 0, 100, 0, 100 );
  TProfile madx3vsn( "madx3vsn",
		     "MAD(dx) vs cluster size;B cluster size [pixels];MAD dx [#mum]",
		     20, 0.5, 20.5, 0, 100 );
  TProfile madu3vsq( "madu3vsq", "MAD(dx) vs Q;B cluster charge [ke];MAD du [#mum]",
		     100, 0, 100, 0, 100 );

  TH1I hxmod( "xmod", "xavg, triplet, Landau peak;x_{avg} mod 50 [#mum];Landau peak hits",
		 50, 0, 50 );

  TH1I hxbmod( "xbmod", "xB, triplet, Landau peak;x_{B} mod 50 [#mum];Landau peak hits",
	      200, 0, 100 );

  TH2I hxbmod_clsz( "xbmod_clsz", "xB, triplet, Landau peak;cluster size;Landau peak hits",
		    8, 1, 9, 200, 0, 100 );

  TH1I hxBmodptch( "xBmodptch", "xB, triplet, Landau peak;x_{B} mod pitch [#mum];Landau peak hits",
	       50, 0, 25 );

  TH1I hxcorrB( "xcorrB", "xB, triplet, Landau peak;x_{B} mod pitch [#mum];Landau peak hits", 34, 0, 17 );
 
  TH1I hdx3cq( "dx3cq", "triplet dx, Landau peak;dx [#mum];Landau peak triplets",
	       400, -100, 100 );
  TH1I hdx3cqi( "dx3cqi",
		"triplet dx, Landau peak, isolated;dx [#mum];isolated Landau peak triplets",
		400, -100, 100 );
  TH1I hdx3cq3( "dx3cq3", "triplet dx, 3 Landau peak;dx [#mum];Landau peak triplets",
		400, -100, 100 );
  TH1I hdx3cq3i( "dx3cq3i",
		 "triplet dx, 3 Landau peak, isolated;dx [#mum];isolated Landau peak triplets",
		 400, -100, 100 );

  TH1I hdx3clsz1( "dx3clsz1", "triplet dx, 3 Landau peak [cluster size 1 events];dx [#mum];Landau peak triplets",
		400, -100, 100 );

  TH1I hdx3clsz2( "dx3clsz2", "triplet dx, 3 Landau peak [cluster size 2 events];dx [#mum];Landau peak triplets",
		 400, -100, 100 );

  TH1I hdx3clsz3( "dx3clsz3", "triplet dx, 3 Landau peak [cluster size 3 or more events];dx [#mum];Landau peak tripets", 400, -100, 100 );

  TH1I hdu3cq3( "du3cq3", "triplet du, 3 Landau peak;du [#mum];Landau peak triplets",
		400, -100, 100 );
  TH1I hdu3cq3l( "du3cq3l", "triplet du, 3 Landau peak;du [#mum];Landau peak triplets",
		400, -100, 100 );
  TH1I hdu3cq3r( "du3cq3r", "triplet du, 3 Landau peak;du [#mum];Landau peak triplets",
		400, -100, 100 );

  TH1I hdv3cq3( "dv3cq3", "triplet dv, 3 Landau peak;dv [#mum];Landau peak triplets",
		400, -100, 100 );

  TH1I hdxfl3cq3( "dxfl3cq3",
		  "triplet dx, 3 Landau peak;first-last dx [#mum];Landau peak triplets",
		  400, -100, 100 );

  TProfile dx3vsev( "dx3vsev", "dx3 vs time;trigger;<dx3> [#mum]",
		    700, 0, 7000*1000, -100, 100 );

  TProfile dx3vsx( "dx3vsx", "dx vs x;x [mm];<dx3> [#mum]", nbr, -4, 4, -100, 100 );
  TProfile dx3vsy( "dx3vsy", "dx vs y;y [mm];<dx3> [#mum]",  nbc, -4, 4, -100, 100 );

  TProfile dx3vsdx( "dx3vsdx", "MAD(dx3) vs dx C-A;C-A dx [#mum];<dx3> [#mum]",
		      100, -100, 100, -50, 50 );

  TProfile dx3vsxm( "dx3vsxm",
		    "dx vs x mod 50 um;x mod 50 [#mum];<dx3> [#mum]", 50, 0, 50, -100, 100 );
  TProfile du3vsxm( "du3vsxm",
		    "du vs x mod 50 um;x mod 50 [#mum];<du3> [#mum]", 50, 0, 50, -50, 50 );
  TProfile dv3vsxm( "dv3vsxm",
		    "dv vs x mod 50 um;x mod 50 [#mum];<dv3> [#mum]", 50, 0, 50, -50, 50 );
  TProfile dx3vsxm2( "dx3vsxm2",
		     "dx vs x mod 25 um;x mod 25 [#mum];<dx3> [#mum]", 25, 0, 25, -100, 100 );
  TProfile du3vsxm2( "du3vsxm2",
		     "du vs x mod 25 um;x mod 25 [#mum];<du3> [#mum]", 25, 0, 25, -100, 100 );

  TProfile dx3vsskwA( "dx3vsskwA",
		      "dx vs skw;A skw [rows];<dx3> [#mum]",  100, -0.25, 0.25, -100, 100 );
  TProfile dx3vsskwB( "dx3vsskwB",
		      "dx vs skw;B skw [rows];<dx3> [#mum]",  100, -0.25, 0.25, -100, 100 );
  TProfile dx3vsskwC( "dx3vsskwC",
		      "dx vs skw;C skw [rows];<dx3> [#mum]",  100, -0.25, 0.25, -100, 100 );

  TProfile du3vsskwA( "du3vsskwA",
		      "du vs A skw;A skw [rows];<du3> [#mum]",  100, -0.25, 0.25, -100, 100 );
  TProfile du3vsskwB( "du3vsskwB",
		      "du vs B skw;B skw [rows];<du3> [#mum]",  100, -0.25, 0.25, -100, 100 );
  TProfile du3vsskwC( "du3vsskwC",
		      "du vs C skw;C skw [rows];<du3> [#mum]",  100, -0.25, 0.25, -100, 100 );

  TProfile madx3vsdx( "madx3vsdx", "MAD(dx3) vs dx C-A;C-A dx [#mum];MAD dx3 [#mum]",
		      100, -100, 100, 0, 25 );

  TH1I hdx3cq3t( "dx3cq3t",
		 "triplet dx, 3 Landau peak, forward;dx [#mum];Landau peak forward triplets",
		 400, -100, 100 );
  TProfile madx3vsx( "madx3vsx", "MAD(dx3) vs x;x [#mum];MAD dx3 [#mum]", nbr, -4, 4, 0, 25 );
  TProfile madx3vsy( "madx3vsy", "MAD(dx3) vs y;y [#mum];MAD dx3 [#mum]",  nbc, -4, 4, 0, 25 );
  TProfile madx3vsxm( "madx3vsxm", "MAD(dx3) vs xmod;x mod 50 [#mum];MAD dx3 [#mum]",
		      50, 0, 50, 0, 25 );
  TProfile madu3vsxm( "madu3vsxm", "MAD(du3) vs xmod;x mod 50 [#mum];MAD du3 [#mum]",
		      50, 0, 50, 0, 25 );
  TProfile madu3vsxm2( "madu3vsxm2",
		       "MAD(du3) vs x mod 25 um;x mod 25 [#mum];MAD du3 [#mum]",
		       25, 0, 25, 0, 25 );
  TProfile madv3vsxm( "madv3vsxm", "MAD(dv3) vs xmod;x mod 50 [#mum];MAD dv3 [#mum]",
		      50, 0, 50, 0, 25 );

  TProfile etavsxmB3( "etavsxmB3", "eta vs xmod;x mod 50 [#mum];B <eta>",
		      50, 0, 50, -1.1, 1.1 );
  TProfile madx3vseta( "madx3vseta", "MAD(dx3) vs eta;eta;MAD dx3 [#mum]",
		       100, -1, 1, 0, 100 );
  TH1I hdx3cq3t2( "dx3cq3t2",
		  "triplet dx, 3 Landau peak, forward, 2-px;dx [#mum];Landau peak forward triplets",
		  400, -100, 100 );

  TH2I * hclmapB3 = new
    TH2I( "clmapB3", "linked cluster map B;col;row;B clusters on tracks",
	  nbc, 0, nbc, nbr, 0, nbr );

  TH1I hxA3( "xA3", "x A linked;x [mm];A clusters on tracks", 100, -5, 5 );
  TH1I hyA3( "yA3", "y A linked;y [mm];A clusters on tracks", 100, -5, 5 );
  TH1I hxB3( "xB3", "x B linked;x [mm];B clusters on tracks", 100, -5, 5 );
  TH1I hyB3( "yB3", "y B linked;y [mm];B clusters on tracks", 100, -5, 5 );
  TH1I hxC3( "xC3", "x C linked;x [mm];C clusters on tracks", 100, -5, 5 );
  TH1I hyC3( "yC3", "y C linked;y [mm];C clusters on tracks", 100, -5, 5 );
 
  TH1I hclszA3( "clszA3", "A cluster size on tracks;cluster size [pixels];Aclusters on tracks",
		40, 0.5, 40.5 );
 
  TH1I hclszA3cut( "clszA3cut", "A cluster size on tracks after cuts ;cluster size [pixels];Aclusters on tracks", 40, 0.5, 40.5 );
  TH1I hclszA3a("clszA3a", "A cluster size on tracks above qR ;cluster size [pixels];Aclusters on tracks", 40, 0.5, 40.5 );
  TH1I hclszB3a("clszB3a", "B cluster size on tracks above qR ; cluster size [pixels];Aclusters on tracks", 40, 0.5, 40.5);
  TH1I hclszC3a("clszC3a", "C cluster size on tracks above qR ; cluster size [pixels];Aclusters on tracks", 40, 0.5, 40.5);

  TH1I hclszB3b("clszB3b", "B cluster size on tracks below qL ; cluster size [pixels];Aclusters on tracks", 40, 0.5, 40.5);

  TH1I hnrowA3cut("nrowA3cut","A cluster size on tracks;cluster size [rows];A clusters on tracks", 20, 0, 20);
  TH1I hnrowB3cut("nrowB3cut","B cluster size on tracks;cluster size [rows];B clusters on tracks", 20, 0, 20);
  TH1I hnrowC3cut("nrowC3cut","C cluster size on tracks;cluster size [rows];C clusters on tracks", 20, 0, 20);
 
  TH1I hclszvsxm("clszvsxm", "clsz vs xm; x mod 50; cluster size [rows]",100, 0,50);

  TH1I hnrowA3( "nrowA3", "A cluster size on tracks;cluster size [rows];A clusters on tracks",
		20, 0.5, 20.5 );
  TH1I hclqA3( "clqA3", "A cluster charge on tracks;cluster charge [ke];A clusters on tracks",
	       800, 0, 80 );

  TH1I hclszB3( "clszB3", "B cluster size on tracks;cluster size [pixels];B clusters on tracks",
		40, 0.5, 40.5 );
  TH1I hclszB3cut( "clszB3cut", "B cluster size on tracks after cuts ;cluster size [pixels];B clusters on tracks",
                40, 0.5, 40.5 );
  TH1I hncolB3( "ncolB3", "B cluster size on tracks;cluster size [columns];B clusters on tracks",
		20, 0.5, 20.5 );
  TH1I hnrowB3( "nrowB3", "B cluster size on tracks;cluster size [rows];B clusters on tracks",
		20, 0.5, 20.5 );
  TH1I hclqB3( "clqB3", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",
	       160, 0, 80 );

  TH1I hncolB3cut( "ncolB3cut", "B cluster size on tracks;cluster size [columns];B clusters on tracks",
                20, 0.5, 20.5 );

  //correlation b/w cluster size and cluster charge:

  TH2I hclqvszB("clqvszB","cluster size vs cluster charge for B;cluster charge[ke];cluster size",
		160, 0, 80, 25, 0.5, 25.5);
  TH2I hclqvszA("clqvszA","cluster size vs cluster charge for A;cluster charge[ke];cluster size",
		160, 0, 80, 25, 0.5, 25.5);
  TH2I hclqvszC("clqvszC","cluster size vs cluster charge for C;cluster charge[ke];cluster size",
		160, 0, 80, 25, 0.5, 25.5);

  //correlation b/w cluster charge in  diff planes:

  TH2I hclqAvsB("clqAvsB","cluster charge in A vs B; cluster charge A [ke];cluster charge B [ke]",
		160, 0, 80, 160, 0, 80);
  TH2I hclqAvsC("clqAvsC","cluster charge in A vs C; cluster charge A [ke];cluster charge C [ke]",
		160, 0, 80, 160, 0, 80);
  TH2I hclqBvsC("clqBvsC","cluster charge in B vs C; cluster charge B [ke];cluster charge C [ke]",
		160, 0, 80, 160, 0, 80);

  TH2I hclszAvsB("clszAvsB","cluster size in A vs B; cluster size A;cluster size B",
		 25, 0.5, 25.5, 25, 0.5, 25.5);
  TH2I hclszAvsC("clszAvsC","cluster size in B vs C; cluster size A;cluster size C",
		 25, 0.5, 25.5, 25, 0.5, 25.5);
  TH2I hclszBvsC("clszBvsC","cluster size in B vs C; cluster size B;cluster size C",
		 25, 0.5, 25.5, 25, 0.5, 25.5);

  TH1I hclqB3i( "clqB3i", "B cluster charge on tracks;cluster charge [ke];B clusters on tracks",
	       80, 0, 20 );
  TH1I hclqB3n( "clqB3n",
		"B cluster charge on tracks, npx < 4;cluster charge [ke];B clusters on tracks, npx < 4",
		160, 0, 80 );

  TH1I hclszC3( "clszC3", "C cluster size on tracks;cluster size [pixels];C clusters on tracks",
		40, 0.5, 40.5 );
  TH1I hclszC3cut("clszC3cut", "C cluster size on tracks after cuts;cluster size [pixels];C clusters on tracks",40, 0.5, 40.5 );
  TH1I hnrowC3( "nrowC3", "C cluster size on tracks;cluster size [rows];C clusters on tracks",
		20, 0.5, 20.5 );
  TH1I hclqC3( "clqC3", "C cluster charge on tracks;cluster charge [ke];C clusters on tracks",
	       160, 0, 80 );

  TProfile nrowvsxmB3( "nrowvsxmB3",
		       "B rows vs xmod;x mod 50 [#mum];<B cluster size [rows]>",
		       50, 0, 50, 0.5, 10.5 );
 
  TProfile nrowvsymB3 ("nrowvsymB3", 
		       "B rows vs ymod; y mod  300; <B cluster size [rows]>", 
		       100, 0, 300, 0.5, 10.5);
  TProfile ncolvsxmB3 ("ncolvsxmB3",
		       "B cols vs xmod; x mod 50; <B cluster size [cols]>",
		       50, 0, 50, 0.5, 10.5);
  TProfile ncolvsymB3 ("ncolvsymB3",
		       "B cols vs ymod; y mod  300; <B cluster size [cols]>",
		       100, 0, 300, 0.5, 10.5); 

  TProfile clqvsxmB3( "clqvsxmB3",
		      "B cluster charge vs xmod;x mod 50 [#mum];<B cluster charge [ke]>",
		      50, 0, 50, 0, 50 );
 
  TH1I hetaA3( "etaA3", "A cluster eta;eta;A 2-pix clusters on tracks",
	       100, -1, 1 );
  TH1I hetaB3( "etaB3", "B cluster eta;eta;B 2-pix clusters on tracks",
	       100, -1, 1 );
  TH1I hetaC3( "etaC3", "C cluster eta;eta;C 2-pix clusters on tracks",
	       100, -1, 1 );

  TH1I hpxqA3( "pxqA3", "A pixel charge;pixel charge [ke];A pixels on tracks",
	       100, 0, 20 );
  TH1I hpxqB3( "pxqB3", "B pixel charge;pixel charge [ke];B pixels on tracks",
	       100, 0, 20 );
  TH1I hpxqC3( "pxqC3", "C pixel charge;pixel charge [ke];C pixels on tracks",
	       100, 0, 20 );

  TH1I hpxq1stB3( "pxq1stB3", "B 1st pixel charge;pixel charge [ke];B `st pixels on tracks",
		  100, 0, 20 );
  TH1I hpxq2ndB3( "pxq2ndB3", "B 2nd pixel charge;pixel charge [ke];B `st pixels on tracks",
		  100, 0, 20 );
  //TH1I hpxq3rdB3( "pxq2ndB3", "B 3rd pixel charge;pixel charge [ke];B `st pixels on tracks",
  //              100, 0, 20 );

  TH1I hrmsA3( "rmsA3", "A cluster rms;rms [rows];A clusters on tracks", 100, 0, 3 );
  TH1I hskwA3( "skwA3", "A cluster skw;skw [rows];A clusters on tracks", 100, -0.25, 0.25 );

  TH1I hrmsC3( "rmsC3", "C cluster rms;rms [rows];C clusters on tracks", 100, 0, 3 );
  TH1I hskwC3( "skwC3", "C cluster skw;skw [rows];C clusters on tracks", 100, -0.25, 0.25 );

  TH1I hrmsB3( "rmsB3", "B cluster rms;rms [rows];B clusters on tracks", 100, 0, 3 );
  TH1I hskwB3( "skwB3", "B cluster skw;skw [rows];B clusters on tracks", 100, -0.25, 0.25 );
  TProfile rmsvsxmB3( "rmsvsxmB3",
		      "B rms vs xmod;x mod 50 [#mum];<B cluster rms> [rows]",
		      50, 0, 50, 0, 5 );
  TProfile skwvsxmB3( "skwvsxmB3",
		      "B skw vs xmod;x mod 50 [#mum];<B cluster skw> [rows]",
		      50, 0, 50, -0.1, 0.1 );

  TProfile effvsdxy( "effvsdxy",
		     "DUT efficiency vs triplet dxy;xy match radius [mm];DUT efficiency",
		     1000, 0, 10, -0.1, 1.1 );

  TProfile2D * effvsxy =
    new TProfile2D( "effvsxy",
		    "DUT efficiency map;x [mm];y[mm];DUT efficiency",
		    nbr, -4, 4, nbc, -4, 4, -0.1, 1.1 );
  TProfile effvsx( "effvsx", "eff vs x;x [mm];DUT efficiency",
		   nbr, -4, 4, -0.1, 1.1 );
  TProfile effvsy( "effvsy", "eff vs y;y [mm];DUT efficiency",
		   nbc, -4, 4, -0.1, 1.1 );
  TProfile effvsxm( "effvsxm", "eff vs x mod 50;x mod 50 [#mum];DUT efficiency",
		    50, 0, 50, -0.1, 1.1 );

  TProfile effvsev( "effvsev", "eff vs time;trigger;DUT efficiency",
		    7000, 0, 7000*1000, -0.1, 1.1 );
  TProfile effvsiev( "effvsiev", "eff vs event;event mod 200;DUT efficiency",
		     100, -0.5, 199.5, -0.1, 1.1 );
  TProfile effvsqA( "effvsqA", "eff vs charge A;cluster charge A [ke];DUT efficiency",
		    100, 0, 100, -0.1, 1.1 );
  TProfile effvstxy( "effvstxy", "eff vs angle;dxy CA [mm];DUT efficiency",
		     100, 0, 0.2, -0.1, 1.1 );

  //-----------------------------------------------------------------------------

  TFile *bias = new TFile( Form( "bias-r%i.root", run ), "READ");

  TH1I *xbBias = new TH1I();

  if( bias->IsOpen() )
    {
    cout << "bias file open!" << endl;
    xbBias = (TH1I*) bias->Get("xBias");
    xbBias->SetName("xbBias");    
    }
  else
    cout << "no bias file included" << endl;

   histoFile -> cd();
   xbBias -> Write();

  //-------------------------------------------------------------------------------

  const double log10 = log(10);
  string ADD {"A"}; // added flag

  // PSI 285 mu:

  double qL = 14; // peak 22 ke
  double qR = 30;

  double qLB = qL;
  double qRB = qR;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s0 = ts.tv_sec; // seconds since 1.1.1970
  long f0 = ts.tv_nsec; // nanoseconds

  clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ); // sum of threads time
  time_t t0 = ts.tv_sec; // seconds
  long n0 = ts.tv_nsec; // nanoseconds

  double zeit1 = 0;
  double zeit2 = 0;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Read and process the run for each plane in parallel:

  list < vector < cluster > > evlistA;
  list < vector < cluster > > evlistB;
  list < vector < cluster > > evlistC;

  //#pragma omp sections // test, not parallel
#pragma omp parallel sections
  {
#pragma omp section
    {
      //cout << "A " << sched_getcpu() << endl << flush; // changes

      evlistA = oneplane( A, runnum, Nev, geo );

    } //omp section

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // B:

#pragma omp section
    {
      //cout << "B " << sched_getcpu() << endl << flush; // different from A

      evlistB = oneplane( B, runnum, Nev, geo );

    } // omp section

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // C:

#pragma omp section
    {
      //cout << "C " << sched_getcpu() << endl << flush; // different from A and B

      evlistC = oneplane( C, runnum, Nev, geo );

    } // omp section

  } // omp parallel sections

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s2 = ts.tv_sec; // seconds since 1.1.1970
  long f2 = ts.tv_nsec; // nanoseconds
  zeit1 += s2 - s0 + ( f2 - f0 ) * 1e-9; // read and cluster

  clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ); // sum threads
  time_t t9 = ts.tv_sec; // seconds
  long n9 = ts.tv_nsec; // nanoseconds

  cout << "time " << zeit1 << " s"
       << " (sum threads " << t9 - t0 + ( n9 - n0 ) * 1e-9 << " s)"
       << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // loop over events, correlate planes:

  auto evA = evlistA.begin();
  auto evB = evlistB.begin();
  auto evC = evlistC.begin();

  auto evinfoA = infoA.begin();
  auto evinfoB = infoB.begin();
  auto evinfoC = infoC.begin();

  unsigned nev = 0;
  uint64_t prevtimeA = 0;
  uint64_t prevtimeB = 0;
  uint64_t prevtimeC = 0;

  int dsyncAB = 0;
  int dsyncCB = 0;
  int dsyncCA = 0;

  cout << endl << "tracking ev" << flush;

  for( ; evinfoB != infoB.end() &&
	 evA != evlistA.end() && evB != evlistB.end() && evC != evlistC.end();
       ++evinfoA, ++evinfoB, ++evinfoC, ++evA, ++evB, ++evC ) {

    vector <cluster> vclA = *evA;
    vector <cluster> vclB = *evB;
    vector <cluster> vclC = *evC;

    ++nev;
    //if( nev%10000 == 0 )
      //cout << " " << nev << flush;

    int iev = evinfoA->iev;

    // if( evinfoA->filled == ADD ) cout << "  " << iev << " add A" << endl;
    //if( evinfoB->filled == ADD ) cout << "  " << iev << " add B" << endl;
    //if( evinfoC->filled == ADD ) cout << "  " << iev << " add C" << endl;

    int64_t dtA = evinfoA->evtime - prevtimeA;
    prevtimeA = evinfoA->evtime;

    int64_t dtB = evinfoB->evtime - prevtimeB;
    prevtimeB = evinfoB->evtime;
    if( iev > 1 && dtB > 0 ) // added events have same time
      hdt.Fill( log(dtB/40e6) / log10 ); // MHz -> s

    int64_t dtC = evinfoC->evtime - prevtimeC;
    prevtimeC = evinfoC->evtime;

    int64_t ddtAB = dtA - dtB;
    int64_t ddtCB = dtC - dtB;
    int64_t ddtCA = dtC - dtA;
    if( iev > 1 && dtB > 0 ) {
      hddtAB.Fill( ddtAB );
      hddtCB.Fill( ddtCB );
      hddtCA.Fill( ddtCA );
      ddtvsdtAB.Fill( log(dtB/40e6) / log10, ddtAB );
      ddtvsevAB.Fill( iev, ddtAB );
      ddtvsevCB.Fill( iev, ddtCB );
      ddtvsevCA.Fill( iev, ddtCA );
    }
    if( abs( ddtAB ) > 999 ) ++dsyncAB;
    if( abs( ddtCB ) > 999 ) ++dsyncCB;
    if( abs( ddtCA ) > 999 ) ++dsyncCA;

    // 2207 AB unsync 7150
    /*
    if( abs( ddtAB ) > 999 )
      cout << "  " << iev
	   << "  AB ddt " << ddtAB
	   << "  A dt " << dtA
	   << "  B dt " << dtB
	   << endl;

    if( abs( ddtCB ) > 999 )
      cout << "  " << iev
	   << "  CB ddt " << ddtCB
	   << "  C dt " << dtC
	   << "  B dt " << dtB
	   << endl;*/

    if( evinfoA->skip ) skipvsevA.Fill( iev );
    if( evinfoB->skip ) skipvsevB.Fill( iev );
    if( evinfoC->skip ) skipvsevC.Fill( iev );


    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // A-B cluster correlations:

    double ptchc = 0.100; // [mm] col size
    double ptchr = 0.025; // [mm] row size
    if( abs(geo) == 50 ) {
      ptchc = 0.050;
      ptchr = 0.050;
    }
    else if( geo == 16 ) {
      ptchc = 0.150;
      ptchr = 0.050/3;
    }
    int nm = 0;

    for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ) {

      double xA = cA->row * ptchr - 4.0 - alignxA; // rot90 Dreimaster
      double yA = cA->col * ptchc - 3.9 - alignyA; // 100 um px

      double xAr = xA*cfA - yA*sfA - xA*aligntA;
      double yAr = xA*sfA + yA*cfA;

      //double xArm = xAr + 4.0;

      hxA.Fill( xAr );
      hyA.Fill( yAr );
      if( cA->iso ) {
	hxAi.Fill( xAr );
	hyAi.Fill( yAr );
	hclqAi.Fill( cA->q );
      }

      for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

	double xB = cB->row * ptchr - 4.0;
	double yB = cB->col * ptchc - 3.9;

	hxxAB->Fill( xAr, xB );
	hyyAB->Fill( yAr, yB );

	double dxAB = xAr - xB;
	double dyAB = yAr - yB;

	if( cA->q > qL  && cA->q < qR &&
	    cB->q > qLB && cB->q < qRB &&
	    cA->iso && cB->iso ) {	  

	  hdxAB.Fill( dxAB );
	  hdyAB.Fill( dyAB );
	  dxvsxAB.Fill( xB, dxAB );
	  dxvsyAB.Fill( yB, dxAB );

	}

	hdxvsev->Fill( iev, dxAB );

	if( fabs( dxAB ) < 40*3E-3 * 5/pp + 0.020 &&
	    fabs( dyAB ) < 40*3E-3 * 5/pp + 0.100 )
	  ++nm;

	// extrapolate to C for eff:

	double xatC = xB + dxAB; // equal spacing
	double yatC = yB + dyAB; // equal spacing

	double dmin = 99;

	for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ) {

	  double xC = cC->row * ptchr - 4.0 - alignxC; // rot90 Dreimaster
	  double yC = cC->col * ptchc - 3.9 - alignyC; // down

	  double xCr = xC*cfC - yC*sfC - xC*aligntC;
	  double yCr = xC*sfC + yC*cfC;

	  double dx = xCr - xatC;
	  double dy = yCr - yatC;

	  if( cA->q > qL  && cA->q < qR &&
	      cB->q > qLB && cB->q < qRB &&
	      cA->iso && cB->iso ) {
	   
	    hdxABC.Fill( dx );
	    hdyABC.Fill( dy );
	    dxvsxABC.Fill( xCr, dx );
	    dxvsyABC.Fill( yCr, dx );

	  }

	  double dxy = sqrt( dx*dx + dy*dy );
	  if( dxy < dmin ) dmin = dxy;

	} // clusters C

	if( cA->q > qL  && cA->q < qR &&
	    cB->q > qLB && cB->q < qRB &&
	    fabs( dxAB ) < 40*3E-3 * 5/pp + 0.020 &&
	    fabs( dyAB ) < 40*3E-3 * 5/pp + 0.100 &&
	    abs(ddtAB) < 400 &&
	    abs(ddtCB) < 400 &&
	    abs(ddtCA) < 400 &&
	    evinfoC->filled != ADD &&
	    !evinfoC->skip &&
	    cA->iso && cB->iso ) {

	  hdminC.Fill( dmin );

	  int eff = 1;
	  if( dmin > 0.5 ) eff = 0;

	  effCvsx.Fill( xatC, eff );
	  effCvsy.Fill( yatC, eff );
	  effCvsxy->Fill( xatC, yatC, eff );
	  effCvsev200.Fill( iev, eff );
	  effCvsev.Fill( iev, eff );
	  effCvsdx.Fill( dxAB, eff );
	  effCvsdy.Fill( dyAB, eff );
	  effCvsncl.Fill( vclA.size()*vclB.size(), eff );

	}

      } // clusters B

    } // cl A

    nmvsevAB.Fill( iev, nm );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // C-B cluster correlations:

    nm = 0;

    for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

      double xB = cB->row * ptchr - 4.0;
      double yB = cB->col * ptchc - 3.9;
     
      hxB.Fill( xB );
     
      hyB.Fill( yB );
      if( cB->iso ) {
	hxBi.Fill( xB );
	hyBi.Fill( yB );
	hclqBi.Fill( cB->q );
      }

      for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ) {

	double xC = cC->row * ptchr - 4.0 - alignxC; // rot90 Dreimaster
	double yC = cC->col * ptchc - 3.9 - alignyC; // down

	double xCr = xC*cfC - yC*sfC - xC*aligntC;
	double yCr = xC*sfC + yC*cfC;

	hxxCB->Fill( xB, xCr );
	hyyCB->Fill( yB, yCr );

	double dxCB = xCr - xB;
	double dyCB = yCr - yB;

	if( cC->q > qL  && cC->q < qR &&
	    cB->q > qLB && cB->q < qRB &&
	    cC->iso && cB->iso ) {

	  hdxCB.Fill( dxCB );
	  hdyCB.Fill( dyCB );
	  dxvsxCB.Fill( xB, dxCB );
	  dxvsyCB.Fill( yB, dxCB );

	}

	if( fabs( dxCB ) < 40*3E-3 * 5/pp + 0.020 &&
	    fabs( dyCB ) < 40*3E-3 * 5/pp + 0.100 )
	  ++nm;

	// extrapolate to A for eff:

	double xatA = xB + dxCB; // equal spacing
	double yatA = yB + dyCB; // equal spacing

	double dmin = 99;

	for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ) {

	  double xA = cA->row * ptchr - 4.0 - alignxA; // rot90 Dreimaster
	  double yA = cA->col * ptchc - 3.9 - alignyA; // down

	  double xAr = xA*cfA - yA*sfA - xA*aligntA;
	  double yAr = xA*sfA + yA*cfA;

	  double dx = xAr - xatA;
	  double dy = yAr - yatA;

	  if( cC->q > qL  && cC->q < qR &&
	      cB->q > qLB && cB->q < qRB &&
	      cC->iso && cB->iso ) {

	    hdxCBA.Fill( dx );
	    hdyCBA.Fill( dy );
	    dxvsxCBA.Fill( xAr, dx );
	    dxvsyCBA.Fill( yAr, dx );

	  }

	  double dxy = sqrt( dx*dx + dy*dy );
	  if( dxy < dmin ) dmin = dxy;

	} // clusters A

	if( cC->q > qL  && cC->q < qR &&
	    cB->q > qLB && cB->q < qRB &&
	    fabs( dxCB ) < 40*3E-3 * 5/pp + 0.020 &&
	    fabs( dyCB ) < 40*3E-3 * 5/pp + 0.100 &&
	    abs(ddtAB) < 400 &&
	    abs(ddtCB) < 400 &&
	    abs(ddtCA) < 400 &&
	    evinfoA->filled != ADD &&
	    !evinfoA->skip &&
	    cC->iso && cB->iso ) {

	  hdminA.Fill( dmin );

	  int eff = 1;
	  if( dmin > 0.5 ) eff = 0;

	  effAvsx.Fill( xatA, eff );
	  effAvsy.Fill( yatA, eff );
	  effAvsxy->Fill( xatA, yatA, eff );
	  effAvsev200.Fill( iev, eff );
	  effAvsev.Fill( iev, eff );
	  effAvsdx.Fill( dxCB, eff );
	  effAvsdy.Fill( dyCB, eff );
	  effAvsncl.Fill( vclA.size()*vclB.size(), eff );

	}

      } // clusters C

    } // cl B

    nmvsevCB.Fill( iev, nm );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // A-C cluster correlations:

    nm = 0;

    for( vector<cluster>::iterator cC = vclC.begin(); cC != vclC.end(); ++cC ) {

      double xC = cC->row * ptchr - 4.0 - alignxC;
      double yC = cC->col * ptchc - 3.9 - alignyC;

      double xCr = xC*cfC - yC*sfC - xC*aligntC;
      double yCr = xC*sfC + yC*cfC;
      //double xCrm = xCr + 4.0;

      hxC.Fill( xCr );
      // if(geo == 16)
      // xCmod -> Fill(fmod(xCr,17));
      // if(geo == 25)
      // xCmod -> Fill(fmod(xCr,25));
      // if(geo == 50)
      // xCmod ->Fill(fmod(xCr,50));

      hyC.Fill( yCr );
      if( cC->iso ) {
	hxCi.Fill( xCr );
	hyCi.Fill( yCr );
	hclqCi.Fill( cC->q );
      }

      double etaC = -2;
      if( cC->size == 2 ) {
	double q0 = cC->vpix[0].q;
	double q1 = cC->vpix[1].q;
	etaC = (q1-q0)/(q1+q0);
      }

      int rowmin = 999;
      int rowmax = 0;
      int colmin = 999;
      int colmax = 0;
      vector <double> qrow(480); // initialized to zero

      for( int ipx = 0; ipx < cC->size; ++ipx ) {

	int row = cC->vpix[ipx].row;
	if( row < rowmin ) rowmin = row;
	if( row > rowmax ) rowmax = row;

	int col = cC->vpix[ipx].col;
	if( col < colmin ) colmin = col;
	if( col > colmax ) colmax = col;

	qrow[row] += cC->vpix[ipx].q; // project cluster onto rows

      }

      int nrowC = rowmax-rowmin+1;
      int ncolC = colmax-colmin+1;

      double sumq = 0;
      double sumrow2 = 0;
      double sumrow3 = 0;

      for( int row = rowmin; row <= rowmax; ++row ) {

	double q = qrow[row];
	sumq += q;

	double drow = row - cC->row; // for central moments
	sumrow2 += drow*drow*q;
	sumrow3 += drow*drow*drow*q;

      } // rows

      double rmsrowC = 0;
      double skwrowC = 0;
      if( nrowC > 1 ) {
	rmsrowC = sqrt( sumrow2/sumq );
	skwrowC = sumrow3/sumq*8/nrowC/nrowC/nrowC; // normalized 3rd moment
      }

      // first-last:

      double flC = ( rowmin*qrow[rowmin] + rowmax*qrow[rowmax] ) / ( qrow[rowmin] + qrow[rowmax] );
      double xflC = flC * ptchr - 4.0 - alignxC;
      double xflCr = xflC*cfC - yC*sfC - xC*aligntC; //rot

      // A:

      for( vector<cluster>::iterator cA = vclA.begin(); cA != vclA.end(); ++cA ) {

	double xA = cA->row * ptchr - 4.0 - alignxA; // rot90 Dreimaster
	double yA = cA->col * ptchc - 3.9 - alignyA; // down

	double xAr = xA*cfA - yA*sfA - xA*aligntA;
	double yAr = xA*sfA + yA*cfA;

	hxxCA->Fill( xAr, xCr );
	hyyCA->Fill( yAr, yCr );

	double dxCA = xCr - xAr;
	double dyCA = yCr - yAr;
	double dxyCA = sqrt( dxCA*dxCA + dyCA*dyCA );
	double txCA = dxCA/  40; // angle
	hdxCA.Fill( dxCA );
	hdyCA.Fill( dyCA );

	//cut on track angle:

	if( fabs( dxCA ) > 3 * 40E-3 * 5/pp + 0.02 ) continue; // includes beam divergence: +-5 sigma

	hdyCAc.Fill( dyCA );
	dyvsyCA.Fill( yAr, dyCA );

	if( fabs( dyCA ) > 3 * 40E-3 * 5/pp + 0.1 ) continue; // [mm]

	++nm;

	dxvsyCA.Fill( yAr, dxCA );
	dxvsxCA.Fill( xAr, dxCA ); // linear trend in run 392, 403: acceptance and beam divergence ?

	double xavg = 0.5 * ( xAr + xCr );
	double yavg = 0.5 * ( yAr + yCr );

	double xmod = fmod( xavg + 8.0125, 0.05 ); // [mm] 0..0.05, like RD53A (8.7.2018)
	double xmod2 = fmod( xavg + 8.0125, 0.025 ); // [mm] 0..0.025

	double ymod = fmod( yavg + 8.0125, 0.3 );
 
	if( abs(geo) == 50 )
	  xmod = fmod( xavg + 8.025, 0.050 ); // [mm] 0..0.05
	else if( geo == 16 )
	  xmod = fmod( xavg + 8.025, 0.100 ); // [mm] 0..0.05

	double etaA = -2;
	if( cA->size == 2 ) {
	  double q0 = cA->vpix[0].q;
	  double q1 = cA->vpix[1].q;
	  etaA = (q1-q0)/(q1+q0);
	}

	int rowmin = 999;
	int rowmax = 0;
	int colmin = 999;
	int colmax = 0;
	vector <double> qrow(480); // initialized to zero

	for( int ipx = 0; ipx < cA->size; ++ipx ) {

	  int row = cA->vpix[ipx].row;
	  if( row < rowmin ) rowmin = row;
	  if( row > rowmax ) rowmax = row;

	  int col = cA->vpix[ipx].col;
	  if( col < colmin ) colmin = col;
	  if( col > colmax ) colmax = col;

	  qrow[row] += cA->vpix[ipx].q; // project cluster onto rows

	}

	int nrowA = rowmax-rowmin+1;
	int ncolA = colmax-colmin+1;

	double sumq = 0;
	double sumrow2 = 0;
	double sumrow3 = 0;

	for( int row = rowmin; row <= rowmax; ++row ) {

	  double q = qrow[row];
	  sumq += q;

	  double drow = row - cA->row; // for central moments
	  sumrow2 += drow*drow*q;
	  sumrow3 += drow*drow*drow*q;

	} // rows

	double rmsrowA = 0;
	double skwrowA = 0;
	if( nrowA > 1 ) {
	  rmsrowA = sqrt( sumrow2/sumq );
	  skwrowA = sumrow3/sumq*8/nrowA/nrowA/nrowA; // normalized 3rd moment
	}

	// first-last:

	double flA = ( rowmin*qrow[rowmin] + rowmax*qrow[rowmax] ) / ( qrow[rowmin] + qrow[rowmax] );
	double xflA = flA * ptchr - 4.0 - alignxA;
	double xflAr = xflA*cfA - yA*sfA - xA*aligntA; // rot

	double xflavg = 0.5 * ( xflAr + xflCr ); // average (at B)

	// B clusters:

	double pdmin = 99;

	for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

	  double xB = cB->row * ptchr - 4.0;
	  double yB = cB->col * ptchc - 3.9;

	  double etaB = -2;
	  if( cB->size == 2 ) {
	    double q0 = cB->vpix[0].q;
	    double q1 = cB->vpix[1].q;
	    etaB = (q1-q0)/(q1+q0);
	  }

	  int rowmin = 999;
	  int rowmax = 0;
	  int colmin = 999;
	  int colmax = 0;
	  vector <double> qrow(480); // initialized to zero

	  for( int ipx = 0; ipx < cB->size; ++ipx ) {

	    int row = cB->vpix[ipx].row;
	    if( row < rowmin ) rowmin = row;
	    if( row > rowmax ) rowmax = row;

	    int col = cB->vpix[ipx].col;
	    if( col < colmin ) colmin = col;
	    if( col > colmax ) colmax = col;

	    qrow[row] += cB->vpix[ipx].q; // project cluster onto rows

	  }

	  int nrowB = rowmax-rowmin+1;
	  int ncolB = colmax-colmin+1;

	  double sumq = 0;
	  double sumrow2 = 0;
	  double sumrow3 = 0;

	  for( int row = rowmin; row <= rowmax; ++row ) {

	    double q = qrow[row];
	    sumq += q;

	    double drow = row - cB->row; // for central moments
	    sumrow2 += drow*drow*q;
	    sumrow3 += drow*drow*drow*q;

	  } // rows

	  double rmsrowB = 0;
	  double skwrowB = 0;
	  if( nrowB > 1 ) {
	    rmsrowB = sqrt( sumrow2/sumq );
	    skwrowB = sumrow3/sumq*8/nrowB/nrowB/nrowB; // normalized 3rd moment
	  }

	  // first-last:

	  double flB =
	    ( rowmin*qrow[rowmin] + rowmax*qrow[rowmax] ) /
	    ( qrow[rowmin] + qrow[rowmax] );
	  double xflB = flB * ptchr - 4.0; // no align for B
	  double dxfl3 = xflB - xflavg;

	  // triplet residual:

	  double dx3 = xB - xavg;
	 
	  double du3 = dx3;

	  if( run >= 4655 && run <= 4672 ) { // du3vsxm->Fit("pol9")
	    du3 -= -0.029*skwrowA; // dx3vsskwA->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.007*skwrowB; // dx3vsskwB->Fit("pol1","","",-0.1,0.08)
	    du3 -= -0.016*skwrowC;
	  }

	  if( run == 4673 ){
	    du3 -=  0.029*skwrowA;
	    du3 -=  -0.031*skwrowB;
	    du3 -=  -0.008*skwrowC;
	  }
	  if( run == 4674 ){
	    du3 -=  0.002*skwrowA;
	    du3 -=  -0.019*skwrowB;
	    du3 -=  -0.006*skwrowC;
	  }
	  if( run == 4675 ){
	    du3 -=  0.004*skwrowA;
	    du3 -=  -0.023*skwrowB;
	    du3 -=  -0.028*skwrowC;
	  }
	  if( run == 4676 ){
	    du3 -=  -0.005*skwrowA;
	    du3 -=  -0.006*skwrowB;
	    du3 -=  -0.005*skwrowC;
	  }
	  if( run == 4677 ){
	    du3 -=  -0.016*skwrowA;
	    du3 -=  0.032*skwrowB;
	    du3 -=  -0.023*skwrowC;
	  }
	  if( run == 4678 ){
	    du3 -=  -0.029*skwrowA;
	    du3 -=  0.034*skwrowB;
	    du3 -=  -0.025*skwrowC;
	  }

	  if( run == 4679 ){
	    du3 -=  -0.025*skwrowA;
	    du3 -=  0.025*skwrowB;
	    du3 -=  -0.024*skwrowC;
	  }
	  if( run == 4680 ){
	    du3 -=  -0.006*skwrowA;
	    du3 -=  -0.013*skwrowB;
	    du3 -=  -0.013*skwrowC;
	  }
	  if( run == 4681 ){
	    du3 -=  0.006*skwrowA;
	    du3 -=  -0.042*skwrowB;
	    du3 -=  -0.010*skwrowC;
	  }
	  if( run == 4682 ){
	    du3 -=  0.016*skwrowA;
	    du3 -=  -0.055*skwrowB;                  
	    du3 -=  -0.012*skwrowC;
	  }
	  if( run == 4683 ){
	    du3 -=  0.008*skwrowA;
	    du3 -=  -0.025*skwrowB;
	    du3 -=  -0.012*skwrowC;
	  }
	  if( run == 4684 ){
	    du3 -=  0.008*skwrowA;
	    du3 -=  -0.027*skwrowB;
	    du3 -=  -0.018*skwrowC;
	  }
	    if( run == 4685 ){
	      du3 -=  0.011*skwrowA;
	      du3 -=  -0.010*skwrowB;
	      du3 -=  -0.015*skwrowC;
	    }
	    if( run == 4686 ){
	      du3 -=  -0.009*skwrowA;
	      du3 -=  0.012*skwrowB;
	      du3 -=  -0.018*skwrowC;
	    }
	    if( run == 4687 ){
	      du3 -=  -0.017*skwrowA;
	      du3 -=  0.017*skwrowB;
	      du3 -=  -0.015*skwrowC;
	    }
	    if( run == 4688 ){
	      du3 -=  -0.014*skwrowA;
	      du3 -=  -0.002*skwrowB;
	      du3 -=  -0.021*skwrowC;
	    }
	    if( run == 4689 ){
	      du3 -=  -0.003*skwrowA;
	      du3 -=  -0.024*skwrowB;
	      du3 -=  -0.018*skwrowC;
	    }

	    if( run == 4690 ){
	      du3 -=  0.001*skwrowA;
	      du3 -=  -0.032*skwrowB;
	      du3 -=  -0.019*skwrowC;
	    }
	    if( run == 4691 ){
	      du3 -=  0.006*skwrowA;
	      du3 -=  -0.035*skwrowB;
	      du3 -=  -0.002*skwrowC;
	    }
	    if( run == 4692 ){
	      du3 -=  0.022*skwrowA;
	      du3 -=  -0.064*skwrowB;
	      du3 -=  -0.014*skwrowC;
	    }
	    if( run == 4693 ){
	      du3 -=  0.032*skwrowA;
	      du3 -=  -0.087*skwrowB;
	      du3 -=  0.018*skwrowC;
	    }
	    if( run == 4694 ){
	      du3 -=  0.044*skwrowA;
	      du3 -=  -0.079*skwrowB;
	      du3 -=  0.019*skwrowC;
	    }
	    if( run == 4695 ){
	      du3 -=  0.044*skwrowA;
	      du3 -=  -0.093*skwrowB;
	      du3 -=  0.027*skwrowC;
	    }

	  
	  if( run == 4698 ) {
	    du3 -= -0.008*skwrowA; // dx3vsskwA->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.014*skwrowB; // dx3vsskwB->Fit("pol1","","",-0.1,0.1)
	    du3 -= -0.016*skwrowC;
	  }

	  if( run == 4701 ) {
	    du3 -=  0.002*skwrowA; // dx3vsskwA->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.011*skwrowB; // dx3vsskwB->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.002*skwrowC;
	  }

	  if( run == 4702 ) {
            du3 -=  0.010*skwrowA; // dx3vsskwA->Fit("pol1","","",-0.1,0.1)
	    du3 -= 0.005*skwrowB; // dx3vsskwB->Fit("pol1","","",-0.1,0.1)
            du3 -=  0.009*skwrowC;
          }
	  if( run == 4703 ) {
	    du3 -=  0.017*skwrowA; // dx3vsskwA->Fit("pol1","","",-0.1,0.1)
	    du3 -= -0.008*skwrowB; // dx3vsskwB->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.019*skwrowC;
	  }

	  if( run == 4704 ) {
            du3 -=  0.024*skwrowA; // dx3vsskwA->Fit("pol1","","",-0.1,0.1)                                                                               
            du3 -= -0.023*skwrowB; // dx3vsskwB->Fit("pol1","","",-0.1,0.1)                                                  
	    du3 -=  0.025*skwrowC;
          }
	  if( run == 4704 ){
	    du3 -=  0.024*skwrowA;
	    du3 -=  -0.023*skwrowB;
	    du3 -=  0.025*skwrowC;
	  }
	  if( run == 4705 ){
	    du3 -=  0.027*skwrowA;
	    du3 -=  -0.036*skwrowB;
	    du3 -=  0.029*skwrowC;
	  }
	  if( run == 4706 ){
	    du3 -=  0.026*skwrowA;
	    du3 -=  -0.047*skwrowB;
	    du3 -=  0.009*skwrowC;
	  }
	  if( run == 4707 ){
	    du3 -=  0.023*skwrowA;
	    du3 -=  -0.048*skwrowB;
	    du3 -=  0.002*skwrowC;
	  }
	  if( run == 4708 ){
	    du3 -=  0.017*skwrowA;
	    du3 -=  -0.044*skwrowB;
	    du3 -=  -0.002*skwrowC;
	  }
	  if( run == 4709 ){
	    du3 -=  0.012*skwrowA;
	    du3 -=  -0.033*skwrowB;
	    du3 -=  -0.004*skwrowC;
	  }
	  if( run == 4710 ) {
	    du3 -=  0.001*skwrowA; // dx3vsskwA->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.011*skwrowB; // dx3vsskwB->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.005*skwrowC;
	  }
	  if( run == 4711 ){
	    du3 -=  -0.005*skwrowA;
	    du3 -=  0.009*skwrowB;
	    du3 -=  -0.0004*skwrowC;
	  }
	  if( run == 4712 ){
	    du3 -=  -0.006*skwrowA;
	    du3 -=  0.004*skwrowB;
	    du3 -=  -0.003*skwrowC;
	  }
	  if( run == 4713 ){
	    du3 -=  -0.006*skwrowA;
	    du3 -=  -0.018*skwrowB;
	    du3 -=  -0.003*skwrowC;
	  }
	  if( run == 4714 ){
	    du3 -=  -0.0008*skwrowA;
	    du3 -=  -0.027*skwrowB;
	    du3 -=  0.011*skwrowC;
	  }
	  if( run == 4715 ){
	    du3 -=  0.002*skwrowA;
	    du3 -=  -0.014*skwrowB;
	    du3 -=  0.022*skwrowC;
	  }
	  if( run == 4716 ){
	    du3 -=  0.008*skwrowA;
	    du3 -=  -0.013*skwrowB;
	    du3 -=  0.013*skwrowC;
	  }
	  if( run == 4717 ){
	    du3 -=  0.007*skwrowA;
	    du3 -=  -0.026*skwrowB;
	    du3 -=  0.011*skwrowC;
	  }
	  if( run == 4718 ){
	    du3 -=  0.007 *skwrowA;
	    du3 -=  -0.043*skwrowB;
	    du3 -=  0.013*skwrowC;
	  }
	  if( run == 4719 ){
	    du3 -=  -0.001*skwrowA;
	    du3 -=  0.005*skwrowB;
	    du3 -=  0.004*skwrowC;
	  }
	  if( run == 4720 ) {
	    du3 -=  0.010*skwrowA; // dx3vsskwA->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.004*skwrowB; // dx3vsskwB->Fit("pol1","","",-0.1,0.1)
	    du3 -=  0.015*skwrowC;
	  }
	  if( run == 4765 ){
	    du3 -=  0.026*skwrowA;
	    du3 -=  -0.039*skwrowB;
	    du3 -=  0.018*skwrowC;
	  }

	  if( run == 4767 ){
	    du3 -=  0.017*skwrowA;
	    du3 -=  -0.024*skwrowB;
	    du3 -=  0.018*skwrowC;
	  }

	  if( run == 4768 ){
	    du3 -=  0.007*skwrowA;
	    du3 -=  -0.006*skwrowB;
	    du3 -=  0.0009*skwrowC;
	  }

	  if( run == 4769 ){
	    du3 -=  -6e-6*skwrowA;
	    du3 -=  0.006*skwrowB;
	    du3 -=  -0.0001*skwrowC;
	  }

	  if( run == 4770 ){
	    du3 -=  -0.004*skwrowA;
	    du3 -=  0.016*skwrowB;
	    du3 -=  -0.010*skwrowC;
	  }

	  if( run == 4771 ){
	    du3 -=  0.006*skwrowA;
	    du3 -=  -5e-6*skwrowB;
	    du3 -=  0.001*skwrowC;
	  }

	  if( run == 4772 ){
	    du3 -=  -0.009*skwrowA;
	    du3 -=  0.023*skwrowB;
	    du3 -=  -0.010*skwrowC;
	  }

	  if( run == 4773 ){
	    du3 -=  -0.009*skwrowA;
	    du3 -=  0.028*skwrowB;
	    du3 -=  -0.012*skwrowC;
	  }

	  if( run == 4774 ){
	    du3 -=  -0.009*skwrowA;
	    du3 -=  0.026*skwrowB;
	    du3 -=  -0.009*skwrowC;
	  }

	  if( run == 4775 ){
	    du3 -=  -0.009*skwrowA;
	    du3 -=  0.015*skwrowB;
	    du3 -=  0.003*skwrowC;
	  }

	  if( run == 4776 ){
	    du3 -=  -0.001*skwrowA;
	    du3 -=  -0.010*skwrowB;
	    du3 -=  0.019*skwrowC;        //wavy dx3vsskw distribution                                                            
	  }
	  if( run == 4777 ){
	    du3 -=  0.023*skwrowA;
	    du3 -=  -0.050*skwrowB;
	    du3 -=  0.025*skwrowC;   //wavy dx3vsskw distribution                                                                 
	  }

	  if( run == 4778 ){
	    du3 -=  0.040*skwrowA;
	    du3 -=  -0.081*skwrowB;
	    du3 -=  0.043*skwrowC;
	  }

	  if( run == 4779 ){
	    du3 -=  0.048*skwrowA;
	    du3 -=  -0.081*skwrowB;
	    du3 -=  0.035*skwrowC;   //wavy dx3vsskw distribution                                                                
	  }
	  if( run == 4780 ){
	    du3 -=  0.037*skwrowA;
	    du3 -=  -0.050*skwrowB;
	    du3 -=  0.029*skwrowC;
	  }

	  if( run == 4781 ){
	    du3 -=  0.024*skwrowA;
	    du3 -=  -0.033*skwrowB;
	    du3 -=  0.023*skwrowC;
	  }

	  if( run == 4782 ){
	    du3 -=  0.019*skwrowA;
	    du3 -=  -0.021*skwrowB;
	    du3 -=  0.015*skwrowC;
	  }
	  if( run == 4783 ){
	    du3 -=  0.039*skwrowA;
	    du3 -=  -0.056*skwrowB;
	    du3 -=  0.028*skwrowC;
	  }
	
	  if( run >=4784 && run <= 4792 ){
	    du3 -=  0.005*skwrowA;
	    du3 -=  -0.0006*skwrowB;
	    du3 -=  0.015*skwrowC;
	  }

	  double dv3 = du3;

	  if( run >= 4655 && run <= 4672 ) { // du3vsxm->Fit("pol9")
	    double p0 =    -0.325071;
	    double p1 =    -0.795303;
	    double p2 =      0.64193;
	    double p3 =    -0.166189;
	    double p4 =    0.0194006;
	    double p5 =  -0.00120772;
	    double p6 =  4.29696e-05;
	    double p7 =  -8.7722e-07;
	    double p8 =    9.566e-09;
	    double p9 = -4.32013e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4673){
	    double p0                        =   -0.0126663;
	    double p1                        =     0.757352;
	    double p2                        =    -0.473103;
	    double p3                        =    0.0885939;
	    double p4                        =  -0.00783263;
	    double p5                        =  0.000377553;
	    double p6                        = -1.03967e-05;
	    double p7                        =  1.61115e-07;
	    double p8                        = -1.27933e-09;
	    double p9                        =  3.88606e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4674){
	    double p0                        =    -0.356872;
	    double p1                        =      1.12405;
	    double p2                        =    -0.636105;
	    double p3                        =     0.126649;
	    double p4                        =   -0.0122946;
	    double p5                        =  0.000662474;
	    double p6                        =  -2.0818e-05;
	    double p7                        =  3.79404e-07;
	    double p8                        = -3.71532e-09;
	    double p9                        =  1.51158e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4675){
	    double p0                        =    -0.631455;
	    double p1                        =      1.13891;
	    double p2                        =    -0.625999;
	    double p3                        =     0.126722;
	    double p4                        =    -0.012737;
	    double p5                        =  0.000719102;
	    double p6                        =  -2.3881e-05;
	    double p7                        =  4.63005e-07;
	    double p8                        =  -4.8498e-09;
	    double p9                        =  2.12073e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4676){
	    double p0                        =     -0.10295;
	    double p1                        =    -0.918908;
	    double p2                        =     0.309195;
	    double p3                        =   -0.0385604;
	    double p4                        =   0.00207624;
	    double p5                        = -3.39482e-05;
	    double p6                        = -1.21931e-06;
	    double p7                        =  6.29723e-08;
	    double p8                        = -1.02005e-09;
	    double p9                        =  5.84991e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4677){
	    double p0                        =      1.04075;
	    double p1                        =      -1.2229;
	    double p2                        =     0.429362;
	    double p3                        =    -0.065433;
	    double p4                        =   0.00504498;
	    double p5                        =  -0.00021153;
	    double p6                        =  4.86515e-06;
	    double p7                        =  -5.6606e-08;
	    double p8                        =  2.34351e-10;
	    double p9                        =  4.04142e-13;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4778){
	    double p0                        =     0.156048;
	    double p1                        =    -0.556646;
	    double p2                        =     0.416335;
	    double p3                        =    -0.103542;
	    double p4                        =    0.0117085;
	    double p5                        = -0.000709323;
	    double p6                        =  2.46458e-05;
	    double p7                        = -4.92721e-07;
	    double p8                        =   5.2737e-09;
	    double p9                        =   -2.342e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4679){
	    double p0                        =     0.174788;
	    double p1                        =   0.00355727;
	    double p2                        =     0.248806;
	    double p3                        =   -0.0955524;
	    double p4                        =    0.0130756;
	    double p5                        = -0.000888917;
	    double p6                        =  3.34903e-05;
	    double p7                        = -7.12126e-07;
	    double p8                        =  8.00707e-09;
	    double p9                        = -3.70325e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4680){
	    double p0                        =    -0.347315;
	    double p1                        =     0.115848;
	    double p2                        =     0.144033;
	    double p3                        =    -0.073581;
	    double p4                        =    0.0111161;
	    double p5                        =  -0.00079976;
	    double p6                        =  3.13345e-05;
	    double p7                        = -6.86301e-07;
	    double p8                        =  7.89994e-09;
	    double p9                        = -3.72431e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4681){
	    double p0                        =     -1.04867;
	    double p1                        =    -0.612022;
	    double p2                        =     0.140041;
	    double p3                        =    -0.031657;
	    double p4                        =   0.00449914;
	    double p5                        = -0.000336704;
	    double p6                        =  1.38066e-05;
	    double p7                        = -3.14267e-07;
	    double p8                        =  3.73068e-09;
	    double p9                        = -1.80222e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4682){
	    double p0                        =     -2.16502;
	    double p1                        =     0.245474;
	    double p2                        =    -0.361046;
	    double p3                        =    0.0793908;
	    double p4                        =  -0.00730342;
	    double p5                        =  0.000350287;
	    double p6                        = -9.33598e-06;
	    double p7                        =  1.35909e-07;
	    double p8                        = -9.62812e-10;
	    double p9                        =  2.28003e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4683){
	    double p0                        =     -1.30657;
	    double p1                        =      1.89642;
	    double p2                        =     -1.09961;
	    double p3                        =     0.206532;
	    double p4                        =   -0.0186466;
	    double p5                        =  0.000933453;
	    double p6                        = -2.72754e-05;
	    double p7                        =  4.62042e-07;
	    double p8                        = -4.19551e-09;
	    double p9                        =  1.57549e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4684){
	    double p0                        =     -1.56654;
	    double p1                        =     0.909558;
	    double p2                        =    -0.647281;
	    double p3                        =     0.138739;
	    double p4                        =   -0.0139438;
	    double p5                        =  0.000770002;
	    double p6                        = -2.47417e-05;
	    double p7                        =  4.61016e-07;
	    double p8                        = -4.61788e-09;
	    double p9                        =  1.92312e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4685){
	    double p0                        =      -0.7415;
	    double p1                        =     0.894868;
	    double p2                        =    -0.580423;
	    double p3                        =     0.127798;
	    double p4                        =   -0.0135964;
	    double p5                        =  0.000801617;
	    double p6                        = -2.75508e-05;
	    double p7                        =  5.48934e-07;
	    double p8                        = -5.87485e-09;
	    double p9                        =   2.6121e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4686){
	    double p0                        =    -0.454652;
	    double p1                        =    -0.355914;
	    double p2                        =     0.028462;
	    double p3                        =    0.0159597;
	    double p4                        =  -0.00328118;
	    double p5                        =  0.000262601;
	    double p6                        = -1.08505e-05;
	    double p7                        =  2.44716e-07;
	    double p8                        = -2.86236e-09;
	    double p9                        =  1.35999e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4687){
	    double p0                        =     -1.24549;
	    double p1                        =     0.450725;
	    double p2                        =    -0.140646;
	    double p3                        =    0.0302067;
	    double p4                        =  -0.00386341;
	    double p5                        =  0.000277075;
	    double p6                        =  -1.1246e-05;
	    double p7                        =  2.56715e-07;
	    double p8                        = -3.07198e-09;
	    double p9                        =   1.4989e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4688){
	    double p0                        =     -1.08743;
	    double p1                        =    -0.243086;
	    double p2                        =     0.221091;
	    double p3                        =   -0.0460998 ;
	    double p4                        =   0.00423092 ;
	    double p5                        =  -0.00020439;
	    double p6                        =  5.51319e-06;
	    double p7                        = -8.17503e-08;
	    double p8                        =  5.97563e-10;
	    double p9                        = -1.51914e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4689){
	    double p0                        =    -0.396292;
	    double p1                        =     -0.10224;
	    double p2                        =     0.104844;
	    double p3                        =   -0.0371985;
	    double p4                        =   0.00515829;
	    double p5                        = -0.000357753;
	    double p6                        =  1.37126e-05;
	    double p7                        = -2.95394e-07;
	    double p8                        =  3.35104e-09;
	    double p9                        = -1.55812e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4690){
	    double p0                        =    -0.416804;
	    double p1                        =     0.808811;
	    double p2                        =    -0.383871;
	    double p3                        =    0.0550709;
	    double p4                        =  -0.00331081;
	    double p5                        =  7.12047e-05;
	    double p6                        =  1.14149e-06;
	    double p7                        =  -8.4119e-08;
	    double p8                        =  1.47364e-09;
	    double p9                        =  -8.8095e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4691){
	    double p0                        =    -0.723034;
	    double p1                        =      1.22617;
	    double p2                        =    -0.615611;
	    double p3                        =     0.109059;
	    double p4                        =   -0.0093066;
	    double p5                        =    0.0004327;
	    double p6                        = -1.14213e-05;
	    double p7                        =  1.67709e-07;
	    double p8                        = -1.23312e-09;
	    double p9                        =  3.27571e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4692){
	    double p0                        =    -0.500521;
	    double p1                        =     0.867962;
	    double p2                        =    -0.606846;
	    double p3                        =     0.130378;
	    double p4                        =   -0.0132175;
	    double p5                        =  0.000735612;
	    double p6                        = -2.37934e-05;
	    double p7                        =  4.46126e-07;
	    double p8                        = -4.49874e-09;
	    double p9                        =    1.888e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4693){
	    double p0                        =    -0.180781;
	    double p1                        =    -0.148336;
	    double p2                        =    -0.141526;
	    double p3                        =    0.0605023;
	    double p4                        =  -0.00855849;
	    double p5                        =  0.000596795;
	    double p6                        = -2.29822e-05;
	    double p7                        =   4.9796e-07;
	    double p8                        = -5.68903e-09;
	    double p9                        =  2.66672e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4694){
	    double p0                        =     0.229323;
	    double p1                        =     -1.47957;
	    double p2                        =     0.708519;
	    double p3                        =    -0.121258;
	    double p4                        =    0.0101635;
	    double p5                        = -0.000471583;
	    double p6                        =  1.26332e-05;
	    double p7                        = -1.92372e-07;
	    double p8                        =  1.51965e-09;
	    double p9                        = -4.68813e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4695){
	    double p0                        =     0.516534;
	    double p1                        =     -1.35567;
	    double p2                        =     0.702961;
	    double p3                        =    -0.130991;
	    double p4                        =    0.0118132;
	    double p5                        = -0.000586666;
	    double p6                        =  1.68534e-05;
	    double p7                        = -2.77904e-07;
	    double p8                        =  2.42725e-09;
	    double p9                        = -8.62774e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4698 ) { // du3vsxm->Fit("pol9")
	    double p0 =   -0.0907373;
	    double p1 =    0.0671383;
	    double p2 =    -0.215689;
	    double p3 =    0.0611098;
	    double p4 =  -0.00736417;
	    double p5 =  0.000469881;
	    double p6 = -1.70595e-05;
	    double p7 =  3.53727e-07;
	    double p8 = -3.90056e-09;
	    double p9 =  1.77445e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4701){
	    double p0                        =    -0.728637;
	    double p1                        =     0.571186;
	    double p2                        =    -0.277172;
	    double p3                        =    0.0447151;
	    double p4                        =  -0.00355599;
	    double p5                        =  0.000169667;
	    double p6                        = -5.18006e-06;
	    double p7                        =  9.89537e-08;
	    double p8                        = -1.06728e-09;
	    double p9                        =  4.91556e-12; 
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4702 ) { // du3vsxm->Fit("pol9")
	    double p0 =    -0.768999;
	    double p1 =      0.70397;
	    double p2 =    -0.193329;
	    double p3 =     0.021694;
	    double p4 =  -0.00139475;
	    double p5 =  6.59411e-05;
	    double p6 = -2.30282e-06;
	    double p7 =  5.14119e-08;
	    double p8 = -6.24417e-10;
	    double p9 =  3.10022e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4703){
	    double p0                        =      0.50879;
	    double p1                        =    -0.125573;
	    double p2                        =    0.0930462;
	    double p3                        =   -0.0259806;
	    double p4                        =   0.00279214;
	    double p5                        = -0.000147135;
	    double p6                        =  4.25394e-06;
	    double p7                        = -6.94304e-08;
	    double p8                        =   6.0605e-10;
	    double p9                        = -2.22688e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }
  
	  if( run == 4704){
	    double p0                        =      1.53221;
	    double p1                        =     -1.41434; 
	    double p2                        =     0.480809; 
	    double p3                        =   -0.0890905; 
	    double p4                        =   0.00904545; 
	    double p5                        = -0.000521497; 
	    double p6                        =  1.76384e-05; 
	    double p7                        = -3.47725e-07; 
	    double p8                        =  3.70614e-09; 
	    double p9                        = -1.65236e-11; 
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4705 ){
	    double p0                        =     0.717189;
	    double p1                        =     -1.40859;
	    double p2                        =     0.619353;
	    double p3                        =    -0.108077;
	    double p4                        =   0.00946668;
	    double p5                        = -0.000473501;
	    double p6                        =  1.43085e-05;
	    double p7                        = -2.60295e-07;
	    double p8                        =  2.63272e-09;
	    double p9                        = -1.13778e-11;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4706 ){
	    double p0                        =     0.242833;
	    double p1                        =     -1.03249;
	    double p2                        =     0.466978;
	    double p3                        =   -0.0745645;
	    double p4                        =   0.00584715;
	    double p5                        = -0.000266713;
	    double p6                        =  7.70054e-06;
	    double p7                        = -1.41841e-07;
	    double p8                        =  1.52433e-09;
	    double p9                        = -7.19107e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4707){
	    double p0                        =     0.102753;
	    double p1                        =     -1.06125;
	    double p2                        =     0.589128;
	    double p3                        =    -0.100086;
	    double p4                        =   0.00815266;
	    double p5                        = -0.000380571;
	    double p6                        =  1.10035e-05;
	    double p7                        = -1.98072e-07;
	    double p8                        =  2.04488e-09;
	    double p9                        = -9.21549e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4709){
	    double p0                        =     -0.51582;
	    double p1                        =    -0.922962;
	    double p2                        =     0.552138;
	    double p3                        =   -0.0934296;
	    double p4                        =   0.00796253;
	    double p5                        = -0.000405751;
	    double p6                        =  1.30167e-05;
	    double p7                        = -2.57172e-07;
	    double p8                        =  2.83904e-09;
	    double p9                        = -1.33103e-11;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4710 ) { // du3vsxm->Fit("pol9")
	    double p0 =    -0.691367;
	    double p1 =     0.711649;
	    double p2 =    -0.338717;
	    double p3 =    0.0528456;
	    double p4 =  -0.00402358;
	    double p5 =  0.000181132;
	    double p6 = -5.19231e-06;
	    double p7 =  9.39374e-08;
	    double p8 = -9.74706e-10;
	    double p9 =   4.3856e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4711 ){
	    double p0                        =     -1.01592;
	    double p1                        =     0.175333;
	    double p2                        =    -0.232156;
	    double p3                        =    0.0414381;
	    double p4                        =  -0.00293591;
	    double p5                        =  0.000105819;
	    double p6                        = -2.09488e-06;
	    double p7                        =  2.26769e-08;
	    double p8                        = -1.23901e-10;
	    double p9                        =  2.70779e-13;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4712){
	    double p0                        =    -0.796769;
	    double p1                        =    -0.418083;
	    double p2                        =   0.00591865;
	    double p3                        =    0.0034689;
	    double p4                        =  0.000249506;
	    double p5                        = -5.05092e-05;
	    double p6                        =  2.52898e-06;
	    double p7                        = -5.81109e-08;
	    double p8                        =  6.41288e-10;
	    double p9                        = -2.74905e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4713){
	    double p0                        =    -0.562817;
	    double p1                        =     -0.94152;
	    double p2                        =     0.226679;
	    double p3                        =   -0.0305434;
	    double p4                        =   0.00302373;
	    double p5                        = -0.000191015;
	    double p6                        =  7.11834e-06;
	    double p7                        = -1.51098e-07;
	    double p8                        =  1.69064e-09;
	    double p9                        = -7.74525e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4714){
	    double p0                        =     -1.10458;
	    double p1                        =   -0.0672636;
	    double p2                        =    0.0990226;
	    double p3                        =   -0.0172004;
	    double p4                        =   0.00153986;
	    double p5                        = -8.68768e-05;
	    double p6                        =  3.14928e-06;
	    double p7                        = -6.94665e-08;
	    double p8                        =  8.35443e-10;
	    double p9                        =  -4.1683e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4715){
	    double p0                        =     -1.99733;
	    double p1                        =      1.54421;
	    double p2                        =    -0.470803;
	    double p3                        =    0.0769461;
	    double p4                        =  -0.00709823;
	    double p5                        =  0.000381292;
	    double p6                        = -1.21303e-05;
	    double p7                        =  2.25078e-07;
	    double p8                        = -2.25179e-09;
	    double p9                        =  9.39392e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4716){
	    double p0                        =    -0.470539;
	    double p1                        =     0.340428;
	    double p2                        =   -0.0562998;
	    double p3                        =     0.013201;
	    double p4                        =  -0.00201669;
	    double p5                        =  0.000150477;
	    double p6                        = -5.92187e-06;
	    double p7                        =  1.27419e-07;
	    double p8                        = -1.42519e-09;
	    double p9                        =  6.50682e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4717){
	    double p0                        =    -0.270259;
	    double p1                        =     0.552108;
	    double p2                        =   -0.0740694;
	    double p3                        =     0.012479;
	    double p4                        =  -0.00196827;
	    double p5                        =   0.00015358;
	    double p6                        = -6.17729e-06;
	    double p7                        =   1.3382e-07;
	    double p8                        = -1.49324e-09;
	    double p9                        =  6.76315e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }
	  if( run == 4718){
	    double p0                        =      1.28403;
	    double p1                        =      -1.5997;
	    double p2                        =     0.724072;
	    double p3                        =    -0.121527;
	    double p4                        =    0.0101559;
	    double p5                        = -0.000487808;
	    double p6                        =  1.42623e-05;
	    double p7                        = -2.52558e-07;
	    double p8                        =  2.49853e-09;
	    double p9                        = -1.06043e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4719){
	    double p0                        =    0.0398055;
	    double p1                        =    -0.103044;
	    double p2                        =  -0.00194322;
	    double p3                        =  -0.00897597;
	    double p4                        =   0.00174924;
	    double p5                        = -0.000122251;
	    double p6                        =  4.20993e-06;
	    double p7                        = -7.70912e-08;
	    double p8                        =  7.18591e-10;
	    double p9                        = -2.67252e-12;
	    double xm = xmod*1e3;
            double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
            dv3 -= bx*1e-3;
	  }

	  if( run == 4720 ) { // du3vsxm->Fit("pol9")
	    double p0 =    -0.101448;
	    double p1 =    0.0293417;
	    double p2 =    0.0349432;
	    double p3 =   -0.0162135;
	    double p4 =   0.00209967;
	    double p5 = -0.000123816;
	    double p6 =  3.89767e-06;
	    double p7 = -6.82496e-08;
	    double p8 =  6.30794e-10;
	    double p9 = -2.41363e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4765){
	    double p0                        =    0.0792245;
	    double p1                        =    -0.458272;
	    double p2                        =     0.250128;
	    double p3                        =   -0.0604359;
	    double p4                        =   0.00721027;
	    double p5                        = -0.000456248;
	    double p6                        =  1.61388e-05;
	    double p7                        = -3.21207e-07;
	    double p8                        =  3.36664e-09;
	    double p9                        = -1.44738e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4767){
	    double p0                        =     0.781998;
	    double p1                        =      -1.1802;
	    double p2                        =     0.451922;
	    double p3                        =   -0.0866637;
	    double p4                        =   0.00913561;
	    double p5                        = -0.000550171;
	    double p6                        =  1.93283e-05;
	    double p7                        = -3.91356e-07;
	    double p8                        =  4.23167e-09;
	    double p9                        = -1.89265e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4768){
	    double p0                        =    -0.365445;
	    double p1                        =     0.631095;
	    double p2                        =    -0.269187;
	    double p3                        =    0.0439505;
	    double p4                        =  -0.00356602;
	    double p5                        =  0.000160385;
	    double p6                        =  -4.1682e-06;
	    double p7                        =  6.17629e-08;
	    double p8                        = -4.77277e-10;
	    double p9                        =  1.45138e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4769){
	    double p0                        =   -0.0314529;
	    double p1                        =     0.412959;
	    double p2                        =    -0.209376;
	    double p3                        =    0.0404968;
	    double p4                        =  -0.00392316;
	    double p5                        =  0.000210362;
	    double p6                        = -6.51966e-06;
	    double p7                        =  1.16171e-07;
	    double p8                        = -1.10382e-09;
	    double p9                        =  4.32698e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4770){
	    double p0                        =     0.804503;
	    double p1                        =    0.0475131;
	    double p2                        =    -0.174121;
	    double p3                        =    0.0448898;
	    double p4                        =  -0.00510426;
	    double p5                        =  0.000309293;
	    double p6                        = -1.06707e-05;
	    double p7                        =  2.10546e-07;
	    double p8                        = -2.21597e-09;
	    double p9                        =  9.66209e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4771){
	    double p0                        =    -0.279496;
	    double p1                        =     0.548389;
	    double p2                        =    -0.238064;
	    double p3                        =     0.041604;
	    double p4                        =  -0.00369872;
	    double p5                        =  0.000183844;
	    double p6                        = -5.31917e-06;
	    double p7                        =  8.89039e-08;
	    double p8                        = -7.94716e-10;
	    double p9                        =  2.93523e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4772){
	    double p0                        =     0.442192;
	    double p1                        =    -0.753329;
	    double p2                        =     0.272866;
	    double p3                        =   -0.0410454;
	    double p4                        =   0.00309286;
	    double p5                        = -0.000134627;
	    double p6                        =  3.66253e-06;
	    double p7                        = -6.29636e-08;
	    double p8                        =  6.31415e-10;
	    double p9                        = -2.80328e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4773){
	    double p0                        =    -0.120288;
	    double p1                        =   0.00419937;
	    double p2                        =   -0.0204995;
	    double p3                        =   0.00488272;
	    double p4                        = -0.000749913;
	    double p5                        =  5.58433e-05;
	    double p6                        = -2.09333e-06;
	    double p7                        =  4.09968e-08;
	    double p8                        = -3.98875e-10;
	    double p9                        =   1.5038e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4774){
	    double p0                        =     0.614147;
	    double p1                        =    -0.114214;
	    double p2                        =     -0.11514;
	    double p3                        =    0.0280403;
	    double p4                        =  -0.00324312;
	    double p5                        =  0.000208186;
	    double p6                        = -7.62135e-06;
	    double p7                        =  1.58241e-07;
	    double p8                        = -1.73745e-09;
	    double p9                        =  7.84733e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4775){
	    double p0                        =    -0.324483;
	    double p1                        =    -0.178035;
	    double p2                        =    0.0778996;
	    double p3                        =   -0.0269085;
	    double p4                        =   0.00340063;
	    double p5                        = -0.000214935;
	    double p6                        =  7.67111e-06;
	    double p7                        = -1.57133e-07;
	    double p8                        =  1.72195e-09;
	    double p9                        = -7.82001e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4776){
	    double p0                        =    -0.259959;
	    double p1                        =    -0.380739;
	    double p2                        =    0.0772621;
	    double p3                        =   -0.0278648;
	    double p4                        =   0.00393767;
	    double p5                        = -0.000265479;
	    double p6                        =  9.77815e-06;
	    double p7                        =  -2.0253e-07;
	    double p8                        =  2.21806e-09;
	    double p9                        = -1.00041e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4777){
	    double p0                        =    0.0395608;
	    double p1                        =    -0.252182;
	    double p2                        =    -0.207962;
	    double p3                        =    0.0474461;
	    double p4                        =  -0.00439984;
	    double p5                        =  0.000226975;
	    double p6                        = -7.04411e-06;
	    double p7                        =  1.30973e-07;
	    double p8                        = -1.34418e-09;
	    double p9                        =  5.84693e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }
	  if( run == 4778){
	    double p0                        =      -0.2559;
	    double p1                        =    -0.386069;
	    double p2                        =    -0.083846;
	    double p3                        =    0.0309716;
	    double p4                        =  -0.00336446;
	    double p5                        =  0.000190922;
	    double p6                        =  -6.3923e-06;
	    double p7                        =  1.27273e-07;
	    double p8                        = -1.39155e-09;
	    double p9                        =  6.41442e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4779){
	    double p0                        =    -0.110582;
	    double p1                        =    -0.240694;
	    double p2                        =    0.0452305;
	    double p3                        =  -0.00445865;
	    double p4                        =  0.000591002;
	    double p5                        = -5.07034e-05;
	    double p6                        =  2.26568e-06;
	    double p7                        = -5.38744e-08;
	    double p8                        =  6.54275e-10;
	    double p9                        = -3.20685e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4780){
	    double p0                        =     -1.10945;
	    double p1                        =     0.822837;
	    double p2                        =    -0.250605;
	    double p3                        =    0.0489609;
	    double p4                        =  -0.00510379;
	    double p5                        =  0.000297547;
	    double p6                        = -1.01431e-05;
	    double p7                        =  2.01286e-07;
	    double p8                        = -2.15553e-09;
	    double p9                        =  9.63223e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4781){
	    double p0                        =     -1.25085;
	    double p1                        =     0.751764;
	    double p2                        =   -0.0735535;
	    double p3                        =   0.00169272;
	    double p4                        =  0.000341035;
	    double p5                        = -4.00291e-05;
	    double p6                        =  1.95178e-06;
	    double p7                        = -4.92364e-08;
	    double p8                        =  6.28208e-10;
	    double p9                        = -3.20353e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4782){
	    double p0                        =     -1.26679;
	    double p1                        =      2.02177;
	    double p2                        =    -0.617859;
	    double p3                        =     0.100116;
	    double p4                        =  -0.00892219;
	    double p5                        =  0.000461451;
	    double p6                        = -1.42826e-05;
	    double p7                        =  2.61191e-07;
	    double p8                        = -2.60441e-09;
	    double p9                        =  1.09238e-11;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run == 4783){
	    double p0                        =     -1.67156;
	    double p1                        =      2.19669;
	    double p2                        =    -0.679859;
	    double p3                        =     0.105404;
	    double p4                        =  -0.00865547;
	    double p5                        =  0.000406127;
	    double p6                        = -1.13221e-05;
	    double p7                        =  1.85343e-07;
	    double p8                        = -1.63955e-09;
	    double p9                        =   6.0164e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }

	  if( run >= 4784 && run <=4792){
	    double p0                        =    -0.244237;
	    double p1                        =    0.0172956;
	    double p2                        =   -0.0250505;
	    double p3                        =   0.00512527;
	    double p4                        = -0.000361566;
	    double p5                        =  4.92431e-06;
	    double p6                        =  4.80918e-07;
	    double p7                        = -2.29317e-08;
	    double p8                        =  3.85625e-10;
	    double p9                        = -2.31482e-12;
	    double xm = xmod*1e3;
	    double bx = p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm ) * xm;
	    dv3 -= bx*1e-3;
	  }


	  double dy3 = yB - yavg;

	  hdx3.Fill( dx3*1E3 );
	  hdy3.Fill( dy3*1E3 );

	  if( fabs( dy3 ) < 3*40*1E-3 * 5/pp + 0.05 ) { // cut on y, look at x, see madxvsy

	    hdx3c.Fill( dx3*1E3 );
	    hdu3c.Fill( du3*1E3 );
	    if( cB->size == 2 )
	      hdu3c2.Fill( du3*1E3 );
	    if( cA->size == 2 && cB->size == 2 && cC->size == 2 )
	      hdu3c222.Fill( du3*1E3 ); // 1020: 2.89

	    htx3.Fill( txCA*1e3 ); // [mrad]

	    if( dx3 > 0.04 && dx3 < -0.06 ) { // side lobe
	      cout << endl;
	      cout << "x: " << xAr << ", " << xB << ", " << xCr << ", dx3 " << dx3 << endl;
	      cout << "A:";
	      for( unsigned icl = 0; icl < vclA.size(); ++icl )
		cout << " (" << vclA[icl].col
		     << ", " << vclA[icl].row
		     << ", " << vclA[icl].q << ")";
	      cout << endl;
	      cout << "B:";
	      for( unsigned icl = 0; icl < vclB.size(); ++icl )
		cout << " (" << vclB[icl].col
		     << ", " << vclB[icl].row
		     << ", " << vclB[icl].q << ")";
	      cout << endl;
	      cout << "C:";
	      for( unsigned icl = 0; icl < vclC.size(); ++icl )
		cout << " (" << vclC[icl].col
		     << ", " << vclC[icl].row
		     << ", " << vclC[icl].q << ")";
	      cout << endl;
	    }

	    if( cB->iso )
	      hdx3ci.Fill( dx3*1E3 );

	    if( cA->iso && cC->iso )
	      hdx3cii.Fill( dx3*1E3 );

	    if( cA->iso && cB->iso && cC->iso ) {

	      hdx3ciii.Fill( dx3*1E3 ); // clean

	      if( cB->size == 1 )
		hdx3c1.Fill( dx3*1E3 ); // r447 4.4
	      if( cB->size == 2 )
		hdx3c2.Fill( dx3*1E3 ); // r447 4.4
	      if( cB->size == 3 )
		hdx3c3.Fill( dx3*1E3 ); // r447 4.8
	      if( cB->size == 4 )
		hdx3c4.Fill( dx3*1E3 ); // r447 6.5
	      if( cB->size == 5 )
		hdx3c5.Fill( dx3*1E3 ); // r447 16.6
	      if( cB->size == 6 )
		hdx3c6.Fill( dx3*1E3 ); // r447 24.5
	      if( cB->size > 6 )
		hdx3c7.Fill( dx3*1E3 ); // r447 39.9

	      if( xB < 0 )
		hdx3m.Fill( dx3*1E3 );
	      else
		hdx3p.Fill( dx3*1E3 );

	      if( fabs( dxCA ) < 40*3E-3 * 5/pp ) { // track angle
		hdx3ct.Fill( dx3*1E3 );
		madx3vsq.Fill( cB->q, fabs(dx3)*1E3 );
		madu3vsq.Fill( cB->q, fabs(du3)*1E3 );
		madx3vsn.Fill( cB->size, fabs(dx3)*1E3 );
	      }

	    } // iso

	    if( cB->q > qLB && cB->q < qRB ) {

	      hdx3cq.Fill( dx3*1E3 );

	      if( cA->iso && cB->iso && cC->iso )
		hdx3cqi.Fill( dx3*1E3 );


	      if( cA->q > qL && cA->q < qR &&
		  cC->q > qL && cC->q < qR ) {

		if(cB->size == 1)
		    hdx3clsz1.Fill(dx3*1E3);
                if(cB->size == 2)
                    hdx3clsz2.Fill(dx3*1E3);
                if(cB->size > 2)
                    hdx3clsz3.Fill(dx3*1E3);

							  
		hclszA3cut.Fill( cA->size );
		hnrowA3cut.Fill(nrowA);
		hclszB3cut.Fill( cB->size );
		hnrowB3cut.Fill(nrowB);
		hclszC3cut.Fill( cC->size );
		hnrowC3cut.Fill(nrowC);
		// rms = dx3*dx3; 
		// noe += noe; 
		hdx3cq3.Fill( dx3*1E3 ); // 1020: 2.85
		hdu3cq3.Fill( du3*1E3 ); // skw corrected 1020: 2.74
		hdv3cq3.Fill( dv3*1E3 ); // bias corrected 4702 2.495
		hdxfl3cq3.Fill( dxfl3*1E3 ); // first-last 1020: 2.845
		if( xmod2 < 0.0125 )
		  hdu3cq3l.Fill( du3*1E3 ); // 1020: 2.43
		else
		  hdu3cq3r.Fill( du3*1E3 ); // 1020: 3.09

		if( cA->iso && cB->iso && cC->iso )
		  hdx3cq3i.Fill( dx3*1E3 );

		dx3vsev.Fill( iev, dx3*1E3 );

		dx3vsx.Fill( xB, dx3*1E3 ); // turn
		dx3vsy.Fill( yB, dx3*1E3 ); // rot
		dx3vsdx.Fill( dxCA*1E3, dx3*1E3 ); // slope -0.015 mu/mu

		for( vector<cluster>::iterator cB = vclB.begin(); cB != vclB.end(); ++cB ) {

		  double xB = cB->row * ptchr - 4.0;
		  double yB = cB->col * ptchc - 3.9;

		  double xBmod = fabs(fmod( xB + 8.0125, 0.1 ));
		  double xBp = fmod( xB + 8.0125, ptchr);
		  if(cB->q > qL && cB->q < qR
		     && cA->q > qL && cA->q < qR
		     && cC->q > qL && cC->q < qR)
		    {
		      hxbmod.Fill( xBmod*1E3 );
		      hncolB3cut.Fill(ncolB);
		      hxBmodptch.Fill(xBp*1E3);

		      for(int i = 1; i < 8; i++)
			{
			  if( nrowB == i && ncolB == 1 )
			    {
			      hxbmod_clsz.Fill( i, xBmod*1E3 );
			    }
			}

		      if( nrowB > 1 )
			{
			  double xBcorr = xBp + ( xbBias->GetBinContent(xbBias->FindBin( xBp*1E3 )) ); 
			  double xBmodcorr = fmod( xBcorr + 8.0125, ptchr );
			  hxcorrB.Fill( xBmodcorr*1E3 );
			}

		    }
		}
	
		hxmod.Fill( xmod*1E3 ); // (A+C)/2
		dx3vsxm.Fill( xmod*1E3, dx3*1E3 );
		du3vsxm.Fill( xmod*1E3, du3*1E3 );
		dv3vsxm.Fill( xmod*1E3, dv3*1E3 );

		dx3vsxm2.Fill( xmod2*1E3, dx3*1E3 );
		du3vsxm2.Fill( xmod2*1E3, du3*1E3 );

		dx3vsskwA.Fill( skwrowA, dx3*1E3 );
		dx3vsskwB.Fill( skwrowB, dx3*1E3 );
		dx3vsskwC.Fill( skwrowC, dx3*1E3 );

		du3vsskwA.Fill( skwrowA, du3*1E3 );
		du3vsskwB.Fill( skwrowB, du3*1E3 );
		du3vsskwC.Fill( skwrowC, du3*1E3 );

		madx3vsdx.Fill( dxCA*1E3, fabs(dx3)*1E3 );

		hxcvsxb.Fill( fmod(xB,1), fmod(xC,1) );// corr

		if( fabs( dxCA ) < 40*3E-3 * 5/pp ) { // track angle

		  hdx3cq3t.Fill( dx3*1E3 ); // 447 4.27 um

		  madx3vsx.Fill( xB, fabs(dx3)*1E3 );
		  madx3vsy.Fill( yB, fabs(dx3)*1E3 );
		  madx3vsxm.Fill( xmod*1E3, fabs(dx3)*1E3 );
		  madu3vsxm.Fill( xmod*1E3, fabs(du3)*1E3 );
		  madv3vsxm.Fill( xmod*1E3, fabs(dv3)*1E3 );
		  madu3vsxm2.Fill( xmod2*1E3, fabs(du3)*1E3 );

		  if( cB->size == 2 ) {
		    etavsxmB3.Fill( xmod*1E3, etaB ); // sine
		    madx3vseta.Fill( etaB, fabs(dx3)*1E3 ); // flat
		    hdx3cq3t2.Fill( dx3*1E3 ); // 447 4.25 um
		  }

		} // angle

	      } // Qa, qC

	    } // qB

	  } // cut dy

	  if( fabs( dx3 ) < 0.07 && // hit on track
	      fabs( dy3 ) < 0.15
	      //&& cA->iso && cB->iso && cC->iso
	      ) {

	    hclmapB3->Fill( cB->col, cB->row );

	    hxA3.Fill( xAr );
	    hyA3.Fill( yAr );
	    hxB3.Fill( xB  );
	    hyB3.Fill( yB  );
	    hxC3.Fill( xCr );
	    hyC3.Fill( yCr );
	
	    if(cA->q > qR)
	      hclszA3a.Fill( cA->size );

	    if(cB->q > qR)
              hclszB3a.Fill( cB->size );

	    if(cB->q < qL)
	      hclszB3b.Fill( cB->size );

	    if(cC->q > qR)
              hclszC3a.Fill( cC->size );

	    hclszA3.Fill( cA->size );
	    hclszB3.Fill( cB->size );
	    hclszC3.Fill( cC->size );
	    hncolB3.Fill( ncolB );
	    hnrowA3.Fill( nrowA );
	    hnrowB3.Fill( nrowB );
	    hnrowC3.Fill( nrowC );
	    hclqvszB.Fill( cB->q, cB->size );
	    hclqvszA.Fill( cA->q, cA->size );
	    hclqvszC.Fill( cC->q, cC->size );

	    hclqAvsB.Fill( cA->q, cB->q );
	    hclqAvsC.Fill( cA->q, cC->q );
	    hclqBvsC.Fill( cB->q, cC->q );

	    hclszAvsB.Fill( cA->size, cB->size );
	    hclszAvsC.Fill( cA->size, cB->size );
	    hclszBvsC.Fill( cB->size, cC->size );

	    hclqA3.Fill( cA->q );
	    hclqB3.Fill( cB->q );
	    hclqB3i.Fill( cB->q );
	    hclqC3.Fill( cC->q );
	    if( cB->size < 4 )
	      hclqB3n.Fill( cB->q );

	    nrowvsxmB3.Fill( xmod*1E3, nrowB );
	    nrowvsymB3.Fill( ymod*1E3, nrowB );
	    ncolvsxmB3.Fill( xmod*1E3, ncolB );
	    ncolvsymB3.Fill( ymod*1E3, ncolB );

	    clqvsxmB3.Fill( xmod*1E3, cB->q );
	    // clqvsxB.Fill(xB, cB->q);
	    // q1stvs2nd.Fill( cB->vpix[0].q , cB->vpix[1].q );

	    if( cA->size == 2 )
	      hetaA3.Fill( etaA );

	    if( cB->size == 2 )
	      hetaB3.Fill( etaB );

	    if( cC->size == 2 )
	      hetaC3.Fill( etaC );

	    for( int ipx = 0; ipx < cA->size; ++ipx ) {
	      hpxqA3.Fill( cA->vpix[ipx].q );
	    }
	    for( int ipx = 0; ipx < cB->size; ++ipx ) {
	      hpxqB3.Fill( cB->vpix[ipx].q );
	    }
	    for( int ipx = 0; ipx < cC->size; ++ipx ) {
	      hpxqC3.Fill( cC->vpix[ipx].q );
	    }

	    if( cB->size == 2 ) {
	      hpxq1stB3.Fill( cB->vpix[0].q );
	      hpxq2ndB3.Fill( cB->vpix[1].q ); // identical
	      //hpx3rdB3.Fill( cB->vpix[2].q );
	    }

	    if( nrowA > 1 ) {
	      hrmsA3.Fill( rmsrowA );
	      hskwA3.Fill( skwrowA );
	    }
	    if( nrowC > 1 ) {
	      hrmsC3.Fill( rmsrowC );
	      hskwC3.Fill( skwrowC );
	    }
	    if( nrowB > 1 ) {
	      hrmsB3.Fill( rmsrowB );
	      hskwB3.Fill( skwrowB );
	      rmsvsxmB3.Fill( xmod*1E3, rmsrowB );
	      skwvsxmB3.Fill( xmod*1E3, skwrowB );
	    }

	    // task: store track

	  } // linked, iso

	  // for eff: nearest pixel
	 
	  for( int ipx = 0; ipx < cB->size; ++ipx ) {

	    double px = cB->vpix[ipx].row * ptchr - 4.0; // rot90
	    double py = cB->vpix[ipx].col * ptchc - 3.9;
	    double pdx = px - xavg;
	    double pdy = py - yavg;
	    double pdxy = sqrt( pdx*pdx + pdy*pdy );
	    if( pdxy < pdmin ) pdmin = pdxy;

	  } // pix

	} // clusters B

	if( fabs( dyCA ) > 40*3E-3/pp*5 ) continue; // clean reference "tracks"

	if( evinfoB->filled == ADD ) continue; // padded event in B
	if( evinfoB->skip ) continue; // fat event in B
	if( abs(ddtAB) > 400 ) continue;
	if( abs(ddtCB) > 400 ) continue;
	if( abs(ddtCA) > 400 ) continue;
	if( cA->iso == 0 ) continue;
	if( cC->iso == 0 ) continue;

	int eff[999] = {0};
	for( int iw = 1; iw < 999; ++iw )
	  if( pdmin < iw*0.010 ) // 10 um bins
	    eff[iw] = 1; // eff

	effvsxy->Fill( xavg, yavg, eff[50] );

	if( yavg > -3.7 && yavg < 3.5 )
	  effvsx.Fill( xavg, eff[50] );

	if( xavg > -3.3 && xavg < 3.6 )
	  effvsy.Fill( yavg, eff[50] );

	if( xavg > -3.3 && xavg < 3.6 &&
	    yavg > -3.7 && yavg < 3.5 ) {

	  for( int iw = 1; iw < 999; ++iw )
	    effvsdxy.Fill( iw*0.010+0.005, eff[iw] );

	  effvsxm.Fill( xmod*1E3, eff[50] ); // bias dot
	  effvsev.Fill( iev, eff[50] );
	  effvsiev.Fill( iev%200, eff[50] );
	  effvsqA.Fill( cA->q, eff[50] );
	  effvstxy.Fill( dxyCA, eff[50] ); // flat

	} // fiducial x, y

      } // clusters A

    } // cl C

    nmvsevCA.Fill( iev, nm );

    // task: track correlations, intersects
  } // events

  cout << endl;

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s3 = ts.tv_sec; // seconds since 1.1.1970
  long f3 = ts.tv_nsec; // nanoseconds
  zeit2 += s3 - s2 + ( f3 - f2 ) * 1e-9; // read and cluster
  cout << "time " << zeit2 << " s" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s9 = ts.tv_sec; // seconds since 1.1.1970
  long f9 = ts.tv_nsec; // nanoseconds

  cout << "full time " << s9 - s0 + ( f9 - f0 ) * 1e-9 << " s"
       << " (read and cluster " << zeit1 << " s, tracking " << zeit2 << " s)"
       << endl;

  rusage usage;
  getrusage( RUSAGE_SELF, &usage );
  cout << endl << "resource usage:"
       << endl << "time " << usage.ru_utime.tv_sec + usage.ru_utime.tv_usec*1E-6 << " s"
       << endl << "max resident set size " << usage.ru_maxrss << " kB"
       << endl << "file inputs " << usage.ru_inblock
       << endl << "file ouputs " << usage.ru_oublock
       << endl << "page faults without I/O " << usage.ru_minflt
       << endl << "page faults with I/O " << usage.ru_majflt
       << endl << "voluntary context switches " << usage.ru_nvcsw
       << endl << "involuntary context switches " << usage.ru_nivcsw
       << endl;

  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment fits:

  double newalignxA = alignxA;

  cout << endl << hdxAB.GetTitle() << " entries " << hdxAB.GetEntries() << endl;

  if( hdxAB.GetEntries() > 999 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
    double xpk = hdxAB.GetBinCenter( hdxAB.GetMaximumBin() );
    fgp0->SetParameter( 0, hdxAB.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, xpk );
    fgp0->SetParameter( 2, hdxAB.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, hdxAB.GetBinContent( hdxAB.FindBin(xpk-1) ) ); // BG
    hdxAB.Fit( "fgp0", "q", "", xpk-1, xpk+1 ); // fit range around peak
    cout << "  Fit Gauss + BG:"
	 << endl << "  area " << fgp0->GetParameter(0)
	 << endl << "  mean " << fgp0->GetParameter(1)
	 << endl << "  sigm " << fgp0->GetParameter(2)
	 << endl << "  offs " << fgp0->GetParameter(3)
	 << endl;
    newalignxA += fgp0->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // y:

  double newalignyA = alignyA;

  cout << endl << hdyAB.GetTitle() << " entries " << hdyAB.GetEntries() << endl;
  if( hdyAB.GetEntries() > 999 ) {
    cout << "  y correction " << hdyAB.GetMean() << endl;
    newalignyA += hdyAB.GetMean();
  }
  else
    cout << "  not enough" << endl;

  // dxvsy -> -rot

  double newalignfA = alignfA;

  cout << endl << dxvsyAB.GetTitle() << " entries " << dxvsyAB.GetEntries() << endl;
  if( aligniteration > 0 && dxvsyAB.GetEntries() > 999 ) {
    dxvsyAB.Fit( "pol1", "q", "", -3, 3 );
    TF1 * fdxvsy = dxvsyAB.GetFunction( "pol1" );
    cout << "  extra rot " << fdxvsy->GetParameter(1) << endl;
    newalignfA += fdxvsy->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // dxvsx -> -turn

  double newaligntA = aligntA;

  cout << endl << dxvsxAB.GetTitle() << " entries " << dxvsxAB.GetEntries() << endl;
  if( aligniteration > 0 && dxvsxAB.GetEntries() > 999 ) {
    dxvsxAB.Fit( "pol1", "q", "", -3, 2 );
    TF1 * fdxvsx = dxvsxAB.GetFunction( "pol1" );
    cout << "  extra turn " << fdxvsx->GetParameter(1) << endl;
    newaligntA += fdxvsx->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // C:

  double newalignxC = alignxC;

  cout << endl << hdxCB.GetTitle() << " entries " << hdxCB.GetEntries() << endl;
  if( hdxCB.GetEntries() > 999 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
    double xpk = hdxCB.GetBinCenter( hdxCB.GetMaximumBin() );
    fgp0->SetParameter( 0, hdxCB.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, xpk );
    fgp0->SetParameter( 2, hdxCB.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, hdxCB.GetBinContent( hdxCB.FindBin(xpk-1) ) ); // BG
    hdxCB.Fit( "fgp0", "q", "", xpk-1, xpk+1 ); // fit range around peak
    cout << "  Fit Gauss + BG:"
	 << endl << "  area " << fgp0->GetParameter(0)
	 << endl << "  mean " << fgp0->GetParameter(1)
	 << endl << "  sigm " << fgp0->GetParameter(2)
	 << endl << "  offs " << fgp0->GetParameter(3)
	 << endl;
    newalignxC += fgp0->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // y:

  double newalignyC = alignyC;

  cout << endl << hdyCB.GetTitle() << " entries " << hdyCB.GetEntries() << endl;
  if( hdyCB.GetEntries() > 999 ) {
    cout << "  y correction " << hdyCB.GetMean() << endl;
    newalignyC += hdyCB.GetMean();
  }
  else
    cout << "  not enough" << endl;

  // dxvsy -> -rot

  double newalignfC = alignfC;

  cout << endl << dxvsyCB.GetTitle() << " entries " << dxvsyCB.GetEntries() << endl;
  if( aligniteration > 0 && dxvsyCB.GetEntries() > 999 ) {
    dxvsyCB.Fit( "pol1", "q", "", -3, 3 );
    TF1 * fdxvsy = dxvsyCB.GetFunction( "pol1" );
    cout << "  extra rot " << fdxvsy->GetParameter(1) << endl;
    newalignfC += fdxvsy->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  // dxvsx -> -turn:

  double newaligntC = aligntC;

  cout << endl << dxvsxCB.GetTitle() << " entries " << dxvsxCB.GetEntries() << endl;
  if( aligniteration > 0 && dxvsxCB.GetEntries() > 999 ) {
    dxvsxCB.Fit( "pol1", "q", "", -3, 2 );
    TF1 * fdxvsx = dxvsxCB.GetFunction( "pol1" );
    cout << "  extra turn " << fdxvsx->GetParameter(1) << endl;
    newaligntC += fdxvsx->GetParameter(1);
  }
  else
    cout << "  not enough" << endl;

  ++aligniteration;
  cout << endl
       << "for " << alignFileName << endl
       << "iteration " << aligniteration << endl
       << "alignxA " << setw(11) << newalignxA << endl
       << "alignyA " << setw(11) << newalignyA << endl
       << "alignfA " << setw(11) << newalignfA << endl
       << "aligntA " << setw(11) << newaligntA << endl
       << "alignxC " << setw(11) << newalignxC << endl
       << "alignyC " << setw(11) << newalignyC << endl
       << "alignfC " << setw(11) << newalignfC << endl
       << "aligntC " << setw(11) << newaligntC << endl
    ;

  cout << "update alignment file? (y/n)" << endl;
  string ans{"y"};
  cin >> ans;
  string YES{"y"};
  if( ans == YES ) {

    ofstream alignFile( alignFileName );

    alignFile << "# alignment for run " << run << endl;
    alignFile << "iteration " << aligniteration << endl;
    alignFile << "alignxA " << setw(11) << newalignxA << endl;
    alignFile << "alignyA " << setw(11) << newalignyA << endl;
    alignFile << "alignfA " << setw(11) << newalignfA << endl;
    alignFile << "aligntA " << setw(11) << newaligntA << endl;
    alignFile << "alignxC " << setw(11) << newalignxC << endl;
    alignFile << "alignyC " << setw(11) << newalignyC << endl;
    alignFile << "alignfC " << setw(11) << newalignfC << endl;
    alignFile << "aligntC " << setw(11) << newaligntC << endl;

    alignFile.close();
  }

  cout << endl << "events " << nev
       << endl << "A-B out of sync " << dsyncAB
       << endl << "C-B out of sync " << dsyncCB
       << endl << "C-A out of sync " << dsyncCA
       << endl;

  cout << endl << histoFile->GetName() << endl;

  // double dx3rms = sqrt(rms)/sqrt(noe);

  cout << endl;

  return 0;
}
