/* -------------------------------------------------------------
 *
 *  file:        command.cpp
 *
 *  description: command line interpreter for Chip/Wafer tester
 *
 *  author:      Beat Meier
 *  modified:    24.7.2017
 *
 * -------------------------------------------------------------
 */

#include "cmd.h"
//#include "datalink.h"

#include <iostream> // cout
#include <string>
#include <fstream> // files
#include <time.h> // clock_gettime

#include <omp.h>// openMP
#include <chrono>
#include <thread>

#include <TFile.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>
#include <TGraph.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>

using namespace std;

#define IMG_WIDTH 157
#define R4S_VALUE_OR 4096

bool line_wise = 1;
unsigned IMG_HEIGHT = 163; // FW 0.6

//------------------------------------------------------------------------------
CMD_PROC(chip)
{
  int chip;
  PAR_INT( chip, 0, 9999 );

  roc.chip = chip;
  roc.haveGain = 0;

  string gainFileName;

  if( chip == 13 )
    gainFileName = "r4s_cal_roc013.dat";

  if( chip == 102 )
    gainFileName = "r102_tb21_cal.dat";

  if( gainFileName.length() > 0 ) {

    ifstream gainFile( gainFileName.c_str() );

    if( gainFile ) {

      roc.haveGain = 1;
      cout << "gainFile: " << gainFileName << endl;

      while( ! gainFile.eof() ) {
	int icol;
	int irow;
	gainFile >> icol;
	gainFile >> irow;
	gainFile >> roc.p0[icol][irow];
	gainFile >> roc.p1[icol][irow];
	gainFile >> roc.p2[icol][irow];
	gainFile >> roc.p3[icol][irow];
      }

    } // gainFile open

  } // gainFileName

}

//------------------------------------------------------------------------------
// inverse Fermi PH -> Vcal mV
double PHtoVcal( double ph, unsigned col, unsigned row )
{
  if( !roc.haveGain )
    return ph;

  if( col > 155 )
    return ph;

  if( row > 160 )
    return ph;

  // r4scal.C

  double U = ( ph - roc.p3[col][row] ) / roc.p2[col][row];

  if( U >= 1 ) {
    U = 0.9999999; // avoid overflow
  }

  return roc.p0[col][row] - roc.p1[col][row] * log( (1-U)/U ); // inverse Fermi
}

//------------------------------------------------------------------------------
class R4sImg
{
  int * data;
  int evNum;

public:
  R4sImg() : data(0) {
    evNum = 0;
  }

  void Clear() { if( data ) delete[] data; data = 0; }

  bool CreateRaw( const vector<uint16_t> &rawdata ); // fills data

  int Get( int x, int y ) {
    if( line_wise )
      return data[y*IMG_WIDTH + x];
    else
      return data[x*IMG_HEIGHT + y];
  } // one pix

  void Print( unsigned int count, unsigned line ); // one row

  void Save( const string &filename ); // map.txt
};

//------------------------------------------------------------------------------
bool R4sImg::CreateRaw( const vector<uint16_t> &rawdata )
{
  Clear();

  if( rawdata.size() < IMG_WIDTH * IMG_HEIGHT )
    return false;

  ++evNum;

  data = new int[IMG_WIDTH * IMG_HEIGHT];

  for( unsigned pos = 0; pos < IMG_WIDTH * IMG_HEIGHT; ++pos ) {

    int value = rawdata[pos];

    if( value & 0x1000 ) // ADC overrange
      value = R4S_VALUE_OR;
    else if( value & 0x0800 ) // negative
      value -= 0x1000;

    data[pos] = value;
  }

  return true;
};

//------------------------------------------------------------------------------
void R4sImg::Print( unsigned int count, unsigned line_to_print )
{
  if( !data ) {
    printf( "R4sImg::Print: empty data\n" );
    return;
  }

  if( count > IMG_WIDTH * IMG_HEIGHT )
    count = IMG_WIDTH * IMG_HEIGHT;

  cout << "Event " << evNum << endl;

  unsigned i0 = line_to_print*IMG_WIDTH;
  if( line_wise )
    cout << "row " << i0/IMG_WIDTH << " from bottom" << endl;
  else {
    i0 = line_to_print*IMG_HEIGHT;
    cout << "col " << i0/IMG_HEIGHT << " from left" << endl;
  }
  for( unsigned int i = i0; i < i0+count; ++i ) {
    int v = data[i];
    if( v == R4S_VALUE_OR )
      printf( "   or" );
    else
      printf( " %4i", v );
  }
  printf( "\n" );

  tb.Daq_Close();
}

//------------------------------------------------------------------------------
void R4sImg::Save( const string &filename )
{
  if( ! data ) {
    printf( "R4sImg::Save: empty data\n" );
    return;
  }

  //FILE * f = fopen( filename.c_str(), "wt" );
  FILE * f = fopen( filename.c_str(), "at" ); // append text
  //FILE * f = fopen( filename.c_str(), "ab" ); // append binary
  //FILE * f = fopen( filename.c_str(), "wb" ); // append binary

  if( !f ) return;

  TFile * histoFile = NULL;
  if( evNum == 1 )
    histoFile = new TFile( "one.root", "RECREATE" );

  TProfile2D mapxy( "mapxy",
		    "ADC map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  fprintf( f, " %i\n", evNum );

  if( evNum == 1 )
    cout << "img " << IMG_WIDTH << "x" << IMG_HEIGHT << endl;

  for( unsigned y = 0; y < IMG_HEIGHT; ++y ) { // start row 0

    //cout << setw(3) << y << ":";

    for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

      int v = Get( x, y );

      if( v == R4S_VALUE_OR )
	fprintf( f, "   or" );
      else
	fprintf( f, " %4i", v );

      mapxy.Fill( x+0.5, y+0.5, v );

      //cout << setw(5) << v;

    }

    fputs( "\n", f ); // end of one row

    //cout << endl;

  } // rows

  fclose(f);

  if( evNum == 1 ) {
    histoFile->Write();
    histoFile->Close();
    cout << histoFile->GetName() << endl;
  }

}

//------------------------------------------------------------------------------
void ReadImage( R4sImg &map )
{
  tb.Daq_Open(50000);

  // prepare ADC:

  tb.SignalProbeADC( PROBEA_SDATA1, GAIN_1 );

  tb.uDelay(500); // [us] Beat 23/11/2017: must be larger 400

  // take data:

  tb.Daq_Start();

  tb.r4s_Start(); // start sequence

  tb.uDelay(3000);
  if( roc.ext )
    tb.uDelay(64*1000); // 64 ms max

  tb.Daq_Stop();

  // read buffer:

  vector<uint16_t> vdata;
  unsigned int ret = tb.Daq_Read( vdata );

  printf( ", read status %u, n = %d = 157x%d \n", ret, vdata.size(), vdata.size()/157 );

  tb.Daq_Close();

  map.CreateRaw( vdata );
}

//------------------------------------------------------------------------------
//                  getimg 24900   (each pix once)
CMD_PROC(getimg) // getimg  1240 2 (take raw data, pulse 2 pix, write to file)
{
  int nevents;
  if( ! PAR_IS_INT( nevents, 1, 100*24800+100 ) ) // 155*160 = 24'800 + 100 ped
    nevents = 1;

  int npx;
  if( ! PAR_IS_INT( npx, 1, 4 ) )
    npx = 1;

  int line_to_print;
  if( ! PAR_IS_INT( line_to_print, 0, 160 ) )
    line_to_print = 7; // for row-wise

  cout << "nevents " << nevents << endl;

  // run number from file:

  unsigned int run = 1;

  if( nevents > 1 ) {

    ifstream irunFile( "runNumber.dat" );
    irunFile >> run;
    irunFile.close();

    ++run;
    cout << "Run " << run << endl;

    ofstream orunFile( "runNumber.dat" );
    orunFile << run;
    orunFile.close();
  }

  R4sImg map;

  if( nevents > 99 )
    tb.r4s_SetPixCal( 156, 161 ); // for pedestal

  for( int iev = 0; iev < nevents; ++iev ) {

    if( iev > 99 ) { // 100 quiet events

      if(      npx == 4 )
	tb.r4s_Set4PixCal( (iev-100)%155, (iev-100)/160 ); // start pulsing at pixel 0 0

      else if( npx == 2 )
	tb.r4s_Set2PixCal( (iev-100)%155, (iev-100)/160 ); // start pulsing at pixel 0 0

      else
	tb.r4s_SetPixCal( (iev-100)%155, (iev-100)/160 ); // start pulsing at pixel 0 0

    }

    cout << "ev " << iev;

    ReadImage(map);

    if( line_wise ) {
      if( nevents == 1 || ( iev - 100 ) / 155 == line_to_print ) map.Print(155,line_to_print); // row with pulses
    }
    else
      if( nevents == 1 || ( iev - 100 ) / 160 == line_to_print ) map.Print(160,line_to_print); // col with pulses

    map.Save( Form( "raw%06i.txt", run ) );

  } // ev

  printf( "Written to raw%06i.txt\n", run );

} // getimg

//------------------------------------------------------------------------------
CMD_PROC(td) // roi data
{
  int Nev;
  if( ! PAR_IS_INT( Nev, 1, 100*1000*1000 ) )
    Nev = 10200;

  int ntrg;
  if( ! PAR_IS_INT( ntrg, 1, 200 ) )
    ntrg = 200;

  int planeNr;    // colorize planes in different ways
  if( ! PAR_IS_INT( planeNr, 1, 3 ) )
    planeNr = 1;

  int col_to_print;
  if( ! PAR_IS_INT( col_to_print, 0, 156 ) )
    col_to_print = 0;

  int fupd;     // every fupd events the plots are updated
  if( ! PAR_IS_INT( fupd, 10, 5001 ) )
    fupd = 200;

  bool pulse{1}; // add random test pulse

  // run number from file:

  unsigned int run = 1;

  ifstream irunFile( "runNumber.dat" );
  irunFile >> run;
  irunFile.close();

  ++run;
  cout << "Run " << run << endl;

  ofstream orunFile( "runNumber.dat" );
  orunFile << run;
  orunFile.close();

  string fileName = Form( "roi%06i.txt", run );

  ofstream outfile;

  outfile.setf( ios::fixed );
  outfile.precision(1);
  outfile.open( fileName.c_str() ); // task: check integrity

  int nPedAvg = 200;  // number of events to average for pedestal and noise = initial fupd
  //double thr = -24;   // [ADC] threshold for hit finding (negative!) GAIN_1
  double thr = -4; // [noise significance]
  //int roiCol = 2;     // +-ROI in col until 13.6.2018
  int roiCol = 3;     // +-ROI in col since 13.6.2018 run 2623
  int roiRow = 3;     // +-ROI in row
  bool fiftyfifty = 0;// 1 = 50x50 sensor, 0 = 100x25

  //line_wise = 0; // force column-wise readout

  if( line_wise == 0 )
    IMG_HEIGHT = 162; // Ali FW 0.7

  cout << "line_wise " << line_wise << endl;

  // write header for output file:

  outfile << "run " << run
          << endl << "pedestal " << nPedAvg
          << endl << "threshold " << thr
          << endl << "row_wise " << line_wise
          << endl << "50_by_50 " << fiftyfifty
          << endl << "roi_col " << roiCol*2 + 1
          << endl << "roi_row " << roiRow*2 + 1
          << endl << "START"
    ;

  roc.print(1); // ROC settings to Log
  Log.printf( "  IA      %5.1f mA\n", tb.GetIA()*1E3 );
  Log.printf( "  ID      %5.1f mA\n", tb.GetID()*1E3 );

  Log.printf( "[Run %i takedata %i]\n", run, ntrg );
  Log.printf( "  pedestal block %i events\n", nPedAvg );
  Log.printf( "  seed threshold %3.1f ADC\n", thr );
  if( line_wise )
    Log.printf( "  readout mode is row-wise\n" );
  else
    Log.printf( "  readout mode is column-wise\n" );
  Log.printf( "  ROI window +-%i columns\n", roiCol );
  Log.printf( "  ROI window +-%i rows\n", roiRow );

  TRandom *rnd3 = new TRandom3();

  TFile * histoFile = new TFile( Form( "roi%06i.root", run ), "RECREATE" );

  int argc = 0;
  char* argv[1];

  TApplication * app;

  if( gApplication ) { // singleton
    cout << "gApplication exists" << endl;
    app = gApplication;
  }
  else {
    cout << "no gApplication" << endl;
    app = new TApplication( "app", &argc, argv);
  }

  Char_t name[100];
  sprintf( name,"P%i - Pulse Heights;PH [ADC];N pixels", planeNr );
  TH1D hadc( "hadc", name, 4097, -2048.5, 2048.5 );

  sprintf( name,"P%i - Pulse Heights (incl. Pedestal Correction);PH [ADC];N pixels",planeNr);
  TH1D hadcSubPed( "hadcSubPed", name, 4097, -2048.5, 2048.5 );

  sprintf( name,"P%i - Pulse Heights Differences (incl. Pedestal Correction);PH [ADC];N pixels",planeNr);
  TH1D hadcNoCoMo( "hadcNoCoMo", name, 4097, -2048.5, 2048.5 );

  sprintf( name,"P%i - dPH/noise;significance;N pixels", planeNr );
  TH1D hsigni( "significance", name, 1000, -100, 100 );

  sprintf( name,"P%i - Pixel over Threshold; N Pixel; entries",planeNr);
  TH1D pixelOverThr("pixelOverThr", name, 161, -0.5, 160.5 );

  sprintf( name,"P%i - RMS distribution of #DeltaPH; RMS [ADC]; entries",planeNr);
  TH1D hdphrms( "hdphrms", name, 200, 0, 50 );

  // profiles:

  sprintf( name,"P%i - Running Pedestal Map;col;row;<ADC>",planeNr);
  TProfile2D pedxy( "pedxy", name, 155, 0, 155, 160, 0, 160, -2222, 2222 );

  sprintf( name,"P%i - RMS Pulse Height (incl. Pedestal Correction);col;row;<ADC>",planeNr);
  TProfile2D phRMS( "phRMS", name, 155, 0, 155, 160, 0, 160, -2222, 2222 );

  sprintf(name,"P%i - RMS Pulse Height Differences (incl. Pedestal Correction);col;row;<ADC>",planeNr);
  TProfile2D phCRMS( "phCRMS", name, 155, 0, 155, 160, 0, 160, -2222, 2222 );

  sprintf(name,"P%i - Pulse Height of Pixel over Threshold (incl. Pedestal Correction);col;row;<ADC>",planeNr);
  TProfile2D hitAmplAvg( "hitAmplAvg",
                         name, 155, 0, 155, 160, 0, 160, -2222., 2222. );

  // time development:

  TGraph* gPedAvg = new TGraph();
  gPedAvg->SetName( "gPedAvg" );
  sprintf( name, "P%i - Average Pedestal",planeNr );
  gPedAvg->SetTitle(name);
  gPedAvg->GetXaxis()->SetTitle( "Event Number" );
  gPedAvg->GetYaxis()->SetTitle( "<ADC>" );
  gPedAvg->SetPoint( 0, 0, 0 );

  TGraph* gHitAvg = new TGraph();
  gHitAvg->SetName( "gHitAvg" );
  sprintf( name, "P%i - Average Number of hit Pixels", planeNr );
  gHitAvg->SetTitle( name );
  gHitAvg->GetXaxis()->SetTitle( "Event Number" );
  gHitAvg->GetYaxis()->SetTitle( "N pixels" );
  gHitAvg->SetPoint( 0, 0, 1 );

  TGraph* grmsdph= new TGraph();
  grmsdph->SetName( "grmsdph" );
  sprintf( name, "P%i - Average RMS of #DeltaPH", planeNr );
  grmsdph->SetTitle( name );
  grmsdph->GetXaxis()->SetTitle( "Event Number" );
  grmsdph->GetYaxis()->SetTitle( "RMS [ADC]" );
  grmsdph->SetPoint( 0, 0, 11 );

  // hitmap:

  sprintf(name,"P%i - Hitmap;col;row;[entries]",planeNr);
  TH2D hitmap( "hitmap", name, 155, 0, 155, 160, 0, 160 );

  // canvas for live draw

  //TCanvas * h1 = new TCanvas( "h1", "h1", 600, 400 );
  //h1->SetLogy();
  //hadc.Draw();

  //TCanvas * h2 = new TCanvas( "h2", "h2", 600, 400 );
  //h2->SetLogy();
  //hadcSubPed.Draw();

  TCanvas * h3 = new TCanvas( "h3", "h3", 600, 400 );
  h3->SetLogy();
  hadcNoCoMo.Draw();

  TCanvas * h4 = new TCanvas( "h4", "h4", 600, 400 );
  h4->SetLogy();
  hsigni.Draw();

  //TCanvas * h5 = new TCanvas( "h5", "h5", 600, 400 );
  //h5->SetLogy();
  //pixelOverThr.Draw();

  TCanvas * h6 = new TCanvas( "h6", "h6", 600, 400 );
  //h6->SetLogy();
  hdphrms.Draw();
  /*
  TCanvas * p1 = new TCanvas( "p1", "p1", 900, 800 );
  //pedxy.Draw( "colz" );

  TCanvas * p2 = new TCanvas( "p2", "p2", 900, 800 );
  //phRMS.Draw( "colz" );

  TCanvas * p3 = new TCanvas( "p3", "p3", 900, 800 );
  //phCRMS.Draw( "colz" );

  TCanvas * p4 = new TCanvas( "p4", "p4", 900, 800 );
  //hitAmplAvg.Draw("colz");
  */
  TCanvas * g1 = new TCanvas( "g1", "g1", 200, 10, 600, 400 );
  gPedAvg->Draw( "AC*" );

  TCanvas * g2 = new TCanvas( "g2", "g2", 200, 10, 600, 400 );
  gHitAvg->Draw( "AC*" );

  TCanvas * g3 = new TCanvas( "g3", "g3", 200, 10, 600, 400 );
  grmsdph->Draw( "AC*" );

  TCanvas * d1 = new TCanvas( "d1", "d1", 900, 800 ); // square
  d1->SetLogz();
  hitmap.Draw( "colz" );

  if( planeNr == 1 ){ // upstream = red

    //h1->SetFillColor( kRed-10 );
    //h2->SetFillColor( kRed-10 );
    h3->SetFillColor( kRed-10 );
    h4->SetFillColor( kRed-10 );
    //h5->SetFillColor( kRed-10 );
    h6->SetFillColor( kRed-10 );
    /*
    p1->SetFillColor( kRed-10 );
    p2->SetFillColor( kRed-10 );
    p3->SetFillColor( kRed-10 );
    p4->SetFillColor( kRed-10 );
    */
    g1->SetFillColor( kRed-10 );
    g2->SetFillColor( kRed-10 );
    g3->SetFillColor( kRed-10 );
    d1->SetFillColor( kRed-10 );
  }

  if( planeNr == 2 ) { // middle = yellow

    //h1->SetFillColor( kYellow-10 );
    //h2->SetFillColor( kYellow-10 );
    h3->SetFillColor( kYellow-10 );
    h4->SetFillColor( kYellow-10 );
    //h5->SetFillColor( kYellow-10 );
    h6->SetFillColor( kYellow-10 );
    /*
    p1->SetFillColor( kYellow-10 );
    p2->SetFillColor( kYellow-10 );
    p3->SetFillColor( kYellow-10 );
    p4->SetFillColor( kYellow-10 );
    */
    g1->SetFillColor( kYellow-10 );
    g2->SetFillColor( kYellow-10 );
    g3->SetFillColor( kYellow-10 );
    d1->SetFillColor( kYellow-10 );
  }

  if( planeNr == 3 ){ //down stream = green

    //h1->SetFillColor( kGreen-10 );
    //h2->SetFillColor( kGreen-10 );
    h3->SetFillColor( kGreen-10 );
    h4->SetFillColor( kGreen-10 );
    //h5->SetFillColor( kGreen-10 );
    h6->SetFillColor( kGreen-10 );
    /*
    p1->SetFillColor( kGreen-10 );
    p2->SetFillColor( kGreen-10 );
    p3->SetFillColor( kGreen-10 );
    p4->SetFillColor( kGreen-10 );
    */
    g1->SetFillColor( kGreen-10 );
    g2->SetFillColor( kGreen-10 );
    g3->SetFillColor( kGreen-10 );
    d1->SetFillColor( kGreen-10 );
  }

  tb.Daq_Open( 48*1024*1024 ); // [words] max 64*1024*1024

  const uint32_t Blocksize = 48 * 1024 * 1024;

  // prepare ADC:

  tb.SignalProbeADC( PROBEA_SDATA1, GAIN_1 ); // 2017 gain_1

  tb.uDelay(500); // [us] Beat 23/11/2017: must be larger 400

  const int iFWVersion = tb.GetFWVersion();

  int iev = 0;
  int wev = 0; // write event
  double hitAvg = 0;
  double dphrmsAvg  = 0;

  // Define matrices to store e.g. pedestal information on the run:

  double savePH[155][160] = {0};
  double runPed[155][160] = {0};
  double sumPHsq[155][160] = {0};
  double sumdPHsq[155][160] = {0};
  double noise[155][160];
  double notHitPH[155][160] = {0};
  double notHitPHC[155][160] = {0};
  for( int col = 0; col < 155; ++col )
    for( unsigned row = 0; row < 160; ++row )
      noise[col][row] = 10; // initial noise irrad gain_1
  //noise[col][row] = 6; // initial noise fresh gain_1
  //noise[col][row] = 14; // safety irrad gain_2

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s0 = ts.tv_sec; // seconds since 1.1.1970
  long f0 = ts.tv_nsec; // nanoseconds

  double zeit1 = 0;
  double zeit2 = 0;

  //while( !keypressed() ) {
  while( iev < Nev ) {

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s1 = ts.tv_sec; // seconds since 1.1.1970
    long f1 = ts.tv_nsec; // nanoseconds

    tb.Daq_Start(); // ADC -> RAM

    for( int itrg = 0; itrg < ntrg; ++itrg ) {

      // random test pulse:

      if( pulse && iev > nPedAvg ) {
	int col = rnd3->Integer(155);
	int row = rnd3->Integer(160);
	//cout << " pulse " << col << "  " << row << endl; /// OK
	tb.r4s_SetPixCal( col, row ); // pulse
	tb.uDelay(100);
      }

      ++iev; // triggered events
      cout << " ev " << iev;

      tb.r4s_Start(); // R4S sequence

      if( roc.ext == 0 ) //  internal trig
	tb.uDelay(700); // 157*163*25 ns = 640 us R4S readout

    }

    if( roc.ext ) tb.uDelay( 64*1000 ); // uint16 goes to 64*1024-1

    tb.Daq_Stop();

    // read DTB memory:

    vector<uint16_t> vdata;
    vdata.reserve( ( IMG_WIDTH *  IMG_HEIGHT + 4 ) * ntrg );

    unsigned int ret = tb.Daq_Read( vdata, Blocksize );

    cout << "  status " << ret << " read " << vdata.size() << " words";

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s2 = ts.tv_sec; // seconds since 1.1.1970
    long f2 = ts.tv_nsec; // nanoseconds
    zeit1 += s2 - s1 + ( f2 - f1 ) * 1e-9; // read and cluster
    cout << " (in " << s2 - s1 + ( f2 - f1 ) * 1e-9 << " s)";

    cout << endl;

    int timestampblocksize = 0;

    int ktrg = vdata.size() / IMG_WIDTH / IMG_HEIGHT; // coarse

    if( roc.ext && iFWVersion > 256 ) // 256 = 1.0
      timestampblocksize = 4*ktrg;

    ktrg = ( vdata.size() - timestampblocksize ) / IMG_WIDTH / IMG_HEIGHT; // fine

    cout << "  unpack " << ktrg << " event blocks";

    if( ktrg < ntrg ) {
      //outfile << endl << iev << " E incomplete"; // no endl here!
      cout << "  missing " << ntrg-ktrg << endl;
      Log.printf( "  trigger %i data size %i incomplete \n", iev, vdata.size() );
    }

    //cout << endl << "  timestamps:";
    cout << endl << "  hits:";

    unsigned pos = 0;

    for( int itrg = 0; itrg < ktrg; ++itrg ) {

      gSystem->ProcessEvents(); // ROOT

      unsigned long timestamp = 0;

      if( roc.ext && iFWVersion > 256 ) { // 256 = 1.0

	unsigned long ts1  = 0;
	unsigned long ts2  = 0;
	int trgid[4] = {0};
	for( size_t i = 0; i < 4; ++i ) {
	  trgid[i] = vdata.at(pos);
	  //std::cout << std::hex << trgid[i] << std::endl;
	  trgid[i] = trgid[i] & ~0xf000;
	  //timestamp = (timestamp & ~0xf000) | ((a & 0x300) >> 6);
	  ++pos;
	} // for i
	//std::cout << "End times tamps>>>>>>>>>>>>>>" << std::endl;

	ts1 =  ( trgid[1] << 12 ) + (trgid[0]);
	ts2 =  ( trgid[3] << 12 ) + (trgid[2]);
	timestamp =  ( ts2 << 24 ) + ts1; // 48 bits
	//cout << " " << timestamp;

      } // FW 1.1

      ++wev;

      double prevPH = 0;
      double pedAvg = 0;

      int nHit = 0; // number of hit pixel per event
      double hit[155][160] = {0};
      double stored[155][160] = {0};

      for( int col = 0; col < IMG_WIDTH; ++col ) { // column-wise

	for( unsigned row = 0; row < IMG_HEIGHT; ++row ) { // row-wise

	  int adc = vdata.at(pos);
	  ++pos;

	  if( adc & 0x1000 ) // ADC overrange
	    adc = R4S_VALUE_OR;
	  else if( adc & 0x0800 ) // negative
	    adc -= 0x1000;

	  double adcD = adc; // avoid type problems

	  // suppress fake pixels:

	  //if( col < 155 && row < 160 ) { // take all
	  if( col < 155 && row < 159 && row > 0 ) { // against noisy top and bot row

	    // Collect pedestal from the first nPedAvg events. CAUTION ASSUMES NO HITS IN THESE

	    if( wev <= nPedAvg )
	      runPed[col][row] += adcD / nPedAvg;

	    else { // hit finding

	      // calculate pulse height difference:

	      double thisPH = adcD - runPed[col][row]; // hits are negative

	      double diffPH = thisPH - prevPH; // first px in clus is negative

	      prevPH = thisPH;

	      savePH[col][row] = thisPH; // store

	      hadc.Fill( adc );
	      hadcSubPed.Fill( thisPH );
	      hadcNoCoMo.Fill( diffPH );

	      double signi = diffPH / noise[col][row];
	      hsigni.Fill( signi );

	      // seed finding using forward difference:

	      bool overThr = 0;

	      //if( diffPH <  thr ) { // thr is negative, leading edge of clus
	      if( signi <  thr ) { // thr is negative, leading edge of clus

		hit[col][row] = 1;
		++nHit;
		overThr = 1;
		hitAmplAvg.Fill( col+0.5, row+0.5, thisPH );

	      }

	      //else if( diffPH > -thr  && row > 0 ) { // thr is negative, trailing edge of clus
	      else if( signi > -thr  && row > 0 ) { // thr is negative, trailing edge of clus

		hit[col][row-1] = 1; // assign to previous
		++nHit;
		overThr = 1;

	      }

	      if( ! overThr ) { // update pedestal and RMS

		//if( !( thisPH < thr ) && !( thisPH > -thr ) ) { // excludes hits inside clusters as well
		if( !( thisPH < thr*noise[col][row] ) &&
		    !( thisPH > -thr*noise[col][row] ) ) { // excludes hits inside clusters as well

		  // Running Pedestal update:

		  runPed[col][row] = adcD / nPedAvg +
		    runPed[col][row] * (nPedAvg - 1.0) / nPedAvg;

		  // for RMS:
		  sumPHsq[col][row] += thisPH * thisPH;

		  notHitPH[col][row] += 1;

		}

		sumdPHsq[col][row] += diffPH  *diffPH;
		notHitPHC[col][row] += 1;

	      } // !overThr

	      if( wev%fupd == 0 ) { // update plots

		if( col==0 && row == 0 ) { //reset ONCE before fill
		  pedxy.Reset();
		  phRMS.Reset();
		  phCRMS.Reset();
		}

		pedxy.Fill( col+0.5, row+0.5, runPed[col][row] );
		phRMS.Fill( col+0.5, row+0.5, sqrt( sumPHsq[col][row] / notHitPH[col][row] ) );

		double dphrms = sqrt( sumdPHsq[col][row] / notHitPHC[col][row] );
		phCRMS.Fill( col+0.5, row+0.5, dphrms );
		hdphrms.Fill( dphrms );
		noise[col][row] = dphrms;

		// so they can be refilled

		sumPHsq[col][row] = 0;
		sumdPHsq[col][row] = 0;
		notHitPH[col][row] = 0;
		notHitPHC[col][row] = 0;

		pedAvg += runPed[col][row];
		dphrmsAvg += dphrms;

	      } // update plots

	    } // hit finding

	  } // valid px

	} // px rows

      } // col

      cout << " " << nHit;

      hitAvg += nHit;

      pixelOverThr.Fill(nHit);

      // task: make list of seeds

      bool hitFlag = 0;

      for( int col = 0; col< 155; ++col ) { // seed

	for( int row = 0; row< 160; ++row ) { // seed

	  if( hit[col][row] == 1 ) { // seed

	    if( hitFlag == 0 ) { // only once per event
	      outfile << endl // end previous event line
		      << wev << " F "
		      << iev
		      << " " << timestamp
		      << endl;
	      hitFlag = 1;
	    }

	    hitmap.Fill( col+0.5, row+0.5 );

	    // region of interest:

	    for( int i = col - roiCol; i <= col + roiCol; ++i ) {

	      for( int j = row - roiRow; j <= row + roiRow; ++j ) {

		// fiducial region:

		if( i >= 0 && i < 155 && j >= 0 && j < 160 ) {

		  // write pixel to file (once):

		  if( ! stored[i][j] ) {
		    outfile << i << " " << j << " " << -savePH[i][j] << " "; // invert PH
		    stored[i][j] = 1; // avoid double writing
		  }

		} // valid px

	      } // roi_rows j

	    } // roi_cols i

	  } // seed hit

	} // seed_row

      } // seed_col

      if( ! hitFlag ) // write empty line
	outfile << endl
		<< wev << " E "
		<< iev
		<< " " << timestamp; // no endl here !

      if( wev%fupd == 0 ) {

	int pntnr = wev/fupd;
	gPedAvg->SetPoint( pntnr, wev, pedAvg/160/155 );
	gHitAvg->SetPoint( pntnr, wev, hitAvg/fupd );
	grmsdph->SetPoint( pntnr, wev, dphrmsAvg/160/155 );
	hitAvg = 0;
	dphrmsAvg = 0;

	//h1->Modified();
	//h1->Update(); // hadc
	//h2->Modified();
	//h2->Update(); // hadcSubPed
	h3->Modified();
	h3->Update(); // hadcNoCoMo
	h4->Modified();
	h4->Update(); // hsigni
	//h5->Modified();
	//h5->Update();
	h6->Modified();
	h6->Update();
	/* slow 2D
	p1->Modified();
	p1->Update();
	p2->Modified();
	p2->Update();
	p3->Modified();
	p3->Update();
	p4->Modified();
	p4->Update();
	*/
	g1->Modified(); // ped vs time
	g1->Update();
	g2->Modified(); // hits vs time
	g2->Update();
	g3->Modified(); // noise vs time
	g3->Update();

	d1->Modified(); // hitmap
	d1->Update();

	if( wev > 1999 )
	  fupd = 500;
	if( wev > 9999 )
	  fupd = 2000;

      } // update

    } // ktrg

    // patch added event A if ktrig < ntrg

    for( int itrg = ktrg; itrg < ntrg; ++itrg )
      outfile << endl
	      << ++wev << " A "
	      << iev; // no endl here!

    cout << endl; // timestamps or hits

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s3 = ts.tv_sec; // seconds since 1.1.1970
    long f3 = ts.tv_nsec; // nanoseconds
    zeit2 += s3 - s2 + ( f3 - f2 ) * 1e-9; // read and cluster
    cout << " (in " << s3 - s2 + ( f3 - f2 ) * 1e-9 << " s)";
    cout << endl;

  } // run

  tb.Daq_Close();

  outfile.close();

  cout << hitmap.GetEntries() << " hits"
       << " => yield " << 1E2*hitmap.GetEntries()/(wev-100) << "%"
       << endl;

  cout << wev << " events written to " << fileName << endl;

  Log.printf( "  triggers %i, written %i\n", iev, wev );
  Log.flush();

  // Final update of some live plots:

  //h1->Modified();
  //h1->Update();
  //h2->Modified();
  //h2->Update();
  h3->Modified();
  h3->Update();
  h4->Modified();
  h4->Update();
  //h5->Modified();
  //h5->Update();
  h6->Modified();
  h6->Update();
  /*
  p1->Modified();
  p1->Update();
  p2->Modified();
  p2->Update();
  p3->Modified();
  p3->Update();
  p4->Modified();
  p4->Update();
  */
  g1->Modified();
  g1->Update();
  g2->Modified();
  g2->Update();
  g3->Modified();
  g3->Update();

  d1->Modified();
  d1->Update();

  // explicitely write TGraph:

  gPedAvg->Write();
  gHitAvg->Write();
  grmsdph->Write();

  // store all histograms to root file:

  histoFile->Write();
  histoFile->Close();

  delete histoFile;
  delete gPedAvg;
  delete gHitAvg;
  //delete h1;
  //delete h2;
  delete h3;
  delete h4;
  //delete h5;
  delete h6;
  /*
  delete p1;
  delete p2;
  delete p3;
  delete p4;
  */
  delete g1;
  delete g2;
  delete g3;
  delete d1;

  //app->Terminate(); // ends R4Stest

  //delete app; // don't delete singletons

  //app.Run(true); // want to display plots after running program?

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s9 = ts.tv_sec; // seconds since 1.1.1970
  long f9 = ts.tv_nsec; // nanoseconds
  double zeit = s9 - s0 + ( f9 - f0 ) * 1e-9; // total
  cout << "  test duration " << zeit << " s"
       << " (read " << zeit1 << ", process " << zeit2 << " s)"
       << " = " << wev / zeit << " Hz"
       << endl;

} // takedata

//------------------------------------------------------------------------------

CMD_PROC(tdp) // roi datataking - tryouts for parallelisation -- finn
{
  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

  omp_set_num_threads(2); // want 2 prongs

  timeval tv; // initialize container for time measurement

  gettimeofday( &tv, NULL ); // measure time
  long s0 = tv.tv_sec; // seconds since 1.1.1970
  long u0 = tv.tv_usec; // microseconds

  cout << s0 << " " << u0 << endl;

  #pragma omp parallel sections // fork
  {
    #pragma omp section // 1st prong
    {
      printf ("id = %d, \n", omp_get_thread_num());

      for( int i = 0; i < 1000; ++i){

        //printf("thread 1 %i \n", i); // print statements appear sequentially

        this_thread::sleep_for(chrono::microseconds(1000));

      }
    }
    #pragma omp section // 2nd prong
    {
      printf ("id = %d, \n", omp_get_thread_num());

      for( int i = 1000; i < 2000; ++i){

        //printf("thread 2 %i \n", i);

        this_thread::sleep_for(chrono::microseconds(1000));

      }
    }
  }

  gettimeofday( &tv, NULL );
  long s1 = tv.tv_sec; // seconds since 1.1.1970
  long u1 = tv.tv_usec; // microseconds

  cout << s1 << " " << u1 << endl;
  cout << s1-s0 << " " << u1-u0 << endl;

} // takedata - tryouts

//------------------------------------------------------------------------------
CMD_PROC(takeraw)
{
  int ntrg;
  if( ! PAR_IS_INT( ntrg, 1, 999 ) )
    ntrg = 1;

  int col_to_print;
  if( ! PAR_IS_INT( col_to_print, 0, 156 ) )
    col_to_print = 0;

  // run number from file:

  unsigned int run = 1;

  ifstream irunFile( "runNumber.dat" );
  irunFile >> run;
  irunFile.close();

  ++run;
  cout << "Run " << run << endl;

  ofstream orunFile( "runNumber.dat" );
  orunFile << run;
  orunFile.close();

  string fileName = Form( "raw%06i.txt", run );
  FILE * f = fopen( fileName.c_str(), "at" ); // append text

  if( !f ) return;

  line_wise = 0;
  IMG_HEIGHT = 162; // FW 0.7
  cout << "line_wise " << line_wise << endl;

  tb.Daq_Open( 48*1024*1024 ); // [words] max 64*1024*1024

  const uint32_t Blocksize = 4 * 1024 * 1024;

  // prepare ADC:

  tb.SignalProbeADC( PROBEA_SDATA1, GAIN_1 );

  tb.uDelay(500); // [us] Beat 23/11/2017: must be larger 400

  int iev = 0;
  int wev = 0; // write event

  timeval tv;
  gettimeofday( &tv, NULL );
  long s0 = tv.tv_sec; // seconds since 1.1.1970
  long u0 = tv.tv_usec; // microseconds

  while( !keypressed() ) {

    gettimeofday( &tv, NULL );
    long s1 = tv.tv_sec; // seconds since 1.1.1970
    long u1 = tv.tv_usec; // microseconds

    tb.Daq_Start(); // ADC -> RAM

    for( int itrg = 0; itrg < ntrg; ++itrg ) {

      ++iev;
      cout << " ev " << iev;

      tb.r4s_Start(); // R4S sequence

      //tb.uDelay(800); // 157*163*25 ns = 640 us R4S readout

    }

    if( roc.ext ) tb.uDelay(64*1000); // uint16 goes to 64k

    tb.Daq_Stop();

    // read buffer:

    vector<uint16_t> vdata;

    unsigned int ret = tb.Daq_Read( vdata, Blocksize );

    cout << "  status " << ret << " read " << vdata.size() << " words";

    gettimeofday( &tv, NULL );
    long s2 = tv.tv_sec; // seconds since 1.1.1970
    long u2 = tv.tv_usec; // microseconds
    cout << " in " << s2 - s1 + ( u2 - u1 ) * 1e-6 << " s";
    cout << endl;

    // unpack data:

    unsigned pos = 0;

    int ktrg  = vdata.size() / IMG_WIDTH /  IMG_HEIGHT;
    cout << "  unpack " << ktrg << " event blocks" << endl;

    for( int itrg = 0; itrg < ktrg; ++itrg ) {

      ++wev;
      fprintf( f, " %i event\n", wev );

      for( int x = 0; x < IMG_WIDTH; ++x ) { // column-wise

	if( x == col_to_print )
	  cout << "ev " << wev << " column " << x << " raw data:" << endl;

	for( unsigned y = 0; y < IMG_HEIGHT; ++y ) { // rows

	  int val = vdata.at(pos);
	  ++pos;

	  if( val & 0x1000 ) // ADC overrange
	    val = R4S_VALUE_OR;
	  else if( val & 0x0800 ) // negative
	    val -= 0x1000;

	  fprintf( f, " %4i", val );
	  if( x == col_to_print )
	    printf( " %4i", val );

	} // px rows

	fprintf( f, " col %i\n", x );
	//fputs( "\n", f ); // end of one col
	if( x == col_to_print )
	  printf( "\n" );

      } // col

    } // trg

  } // active

  tb.Daq_Close();

  fclose(f);
  cout << wev << " events written to " << fileName << endl;

  gettimeofday( &tv, NULL );
  long s9 = tv.tv_sec; // seconds since 1.1.1970
  long u9 = tv.tv_usec; // microseconds
  cout << "  test duration " << s9 - s0 + ( u9 - u0 ) * 1e-6 << " s" << endl;

} // takeraw

//------------------------------------------------------------------------------
CMD_PROC(getped)
{
  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

  TFile * histoFile = new TFile( "ped.root", "RECREATE" );

  TH1I hped( "ped",
	     "ped;ped;pixels",
	     100, 0, 1000 );

  TProfile2D pedxy( "pedxy",
		    "pedestal map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  TH2I hsuma( "suma",
	      "suma map;col;row;sum ADC",
	      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  TH2I hsumaa( "sumaa",
	       "sumaa map;col;row;sum ADC^2",
	       IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  TH2I hab( "ab",
	    "correlation 76 77 in row 88;pix 88.77 [ADC];pix 88.76 [ADC];events",
	    600, 100, 700, 600, 100, 700 );

  TH2I h8971( "c8971",
	      "correlation 89.71 to 89.70;pix 89.70 [ADC];pix 89.71 [ADC];events",
	      600, 100, 700, 600, 100, 700 );

  TH2I h88120( "c88120",
	       "correlation 88.120 89.70;pix 88.120 [ADC];pix 89.70 [ADC];events",
	       600, 100, 700, 600, 100, 700 );

  TH1I hamb( "amb",
	     "a-b;neighbour difference [ADC];pixel pairs",
	     101, -50.5, 50.5 );
  TH1I hapb( "apb",
	     "a+b;neighbour sum [ADC];pixel pairs",
	     101, 449.5, 550.5 );
  TH1I ha( "a",
	   "a;amplitude [ADC];pixels",
	   101, 199.5, 300.5 );

  R4sImg map;

  tb.r4s_SetPixCal( 159, 159 );

  if( line_wise )
    tb.r4s_SetSeqReadout( roc.ext ); // pedestal
  else
    tb.r4s_SetSeqReadCol( roc.ext );

  unsigned N = 99;

  for( unsigned i = 0; i < N; ++i ) { // events

    cout << "event " << i;

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) { // start top row

      for( unsigned x = 0; x < IMG_WIDTH; ++x ) { // horizontal = readout sequence = readout time

	int a = map.Get( x, y ); // [ADC]
	if( y < 160 && x < 155 )
	  hped.Fill( a );
	pedxy.Fill( x+0.5, y+0.5, a );
	hsuma.Fill( x+0.5, y+0.5, a );

	hsumaa.Fill( x+0.5, y+0.5, a*a );

	int b = a;
	if( x < 154 )
	  b = map.Get( x+1, y ); // [ADC]

	if( x == 77 && y == 88 ) {
	  hab.Fill( a+0.5, b+0,5 );
	  hamb.Fill( a-b );
	  hapb.Fill( a+b );
	  ha.Fill( a );
	}
	if( x == 71 && y == 89 ) {
	  h8971.Fill( a+0.5, b+0,5 );
	}
	if( x == 70 && y == 89 ) {
	  b = map.Get( 120, 88 ); // [ADC]
	  h88120.Fill( a+0.5, b+0,5 );
	}
      }
    }
  } // events

  TProfile2D rmsxy( "rmsxy",
		    "rms map;col;row;RMS [ADC]",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );

  TProfile rms88( "rms88",
		  "rms row 88;pixel;RMS [ADC]",
		  2*155, 0, 2*155, 0, 2222 );

  TH1I hrms( "rms",
	     "rms;rms [ADC];pixels",
	     100, 0, 20 );

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) { // start top row

    for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

      double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
      double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
      double rms = sqrt( sumaa/N - suma/N*suma/N );
      //cout << " " << rms;
      rmsxy.Fill( x+0.5, y+0.5, rms );
      if( y < 160 && x < 155 )
	hrms.Fill( rms );
      if( y == 88 && x < 155 )
	rms88.Fill( x+0.5, rms );
      if( y == 89 && x < 155 )
	rms88.Fill( 155+x+0.5, rms );
    }
    //cout << endl;
  }
  cout << "mean rms " << hrms.GetMean() << endl;
  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

}

//------------------------------------------------------------------------------
CMD_PROC(getcal) // ph-ped map
{
  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

  TFile * histoFile = new TFile( "cal.root", "RECREATE" );

  TProfile2D pedxy( "pedxy",
		    "pedestal map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  TH2I hsuma( "suma",
	      "suma map;col;row;sum ADC",
	      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  TH2I hsumaa( "sumaa",
	       "sumaa map;col;row;sum ADC^2",
	       IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  R4sImg map;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // take pedestal map:

  //tb.r4s_SetPixCal( 159, 159 );
  //tb.r4s_SetSeqReadout(roc.ext); // pedestal

  tb.r4s_SetVcal(0); // baseline?
  tb.r4s_SetSeqCalScan(); // Cal

  line_wise = 1;
  IMG_HEIGHT = 163; // Beat FW 0.6

  unsigned Np = 99;

  for( unsigned i = 0; i < Np; ++i ) { // events

    cout << "event " << i;

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

	int a = map.Get( x, y );
	pedxy.Fill( x+0.5, y+0.5, a );
	hsuma.Fill( x+0.5, y+0.5, a );
	hsumaa.Fill( x+0.5, y+0.5, a*a );

      }

  } // events

  TProfile2D rmsxy( "rmsxy",
		    "rms map;col;row;RMS [ADC]",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );

  TH1I hrms( "rms",
	     "rms;rms [ADC];pixels",
	     100, 0, 100 );

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

    for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

      double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
      double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
      double rms = sqrt( sumaa/Np - suma/Np*suma/Np );
      rmsxy.Fill( x+0.5, y+0.5, rms );
      if( y < 160 && x < 155 )
	hrms.Fill( rms );
    }

  cout << "mean rms " << hrms.GetMean() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tb.r4s_SetVcal(roc.Vcal);
  tb.r4s_SetSeqCalScan(); // Cal

  TProfile2D adcxy( "adcxy",
		    "pulse height;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  // taka data:

  unsigned Nc = 16;

  for( unsigned i = 0; i < Nc; ++i ) { // events

    cout << "event " << i;

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( unsigned x = 0; x < IMG_WIDTH; ++x )

	adcxy.Fill( x+0.5, y+0.5, map.Get( x, y ) );

  } // events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // subtract pedestal:

  TProfile2D phxy( "phxy",
		   "ph;col;row;<ADC> - ped",
		   IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  TH1I hph( "ph",
	     "ph;ph [ADC-ped];pixels",
	     100, -100, 900 );

  TH1I hcal( "cal",
	     "Vcal;Vcal [mV];pixels",
	     200, -100, 1900 );

  TProfile2D s2nxy( "s2nxy",
		    "pulse height / RMS;col;row;<PH>/RMS",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );
  TH1I hs2n( "s2n",
	     "signal.noise;PH/RMS;pixels",
	     100, 0, 200 );

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

    for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

      double a = adcxy.GetBinContent( x+1, y+1 ); // bins start at 1
      double p = pedxy.GetBinContent( x+1, y+1 ); // bins start at 1
      double s = rmsxy.GetBinContent( x+1, y+1 ); // bins start at 1
      phxy.Fill( x+0.5, y+0.5, p-a );
      s2nxy.Fill( x+0.5, y+0.5, (p-a)/s );
      if( y < 160 && x < 155 ) {
	hph.Fill( p-a ); // positive
	hs2n.Fill( (p-a)/s ); // positive
	hcal.Fill( PHtoVcal( p-a, x, y ) );
      }
    }

  cout << hph.GetName() << " mean " << hph.GetMean() << endl;
  cout << hcal.GetName() << " mean " << hcal.GetMean() << endl;

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

  // restore:

  if( line_wise )
    tb.r4s_SetSeqReadout( roc.ext ); // pedestal
  else
    tb.r4s_SetSeqReadCol( roc.ext );

}

//------------------------------------------------------------------------------
CMD_PROC(getcalcol) // ph-ped map column-wise
{
  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

  TFile * histoFile = new TFile( "cal.root", "RECREATE" );

  TProfile2D pedxy( "pedxy",
		    "pedestal map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  TH2I hsuma( "suma",
	      "suma map;col;row;sum ADC",
	      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  TH2I hsumaa( "sumaa",
	       "sumaa map;col;row;sum ADC^2",
	       IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  R4sImg map;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // take pedestal map:

  //tb.r4s_SetPixCal( 159, 159 );
  //tb.r4s_SetSeqReadout(roc.ext); // pedestal

  tb.r4s_SetVcal(0); // baseline?
  tb.r4s_SetSeqCalScanCol(); // column-wise cal sequence

  line_wise = 0;
  IMG_HEIGHT = 162; // FW >= 0.7

  unsigned Np = 99;

  for( unsigned i = 0; i < Np; ++i ) { // events

    cout << "event " << i;

    ReadImage( map );

    for( int x = IMG_WIDTH-1; x >= 0; --x ) // start last col

      for( unsigned y = 0; y < IMG_HEIGHT; ++y ) {

	int a = map.Get( x, y );
	pedxy.Fill( x+0.5, y+0.5, a );
	hsuma.Fill( x+0.5, y+0.5, a );
	hsumaa.Fill( x+0.5, y+0.5, a*a );

      }

  } // events

  TProfile2D rmsxy( "rmsxy",
		    "rms map;col;row;RMS [ADC]",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );

  TH1I hrms( "rms",
	     "rms;rms [ADC];pixels",
	     100, 0, 100 );

  for( int x = IMG_WIDTH-1; x >= 0; --x ) // start top row

    for( unsigned y = 0; y < IMG_HEIGHT; ++y ) {

      double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
      double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
      double rms = sqrt( sumaa/Np - suma/Np*suma/Np );
      rmsxy.Fill( x+0.5, y+0.5, rms );
      if( y < 160 && x < 155 )
	hrms.Fill( rms );
    }

  cout << "mean rms " << hrms.GetMean() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tb.r4s_SetVcal(roc.Vcal);
  tb.r4s_SetSeqCalScanCol(); // Cal

  TProfile2D adcxy( "adcxy",
		    "pulse height;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  // taka data:

  unsigned Nc = 16;

  for( unsigned i = 0; i < Nc; ++i ) { // events

    cout << "event " << i;

    ReadImage( map );

    for( int x = IMG_WIDTH-1; x >= 0; --x ) // start top row

      for( unsigned y = 0; y < IMG_HEIGHT; ++y )

	adcxy.Fill( x+0.5, y+0.5, map.Get( x, y ) );

  } // events

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // subtract pedestal:

  TProfile2D phxy( "phxy",
		   "ph;col;row;<ADC> - ped",
		   IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  TH1I hph( "ph",
	     "ph;ph [ADC-ped];pixels",
	     100, -100, 900 );

  TH1I hcal( "cal",
	     "Vcal;Vcal [mV];pixels",
	     200, -100, 1900 );

  TProfile2D s2nxy( "s2nxy",
		    "pulse height / RMS;col;row;<PH>/RMS",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );
  TH1I hs2n( "s2n",
	     "signal.noise;PH/RMS;pixels",
	     100, 0, 200 );

  for( int x = IMG_WIDTH-1; x >= 0; --x ) // start top row

    for( unsigned y = 0; y < IMG_HEIGHT; ++y ) {

      double a = adcxy.GetBinContent( x+1, y+1 ); // bins start at 1
      double p = pedxy.GetBinContent( x+1, y+1 ); // bins start at 1
      double s = rmsxy.GetBinContent( x+1, y+1 ); // bins start at 1
      phxy.Fill( x+0.5, y+0.5, p-a );
      s2nxy.Fill( x+0.5, y+0.5, (p-a)/s );
      if( y < 160 && x < 155 ) {
	hph.Fill( p-a ); // positive
	hs2n.Fill( (p-a)/s ); // positive
	hcal.Fill( PHtoVcal( p-a, x, y ) );
      }
    }

  cout << hph.GetName() << " mean " << hph.GetMean() << endl;
  cout << hcal.GetName() << " mean " << hcal.GetMean() << endl;

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

  // restore:

  if( line_wise )
    tb.r4s_SetSeqReadout( roc.ext ); // pedestal
  else
    tb.r4s_SetSeqReadCol( roc.ext );

}

//------------------------------------------------------------------------------
CMD_PROC(scancal) // scan Vcal
{
  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

  TFile * histoFile = new TFile( "scancal.root", "RECREATE" );

  TProfile2D pedxy( "pedxy",
		    "pedestal map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  TH2I hsuma( "suma",
	      "suma map;col;row;sum ADC",
	      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  TH2I hsumaa( "sumaa",
	       "sumaa map;col;row;sum ADC^2",
	       IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  TProfile2D rmsxy( "rmsxy",
		    "rms map;col;row;RMS [ADC]",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );

  TH1I hrms( "rms",
	     "rms;rms [ADC];pixels",
	     100, 0, 20 );

  // for pedestal:

  R4sImg map;

  //tb.r4s_SetPixCal( 159, 159 );
  //tb.r4s_SetSeqReadout(roc.ext); // pedestal

  tb.r4s_SetVcal(0); // baseline?
  tb.r4s_SetSeqCalScan(); // Cal

  line_wise = 1;
  IMG_HEIGHT = 163; // FW 0.6

  unsigned iev = 0;
  vector <double> via;
  via.reserve(999);

  unsigned Np = 99; // pedestals

  for( unsigned i = 0; i < Np; ++i ) { // events

    double ia = tb.GetIA();
    via.push_back(ia);
    ++iev;
    cout << "event " << iev << ", ia " << ia*1E3 << " mA";

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

	int a = map.Get( x, y );
	pedxy.Fill( x+0.5, y+0.5, a );
	hsuma.Fill( x+0.5, y+0.5, a );
	hsumaa.Fill( x+0.5, y+0.5, a*a );

      }

  } // events

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

    for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

      double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
      double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
      double rms = sqrt( sumaa/Np - suma/Np*suma/Np );
      rmsxy.Fill( x+0.5, y+0.5, rms );
      if( y < 160 && x < 155 )
	hrms.Fill( rms );
    }

  cout << "mean rms " << hrms.GetMean() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tb.r4s_SetSeqCalScan(); // Cal

  TProfile phvscal( "phvscal", "PH vs Vcal pix 77 88;Vcal [mV];PH-ped [ADC]",
		    240, 0, 2400 );

  int vstp = 10; // [mV]

  for( unsigned vcal = 0; vcal < 2402; vcal += vstp ) { // [mV]

    cout << "Vcal " << vcal << endl;

    tb.r4s_SetVcal(vcal);

    TProfile2D adcxy( Form( "adcxy_v%i", vcal ),
		      Form( "pulse height Vcal %i;col;row;<ADC>", vcal ),
		      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

    unsigned Nc = 16;
    for( unsigned i = 0; i < Nc; ++i ) { // events

      double ia = tb.GetIA();
      via.push_back(ia);
      ++iev;
      cout << "event " << iev << ", ia " << ia*1E3 << " mA";

      ReadImage( map );

      for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row
	for( unsigned x = 0; x < IMG_WIDTH; ++x )
	  adcxy.Fill( x+0.5, y+0.5, map.Get( x, y ) );

    } // events

    TProfile2D calxy( Form( "calxy_v%i", vcal ),
		      Form( "cal Vcal %i;col;row;<ADC> - ped", vcal ),
		      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );
    TH1I hcal( Form( "cal_v%i", vcal ),
	       Form( "cal Vcal %i;cal [ADC-ped];pixels", vcal ),
	       100, -100, 900 );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

	double a = adcxy.GetBinContent( x+1, y+1 ); // bins start at 1
	double p = pedxy.GetBinContent( x+1, y+1 ); // bins start at 1
	calxy.Fill( x+0.5, y+0.5, p-a );
	if( y < 160 && x < 155 )
	  hcal.Fill( p-a ); // positive
	if( x == 77 && y == 88 )
	  phvscal.Fill( vcal, p-a );
      }

    adcxy.Write();
    calxy.Write();
    hcal.Write();

    if( vcal == 200 ) vstp = 20;
    if( vcal == 800 ) vstp = 50;

  } // Vcal

  cout << via.size() << " events" << endl;

  TProfile iavsev( "iavsev",
		   "analog current;events;<IA> [mA]",
		   via.size(), 0, via.size(), 0, 500 );

  for( unsigned i = 0; i < via.size(); ++i )
    iavsev.Fill( i+0.5, via[i]*1E3 );

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

  // restore:

  cout << "Vcal back to " << roc.Vcal << endl;
  tb.r4s_SetVcal(roc.Vcal);

  if( line_wise )
    tb.r4s_SetSeqReadout( roc.ext ); // pedestal
  else
    tb.r4s_SetSeqReadCol( roc.ext );

} // scancal

//------------------------------------------------------------------------------
CMD_PROC(scanhold) // scan Vcal
{
  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

  TFile * histoFile = new TFile( "scanhold.root", "RECREATE" );

  TProfile2D pedxy( "pedxy",
		    "pedestal map;col;row;<ADC>",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

  TH2I hsuma( "suma",
	      "suma map;col;row;sum ADC",
	      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  TH2I hsumaa( "sumaa",
	       "sumaa map;col;row;sum ADC^2",
	       IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

  TProfile2D rmsxy( "rmsxy",
		    "rms map;col;row;RMS [ADC]",
		    IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );

  TH1I hrms( "rms",
	     "rms;rms [ADC];pixels",
	     100, 0, 20 );

  // for pedestal:

  R4sImg map;

  //tb.r4s_SetPixCal( 159, 159 );
  //tb.r4s_SetSeqReadout(roc.ext); // pedestal

  int Vcal = roc.Vcal;
  tb.r4s_SetVcal(0); // baseline?
  tb.r4s_SetSeqCalScan(); // Cal

  line_wise = 1;
  IMG_HEIGHT = 163; // FW 0.6

  unsigned iev = 0;
  vector <double> via;
  via.reserve(999);

  unsigned Np = 99; // pedestal

  for( unsigned i = 0; i < Np; ++i ) { // events

    double ia = tb.GetIA();
    via.push_back(ia);
    ++iev;
    cout << "event " << iev << ", ia " << ia*1E3 << " mA";

    ReadImage( map );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

	int a = map.Get( x, y );
	pedxy.Fill( x+0.5, y+0.5, a );
	hsuma.Fill( x+0.5, y+0.5, a );
	hsumaa.Fill( x+0.5, y+0.5, a*a );

      }

  } // events

  for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

    for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

      double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
      double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
      double rms = sqrt( sumaa/Np - suma/Np*suma/Np );
      rmsxy.Fill( x+0.5, y+0.5, rms );
      if( y < 160 && x < 155 )
	hrms.Fill( rms );
    }

  cout << "mean rms " << hrms.GetMean() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  tb.r4s_SetVcal(Vcal);
  tb.r4s_SetSeqCalScan(); // Cal

  TProfile phvshld( "phvshld",
		    "PH vs hold all pix;hold [6.25 ns];<PH-ped> [ADC]",
		    256, -0.5, 255.5, -2222, 2222 );

  int stp = 1; // [6.25 ns]

  for( unsigned hld = 0; hld < 256; hld += stp ) {

    cout << "hold " << hld << endl;

    tb.r4s_SetHoldPos(hld);

    TProfile2D adcxy( Form( "adcxy_hld%i", hld ),
		      Form( "pulse height hold %i;col;row;<ADC>", hld ),
		      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

    unsigned Nc = 16;
    for( unsigned i = 0; i < Nc; ++i ) { // events

      double ia = tb.GetIA();
      via.push_back(ia);
      ++iev;
      cout << "event " << iev << ", ia " << ia*1E3 << " mA";

      ReadImage( map );

      for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row
	for( unsigned x = 0; x < IMG_WIDTH; ++x )
	  adcxy.Fill( x+0.5, y+0.5, map.Get( x, y ) );

    } // events

    TProfile2D calxy( Form( "calxy_hld%i", hld ),
		      Form( "cal hold %i;col;row;<ADC> - ped", hld ),
		      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );
    TH1I hcal( Form( "cal_hld%i", hld ),
	       Form( "cal hold %i;cal [ADC-ped];pixels", hld ),
	       100, -100, 900 );

    for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

      for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

	double a = adcxy.GetBinContent( x+1, y+1 ); // bins start at 1
	double p = pedxy.GetBinContent( x+1, y+1 ); // bins start at 1
	calxy.Fill( x+0.5, y+0.5, p-a );
	if( y < 160 && x < 155 ) {
	  hcal.Fill( p-a ); // positive
	  phvshld.Fill( hld, p-a );
	}
      }

    //adcxy.Write();
    calxy.Write();
    hcal.Write();

    if( hld ==  50 ) stp = 2;
    if( hld == 100 ) stp = 4;

  } // hold

  cout << via.size() << " events" << endl;

  TProfile iavsev( "iavsev",
		   "analog current;events;<IA> [mA]",
		   via.size(), 0, via.size(), 0, 500 );

  for( unsigned i = 0; i < via.size(); ++i )
    iavsev.Fill( i+0.5, via[i]*1E3 );

  cout << "hold back to " << roc.Hold << endl;
  tb.r4s_SetHoldPos(roc.Hold);

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

  // restore:

  if( line_wise )
    tb.r4s_SetSeqReadout( roc.ext ); // pedestal
  else
    tb.r4s_SetSeqReadCol( roc.ext );

} // scanhold

//------------------------------------------------------------------------------
CMD_PROC(scanhold2d) // scan hold for various VgPr and VgSh -- finn
{
  // scanhold, with additional loops over vgpr and vgsh
  // and without the map and histogram at every delay to save some disk space

  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

    //get current settings to restore after scan
  int svVgPr = roc.VgPr;
  int svVgSh = roc.VgSh;

  for( int vpr = 400; vpr < 801; vpr += 100 ){   // for real
    for( int vsh = 400; vsh < 901; vsh += 100 ){
  //for( int vpr = 700; vpr < 801; vpr += 100 ){  // for testing
  //  for( int vsh = 700; vsh < 801; vsh += 100 ){

      tb.r4s_SetRgpr(vpr);
      roc.VgPr = vpr;
      tb.r4s_SetRgsh(vsh);
      roc.VgSh = vsh;

      TFile * histoFile = new TFile( Form( "schld2d/c%i_scanhold_pr%i_sh%i.root" , roc.chip, vpr, vsh ), "RECREATE" );

      TProfile2D pedxy( "pedxy",
            "pedestal map;col;row;<ADC>",
            IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

      TH2I hsuma( "suma",
            "suma map;col;row;sum ADC",
            IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

      TH2I hsumaa( "sumaa",
            "sumaa map;col;row;sum ADC^2",
            IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT );

      TProfile2D rmsxy( "rmsxy",
            "rms map;col;row;RMS [ADC]",
            IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, 0, 2222 );

      TH1I hrms( "rms",
          "rms;rms [ADC];pixels",
          100, 0, 20 );

      // for pedestal:

      R4sImg map;

      tb.r4s_SetVcal(0); // baseline?
      tb.r4s_SetSeqCalScan(); // Cal

      line_wise = 1;
      IMG_HEIGHT = 163; // FW 0.6

      unsigned iev = 0;
      vector <double> via;
      via.reserve(999);

      unsigned Np = 99;

      for( unsigned i = 0; i < Np; ++i ) { // events

        double ia = tb.GetIA();
        via.push_back(ia);
        ++iev;
        cout << "event " << iev << ", ia " << ia*1E3 << " mA";

        ReadImage( map );

        for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

          for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

            int a = map.Get( x, y );
            pedxy.Fill( x+0.5, y+0.5, a );
            hsuma.Fill( x+0.5, y+0.5, a );
            hsumaa.Fill( x+0.5, y+0.5, a*a );

          }

      } // events

      for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

        for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

          double suma = hsuma.GetBinContent( x+1, y+1 ); // bins start at 1
          double sumaa = hsumaa.GetBinContent( x+1, y+1 ); // bins start at 1
          double rms = sqrt( sumaa/Np - suma/Np*suma/Np );
          rmsxy.Fill( x+0.5, y+0.5, rms );
          if( y < 160 && x < 155 )
            hrms.Fill( rms );

        }

      cout << "mean rms " << hrms.GetMean() << endl;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      tb.r4s_SetVcal(400);
      tb.r4s_SetSeqCalScan(); // Cal

      TProfile phvshld( "phvshld",
            "PH vs hold all pix;hold [6.25 ns];<PH-ped> [ADC]",
            256, -0.5, 255.5, -2222, 2222 );

      int stp = 2; // [12.5 ns]

      for( unsigned hld = 0; hld < 256; hld += stp ) {

        cout << "hold " << hld << endl;

        tb.r4s_SetHoldPos(hld);

        TProfile2D adcxy( Form( "adcxy_hld%i", hld ),
              Form( "pulse height hold %i;col;row;<ADC>", hld ),
              IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );

        unsigned Nc = 16;
        for( unsigned i = 0; i < Nc; ++i ) { // events

          double ia = tb.GetIA();
          via.push_back(ia);
          ++iev;
          cout << "event " << iev << ", ia " << ia*1E3 << " mA";

          ReadImage( map );

          for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row
            for( unsigned x = 0; x < IMG_WIDTH; ++x )
              adcxy.Fill( x+0.5, y+0.5, map.Get( x, y ) );

        } // events

          // save some disc space
        //TProfile2D calxy( Form( "calxy_hld%i", hld ),
        //      Form( "cal hold %i;col;row;<ADC> - ped", hld ),
        //      IMG_WIDTH, 0, IMG_WIDTH, IMG_HEIGHT, 0, IMG_HEIGHT, -2222, 2222 );
        //TH1I hcal( Form( "cal_hld%i", hld ),
        //    Form( "cal hold %i;cal [ADC-ped];pixels", hld ),
        //    100, -100, 900 );

        for( int y = IMG_HEIGHT-1; y >= 0; --y ) // start top row

          for( unsigned x = 0; x < IMG_WIDTH; ++x ) {

            double a = adcxy.GetBinContent( x+1, y+1 ); // bins start at 1
            double p = pedxy.GetBinContent( x+1, y+1 ); // bins start at 1
            //calxy.Fill( x+0.5, y+0.5, p-a );
            if( y < 160 && x < 155 ) {
              //hcal.Fill( p-a ); // positive
              phvshld.Fill( hld, p-a );
            }
          }

        //adcxy.Write();
        //calxy.Write();
        //hcal.Write();

        if( hld ==  80 ) stp =  4;
        if( hld == 120 ) stp = 10;

      } // hold

      cout << via.size() << " events" << endl;

      TProfile iavsev( "iavsev",
          "analog current;events;<IA> [mA]",
          via.size(), 0, via.size(), 0, 500 );

      for( unsigned i = 0; i < via.size(); ++i )
        iavsev.Fill( i+0.5, via[i]*1E3 );

      cout << "hold back to " << roc.Hold << endl;
      tb.r4s_SetHoldPos(roc.Hold);

      histoFile->Write();
      histoFile->Close();
      cout << histoFile->GetName() << endl;

      // restore:

      if( line_wise )
        tb.r4s_SetSeqReadout( roc.ext ); // pedestal
      else
        tb.r4s_SetSeqReadCol( roc.ext );

    }
  } // VgPr and VgSh

    // restore VgPr and VgSh settings
  tb.r4s_SetRgpr( svVgPr );
  roc.VgPr = svVgPr;
  tb.r4s_SetRgsh( svVgSh );
  roc.VgSh = svVgSh;

} // scanhold2d

//------------------------------------------------------------------------------
CMD_PROC(seqreadout)
{
  int ext;
  if( !PAR_IS_INT( ext, 0, 1 ) ) // external trigger flag
    ext = 0; // default = internal

  tb.r4s_SetSeqReadout(ext);
  roc.ext = ext;
  line_wise = 1;
  IMG_HEIGHT = 163; // FW 0.6
  cout << "line_wise " << line_wise << endl;
  DO_FLUSH;
}

//------------------------------------------------------------------------------
CMD_PROC(seqreadcol) // FW 0.7
{
  int ext;
  if( !PAR_IS_INT( ext, 0, 1 ) ) // external trigger flag
    ext = 0; // default = internal

  tb.r4s_SetSeqReadCol(ext);
  roc.ext = ext;
  line_wise = 0;
  IMG_HEIGHT = 162; // FW 0.7
  cout << "line_wise " << line_wise << endl;
  DO_FLUSH;
}

//------------------------------------------------------------------------------
CMD_PROC(seqcalscan)
{
  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

  tb.r4s_SetSeqCalScan();
  line_wise = 1;
  IMG_HEIGHT = 163; // FW 0.6
  cout << "line_wise " << line_wise << endl;
  DO_FLUSH;
}

//------------------------------------------------------------------------------
CMD_PROC(scanva) // IA vs VA
{
  int dummy; if( !PAR_IS_INT( dummy, 0, 1 ) ) dummy = 80;

  TFile * histoFile = new TFile( "scanva.root", "RECREATE" );

  TProfile iavsva( "iavsva",
		   "analog current;VA [mV];<IA> [mA]",
		   150, 1000, 2500, 0, 500 );

  for( int va = 1005; va < 2502; va += 10 ) { // [mV]

    tb._SetVA(va);

    cout << "VA " << va;

    for( int i = 0; i < 9; ++i ) {

      tb.uDelay(10000);
      double ia = tb.GetIA();
      iavsva.Fill( va, ia*1E3 );
      cout << " " << ia*1E3;
    }

    cout << ", VA " << tb.GetVA() << endl;

  } // va

  histoFile->Write();
  histoFile->Close();
  cout << histoFile->GetName() << endl;

} // scanva
//------------------------------------------------------------------------------
CMD_PROC(asetvana) // auto set vana -- finn
{
  int iath; if( !PAR_IS_INT( iath, 0, 200 ) ) iath = 124;
  double iath_d = (double)iath;

  cout << "target Iana: " << iath_d << endl;

  int good = 0;
  for( int va = 1900; va < 2300; va += 5 ) { // [mV]

    tb.r4s_SetVana(va);
    roc.Vana = va;

    cout << "VA " << va << " - ";

    double iavg = 0;

    for( int i = 0; i < 5; ++i ) {

      tb.uDelay(10000);
      double ia = tb.GetIA();
      iavg += ia*1E3/5;
      cout << " " << ia*1E3;

    }

    cout << " - avg. Iana " << iavg << " mA" << endl;

    if( iavg > iath_d ){
      good = 1;
    }
    if( iavg > iath_d + 10){
      good = 2;
    }
    if( good )
      break;

  } // va

  if( good == 0 )
    cout << "FAIL: Target Iana not reached! Set vana manually!" << endl;
  if( good == 1 )
    cout << "Reached target Iana!" << endl;
  if( good == 2 )
    cout << "FAIL: Iana above target! Set vana manually!" << endl;

} // asetvana

/*
//------------------------------------------------------------------------------
CMD_PROC(gui)
{
printf("Connect to GUI ...");
CDataLink link(1024, 65536);
//	ShellExecute(NULL, "open", "file.exe", NULL, NULL, SW_SHOWDEFAULT);

system ("start roc4sens_view\\bin\\Release\\roc4sens_view.exe");
//	system ("start roc4sens_view.exe");

link.Run();
printf("Connection to GUI closed.\n");
}
*/
