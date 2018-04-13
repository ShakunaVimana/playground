// This code loops over all of the events in a given root file and makes plots
// It is run as a function analysisSolar.cc+(path/inputFileName.root, path/outputFileName.root)
//#include <TH1D.h>
//#include <TH2D.h>
//#include <TH3D.h>
//#include <TFile.h>
//#include <TNtuple.h>
//#include <TTree.h>
//#include <fstream>
//#include <sstream>
//#include <iostream>
//#include <stdlib.h> 
//#include <string>
//#include <vector>
//void AnalysisRatN16()
{
  
  ifstream in0, in1;
  ofstream out,outputbadfiles;
  ostringstream oss;
  in0.open("dataFile.dat");
  int nlines = 0;
  double firstX = 368.443633;
  double lastX = -128.443668;
  int bin = floor(firstX-lastX);
  TFile *f = new TFile("NMRdata.root","recreate");
  double point,yintensity,xHz,xppm;
  TH1D *hspec = new TH1D("hspec","NMR",bin,firstX*-1,lastX*-1); 
  while(1)
  {
    in0>>point>>yintensity>>xHz>>xppm;
    xppm = xppm*-1;
    cout<<xppm<<" "<<yintensity<<endl;
    if(!in0.good())break;
    hspec->Fill(xppm,yintensity);
    nlines++;
  }
  in0.close();
  hspec->Draw();
  f->Write();
  f->Close();
}
