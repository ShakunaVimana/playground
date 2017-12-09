//Fit peak with Gaussian and look left and right by 3 sigma; 
//throw away tails until fit well with mean+-3*sigma

#include <vector>
#include "TH2.h"
#include "TH1.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TTree.h"
#include "iostream"
#include "TString.h"

using namespace std;

double calXRMS(TH1F* hTempG);
double calYRMS(TH1F* hTempG);
double calZRMS(TH1F* hTempG);
const double SIGMA = 3;//what sigma you want to cut? 1 sigma or 3 sigma?
const double STEP = 50;//new rms to old rms step [mm]

void cut_tail()
{
  TString infile_name = "ResolPmtSelect_FitJPMPWpdf1_mode-10to100_SNOP_100934_000.root";
  TString s1 = infile_name(0,5);
  TString s2 = infile_name(33,6);
  if(s1!="Resol") {std::cout<<"Error: wrong file. Using Resol* "<<std::endl; exit(0);}
  TFile *fp = new TFile(infile_name,"READ");
  
  TH1F *hfitX=(TH1F*)fp->Get("hfitX");
  TH1F *hfitY=(TH1F*)fp->Get("hfitY");
  TH1F *hfitZ=(TH1F*)fp->Get("hfitZ");
  
  TString out_name = "CheckRMS_"+ s2 +".root";
  TFile *fRMS = new TFile(out_name,"recreate");
  
  const int bin_tot = 2000;//originally, 2000 bins, -10000 to 10000
  std::cout<<"File: "<<infile_name<<std::endl;
  std::cout<<"x   xRMS  y   yRMS	z	zRMS"<<std::endl;
  double xmean = hfitX->GetMean();double xrms= hfitX->GetRMS();
  double ymean = hfitY->GetMean();double yrms= hfitY->GetRMS();
  double zmean = hfitZ->GetMean();double zrms= hfitZ->GetRMS();
  //double nhitMean = hNhit->GetMean();
  std::cout<<xmean<<" "<<xrms<<" "<<ymean<<" "<<yrms<<" "<<zmean<<" "<<zrms<<std::endl;
  
  double oldRMS, newRMS;
  TAxis *xaxis = hfitX->GetXaxis();
  int bin_mean = xaxis->FindBin(xmean);
  int bin_lo = xaxis->FindBin(xmean-SIGMA*xrms);
  int bin_hi = xaxis->FindBin(xmean+SIGMA*xrms);
  double lo = xmean-SIGMA*xrms;
  double hi = xmean+SIGMA*xrms;
 // std::cout<<bin_mean<<"   "<<bin_lo<<"   "<<bin_hi<<"   "<<std::endl;
  
  int binNew_tot =bin_hi-bin_lo+1;
  TH1F *hXtemp =  new TH1F("htempX","temptX",binNew_tot,lo,hi) ;
  
  for(int m=0; m<binNew_tot; m++) {
    hXtemp->SetBinContent(m+1,hfitX->GetBinContent(m+bin_lo));
  }
  //check rms converge
  TH1F *hh = calXRMS(hXtemp)->Clone();
  double oldXrms = hXtemp->GetRMS();
  double newXrms = hh->GetRMS();
  while(1)
  {  
	  //cout<<"new X rms "<<newXrms<<endl; 
	  if(abs(newXrms-oldXrms)<STEP) 
	  {
		  //cout<<"rms converge: "<<newXrms<<endl;
		  TH1F *hXsmooth =(TH1F*) hh->Clone();
		  hXsmooth->SetName("hXsmooth");hXsmooth->SetTitle("fit X after smooth");
		  delete hh;
		  break;}
	  else{
		 oldXrms = newXrms; 
	     TH1F *hh = calXRMS(hh); 
	     newXrms = hh->GetRMS();	  
		 flag = 1 ;
		 continue;  
	  } 
  }
  
  double oldRMS1, newRMS1;
  TAxis *xaxis1 = hfitY->GetXaxis();
  int bin_meanY = xaxis->FindBin(ymean);
  int bin_loY = xaxis->FindBin(ymean-SIGMA*yrms);
  int bin_hiY = xaxis->FindBin(ymean+SIGMA*yrms);
  double loY = ymean-SIGMA*yrms;
  double hiY = ymean+SIGMA*yrms;
  int binNew_totY =bin_hiY-bin_loY+1;
  TH1F *hYtemp =  new TH1F("htempY","temptY",binNew_totY,loY,hiY) ;
  
  for(int m=0; m<binNew_totY; m++) {
    hYtemp->SetBinContent(m+1,hfitY->GetBinContent(m+bin_loY));
  }
  //check rms converge
  TH1F *hhY = calYRMS(hYtemp)->Clone();
  double oldYrms = hYtemp->GetRMS();
  double newYrms = hhY->GetRMS();
  while(1)
  {  
	  //cout<<"new Y rms "<<newYrms<<endl; 
	  if(abs(newYrms-oldYrms)<STEP) 
	  {
		  //cout<<"rms converge: "<<newXrms<<endl;
		  TH1F *hYsmooth =(TH1F*) hhY->Clone();
		  hYsmooth->SetName("hYsmooth");hYsmooth->SetTitle("fit Y after smooth");
		  delete hhY;
		  break;}
	  else{
		 oldYrms = newYrms; 
	     TH1F *hhY = calYRMS(hhY); 
	     newYrms = hhY->GetRMS();	  
		 flag = 1 ;
		 continue;  
	  } 
  }

  double oldRMS2, newRMS2;
  TAxis *xaxis2 = hfitZ->GetXaxis();
  int bin_meanZ = xaxis->FindBin(zmean);
  int bin_loZ = xaxis->FindBin(zmean-SIGMA*zrms);
  int bin_hiZ = xaxis->FindBin(zmean+SIGMA*zrms);
  double loZ = zmean-SIGMA*zrms;
  double hiZ = zmean+SIGMA*zrms;  
  int binNew_totZ =bin_hiZ-bin_loZ+1;
  TH1F *hZtemp =  new TH1F("htempZ","temptZ",binNew_totZ,loZ,hiZ) ;
  
  for(int m=0; m<binNew_tot; m++) {
    hZtemp->SetBinContent(m+1,hfitZ->GetBinContent(m+1+bin_loZ));
  }
  //check rms converge
  TH1F *hhZ = calZRMS(hZtemp)->Clone();
  double oldZrms = hZtemp->GetRMS();
  double newZrms = hhZ->GetRMS();
  while(1)
  {  
    if(abs(newZrms-oldZrms)<STEP) 
    {
      //cout<<"rms converge: "<<newZrms<<endl;
      TH1F *hZsmooth =(TH1F*) hhZ->Clone();
      hZsmooth->SetName("hZsmooth");  hZsmooth->SetTitle("fit Z after smooth");
      delete hhZ;
      break;
     }
    else{
      oldZrms = newZrms; 
      TH1F *hhZ = calZRMS(hhZ); 
      newZrms = hhZ->GetRMS();	  
      flag = 1 ;
      continue;  
    } 
  }
  
  delete hXtemp;
  delete hYtemp;
  delete hZtemp;
//Here just use histo RMS and Mean  
//  double xMu = hXsmooth->GetMean();double yMu = hYsmooth->GetMean();double zMu = hZsmooth->GetMean();
//  double xrms_new = hXsmooth->GetRMS();double yrms_new = hYsmooth->GetRMS();double zrms_new = hZsmooth->GetRMS();
//  cout<<xMu<<" "<<xrms_new<<" "<<yMu<<" "<<yrms_new<<" "<<" "<<zMu<<" "<<zrms_new<<endl;

//now use Gaussian Fit !!
 
   hXsmooth->Fit("gaus","q");
   hYsmooth->Fit("gaus","q");
   hZsmooth->Fit("gaus","q");
 
   TF1 *fitGausX = hXsmooth->GetFunction("gaus");   
   TF1 *fitGausY = hYsmooth->GetFunction("gaus");   
   TF1 *fitGausZ = hZsmooth->GetFunction("gaus");   
   
   double xMu = fitGausX->GetParameter(1);double xrms_new = fitGausX->GetParameter(2);
   double yMu = fitGausY->GetParameter(1);double yrms_new = fitGausY->GetParameter(2);
   double zMu = fitGausZ->GetParameter(1);double zrms_new = fitGausZ->GetParameter(2);

//for table fill-in
  std::cout<<xmean<<" "<<xrms<<" "<<ymean<<" "<<yrms<<" "<<zmean<<" "<<zrms<<" "<<xMu<<" "<<xrms_new<<" "<<yMu<<" "<<yrms_new<<" "<<zMu<<" "<<zrms_new<<std::endl;

  fRMS->cd();
  hfitX->Write();hXsmooth->Write(); 
  hfitY->Write();hYsmooth->Write(); 
  hfitZ->Write();hZsmooth->Write(); 
 
}

TH1F* calXRMS(TH1F* hXTempG)
{
    double Xrms;
    TAxis *xaxis1 = hXTempG->GetXaxis();
    double tempXmean = hXTempG->GetMean();
 //   cout<<" tempX mean "<<tempXmean<<endl;
    double oldXrms = hXTempG->GetRMS();
    double temp_lo = tempXmean-SIGMA*oldXrms;
    double temp_hi = tempXmean+SIGMA*oldXrms;
    int binTemp_lo = xaxis1->FindBin(tempXmean-SIGMA*oldXrms);
    int binTemp_hi = xaxis1->FindBin(tempXmean+SIGMA*oldXrms); 
    int binTemp_tot = binTemp_hi-binTemp_lo+1; 
    TH1F *hXtemp_new =  new TH1F("htemp_new","tempt new",binTemp_tot,temp_lo,temp_hi); 
    for(int m=0; m<binTemp_tot; m++) {
     hXtemp_new->SetBinContent(m+1,hXTempG->GetBinContent(m+binTemp_lo));
    } 	  
    Xrms = hXtemp_new->GetRMS(); 
    return hXtemp_new;
}

TH1F* calYRMS(TH1F* hYTempG)
{
    double Yrms;
    TAxis *xaxis1 = hYTempG->GetXaxis();
    double tempXmean = hYTempG->GetMean();
    //cout<<" tempX mean "<<tempXmean<<endl;
    double oldXrms = hYTempG->GetRMS();
    double temp_lo = tempXmean-SIGMA*oldXrms;
    double temp_hi = tempXmean+SIGMA*oldXrms;
    int binTemp_lo = xaxis1->FindBin(tempXmean-SIGMA*oldXrms);
    int binTemp_hi = xaxis1->FindBin(tempXmean+SIGMA*oldXrms); 
    int binTemp_tot = binTemp_hi-binTemp_lo+1; 
    TH1F *hYtemp_new =  new TH1F("htemp_new","tempt new",binTemp_tot,temp_lo,temp_hi); 
    //cout<<binTemp_tot<<"   "<<temp_lo<<"  "<<temp_hi<<endl;
    for(int m=0; m<binTemp_tot; m++) {
     hYtemp_new->SetBinContent(m+1,hYTempG->GetBinContent(m+binTemp_lo));
    } 	  
    Yrms = hYtemp_new->GetRMS(); 
    return hYtemp_new;
}

TH1F* calZRMS(TH1F* hZTempG)
{
    double Zrms;
    TAxis *xaxis1 = hZTempG->GetXaxis();
    double tempZmean = hZTempG->GetMean();
   // cout<<" tempX mean "<<tempXmean<<endl;
    double oldZrms = hZTempG->GetRMS();
    double temp_lo = tempZmean-SIGMA*oldZrms;
    double temp_hi = tempZmean+SIGMA*oldZrms;
    int binTemp_lo = xaxis1->FindBin(tempZmean-SIGMA*oldZrms);
    int binTemp_hi = xaxis1->FindBin(tempZmean+SIGMA*oldZrms); 
    int binTemp_tot = binTemp_hi-binTemp_lo+1; 
    TH1F *hZtemp_new =  new TH1F("htemp_new","tempt new",binTemp_tot,temp_lo,temp_hi); 
    for(int m=0; m<binTemp_tot; m++) {
     hZtemp_new->SetBinContent(m+1,hZTempG->GetBinContent(m+binTemp_lo));
    } 	  
    Zrms = hZtemp_new->GetRMS(); 
    return hZtemp_new;
}
