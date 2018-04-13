/// \file
/// \ingroup tutorial_spectrum
/// \macro_image
/// \macro_code
/// \author Shakuna Vimana
/// \Modified from ROOT Spectrum/peaks.C for NMR analysis 
#include "TCanvas.h"
#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

Int_t npeaks = 30;
Double_t fpeaks(Double_t *x, Double_t *par) {
   Double_t result = par[0] + par[1]*x[0];
   for (Int_t p=0;p<npeaks;p++) {
      Double_t norm  = par[3*p+2];
      Double_t mean  = par[3*p+3];
      Double_t sigma = par[3*p+4];
      result += norm*TMath::Gaus(x[0],mean,sigma);
   }
   return result;
}
void peaks(Int_t np=10) {
   npeaks = TMath::Abs(np);
   double region1_lo = -50, region1_hi = 0;
   double region2_lo = -110, region2_hi = -50;
   double region3_lo = -160, region3_hi = -110;
   double region4_lo = -220, region4_hi = -160;

   ifstream in0, in1;
   ofstream out,outputbadfiles;
   ostringstream oss;
   TString filename = "dataFile8.txt";
   in0.open(filename);
   int nlines = 0;
   double firstX = 368.443633;
   double lastX = -128.443668;
   int bin = floor(firstX-lastX);
   //cout<<bin<<endl;
   TString numb = filename(8); 
   TFile *file = new TFile("NMRdata"+numb+".root","recreate");
   double point,yintensity,xHz,xppm;
   TH1D *h = new TH1D("hspec","NMR",500,firstX*-1,lastX*-1);
   TH1D *hsimple = new TH1D("hspec_1","NMR",500,firstX*-1,lastX*-1); 
   while(1)
   {
     in0>>point>>yintensity>>xHz>>xppm;
     xppm = xppm*-1;
     //cout<<xppm<<" "<<yintensity<<endl;
     if(!in0.good())break;
     h->Fill(xppm,yintensity);
     hsimple->Fill(xppm,yintensity);
     nlines++;
   }
   in0.close();

   TAxis *xaxis = h->GetXaxis();
   int binx1_lo = xaxis->FindBin(region1_lo);
   int binx1_hi = xaxis->FindBin(region1_hi);
   int binx2_lo = xaxis->FindBin(region2_lo);
   int binx2_hi = xaxis->FindBin(region2_hi);
   int binx3_lo = xaxis->FindBin(region3_lo);
   int binx3_hi = xaxis->FindBin(region3_hi);
   int binx4_lo = xaxis->FindBin(region4_lo);
   int binx4_hi = xaxis->FindBin(region4_hi);
   TH1F *hregion1 = new TH1F("hregion1","sub region1:alkylcarbon",int(binx1_hi-binx1_lo+1),region1_lo,region1_hi);
   TH1F *hregion2 = new TH1F("hregion2","sub region2:alkoxycarbonyl",int(binx2_hi-binx2_lo+1),region2_lo,region2_hi);
   TH1F *hregion3 = new TH1F("hregion3","sub region3:aromaticcarbon",int(binx3_hi-binx3_lo+1),region3_lo,region3_hi);
   TH1F *hregion4 = new TH1F("hregion4","sub region4:carbonylcarbon",int(binx4_hi-binx4_lo+1),region4_lo,region4_hi);

   for(int i =0;i<int(binx1_hi-binx1_lo)+1;i++) 
   {
    hregion1->SetBinContent(i+1,h->GetBinContent(binx1_lo+1 + i));
   }

   for(int i =0;i<int(binx2_hi-binx2_lo)+1;i++)
   {
    hregion2->SetBinContent(i+1,h->GetBinContent(binx2_lo+1 + i));
   }

   for(int i =0;i<int(binx3_hi-binx3_lo)+1;i++)
   {
    hregion3->SetBinContent(i+1,h->GetBinContent(binx3_lo+1 + i));
   } 
  
   for(int i =0;i<int(binx4_hi-binx4_lo)+1;i++)
   {
    hregion4->SetBinContent(i+1,h->GetBinContent(binx4_lo + i));
   }
     
   //generate n peaks at random
   Double_t par[3000];
   par[0] = 0.8;
   par[1] = -0.6/1000;
   Int_t p;
   for (p=0;p<npeaks;p++) {
      par[3*p+2] = 1;
      par[3*p+3] = 10+gRandom->Rndm()*980;
      par[3*p+4] = 3+2*gRandom->Rndm();
   }
   TF1 *f = new TF1("f",fpeaks,0,1000,2+3*npeaks);
   f->SetNpx(1000);
   f->SetParameters(par);
   TCanvas *c1 = new TCanvas("c1","c1",10,10,1000,900);
   c1->Divide(1,2);
   c1->cd(1);
   h->Draw();
   TH1F *h2 = (TH1F*)h->Clone("h2");
   //Use TSpectrum to find the peak candidates
   TSpectrum *s = new TSpectrum(2*npeaks);
   Int_t nfound = s->Search(h,6,"",0.1);
   printf("Found %d candidate peaks to fit\n",nfound);
   //Estimate background using TSpectrum::Background
   TH1 *hb = s->Background(h,20,"same");
   if (hb) c1->Update();
   if (np <0) return;

   //estimate linear background using a fitting method
   c1->cd(2);
   TF1 *fline = new TF1("fline","pol1",-1*firstX,-1*lastX);
   h->Fit("fline","qn");
   //Loop on all found peaks. Eliminate peaks at the background level
   par[0] = fline->GetParameter(0);
   par[1] = fline->GetParameter(1);
   npeaks = 0;
   Double_t *xpeaks = s->GetPositionX();
   for (p=0;p<nfound;p++) {
      Double_t xp = xpeaks[p];
      Int_t bin = h->GetXaxis()->FindBin(xp);
      Double_t yp = h->GetBinContent(bin);
      if (yp-TMath::Sqrt(yp) < fline->Eval(xp)) continue;
      par[3*npeaks+2] = yp;
      par[3*npeaks+3] = xp;
      par[3*npeaks+4] = 3;
      npeaks++;
   }

  //cut the background
   TH1D *h_cut = new TH1D("hspec_cut","NMR",500,firstX*-1,lastX*-1);
   for(int i=0;i<500;i++)
   {
    double xi = h->GetXaxis()->GetBinCenter(i+1);
    double yi = h->GetBinContent(i+1);
    if(yi>fline->Eval(xi))
    h_cut->SetBinContent(i+1,h->GetBinContent(i+1));
    else 
    h_cut->SetBinContent(i+1,0);
   }

   double x1 = h->Integral(binx1_lo,binx1_hi);
   double x2 = h->Integral(binx2_lo,binx2_hi);
   double x3 = h->Integral(binx3_lo,binx3_hi);
   double x4 = h->Integral(binx4_lo,binx4_hi);
   double all = h->Integral();

   double x1_cut = h_cut->Integral(binx1_lo,binx1_hi);
   double x2_cut = h_cut->Integral(binx2_lo,binx2_hi);
   double x3_cut = h_cut->Integral(binx3_lo,binx3_hi);
   double x4_cut = h_cut->Integral(binx4_lo,binx4_hi);
   double all_cut = h_cut->Integral();
   cout<<"样品"<<numb<<endl;
   cout<<"去除基线前: "<<endl;
   cout<<"烷基碳比例: "<<x1/all<<" 烷氧碳比例: "<<x2/all<<" 芳香碳比例: "<<x3/all<<" 羰基碳比例: "<<x4/all<<endl;
   cout<<"去除基线后: "<<endl;
   cout<<"烷基碳比例: "<<x1_cut/all_cut<<" 烷氧碳比例: "<<x2_cut/all_cut<<" 芳香碳比例: "<<x3_cut/all_cut<<" 羰基碳比例: "<<x4_cut/all_cut<<endl;
   file->cd();
   h->Write();hsimple->Write();hregion1->Write();hregion2->Write();hregion3->Write();hregion4->Write();
   h_cut->Write();
   printf("Found %d useful peaks to fit\n",npeaks);
   printf("Now fitting: Be patient\n");
   TF1 *fit = new TF1("fit",fpeaks,0,1000,2+3*npeaks);
   //we may have more than the default 25 parameters
   TVirtualFitter::Fitter(h2,10+3*npeaks);
   fit->SetParameters(par);
   fit->SetNpx(1000);
   h2->Fit("fit");
}
