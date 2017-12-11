#include "TH1D.h"
#include <TVector3.h>
#include <TMath.h>
#include <TROOT.h>
#include "TFile.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TCanvas.h"
#include "TPad.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"

double posResol(double *x, double *par)
{   //par[0]=aE,par[1]=sigma_P, par[2]=mu_P, par[3]=tau_P, par[4]=scale
	double funcPosResol =(1.0-par[0])/(TMath::Sqrt(2*TMath::Pi())*par[1])*TMath::Exp(-0.5*pow((x[0]-par[2])/(par[1]),2))+par[0]/(2*par[3])*TMath::Exp(-TMath::Abs(x[0]-par[2])/(par[3]));
	
    return funcPosResol*par[4];
}

double angularResol(double *x, double *par)
{
    //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale
	double jillings = par[0]*par[1]*TMath::Exp(par[1]*(x[0]-1))/(1-TMath::Exp(-2*par[1]))+(1-par[0])*par[2]*TMath::Exp(par[2]*(x[0]-1))/(1-TMath::Exp(-2*par[2]));
    return jillings*par[3];	
}

/// for calculating cosTheta at xx%
double inte_angularResol(double *x, double *par)//an integral of the angularResol, to calculate cosTheta
{
    //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale
	double inte_jillings = par[0]*TMath::Exp(par[1]*(x[0]-1))/(1-TMath::Exp(-2*par[1]))+(1-par[0])*TMath::Exp(par[2]*(x[0]-1))/(1-TMath::Exp(-2*par[2]));
    return inte_jillings*par[3];	
}
double fSolve(double *x, double *par)//an integral of the angularResol, to calculate cosTheta at xx%
{
   //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale, par[4]=F1, par[5]=percent1, par[6] = integral
	double solve_jillings = par[3]*par[0]*TMath::Exp(par[1]*(x[0]-1))/(1-TMath::Exp(-2*par[1]))+(1-par[0])*TMath::Exp(par[2]*(x[0]-1))/(1-TMath::Exp(-2*par[2]));
    double result = solve_jillings-par[4]+par[5]*par[6];
    return result;	
}

///

//Energy fitting function
double fEresol(double *x, double *par)
{
	
	
}

//the exact function Boulay used in thesis
double boulayPosResol(double *x, double *par)
{   //
	double funcPosResol =(1.0-par[0])/(TMath::Sqrt(2*TMath::Pi())*par[1])*TMath::Exp(-0.5*pow((x[0]-par[2])/(par[1]),2))+par[0]/(2*par[3])*TMath::Exp(-TMath::Abs(x[0]-par[2])/(par[3]));
	
    return funcPosResol*par[4];
}

double boulayAngularResol(double *x, double *par)
{
    //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale

	double jillings = par[0]*par[1]*TMath::Exp(par[1]*(x[0]-1))/(1-TMath::Exp(-2*par[1]))+(1-par[0])*par[2]*TMath::Exp(par[2]*(x[0]-1))/(1-TMath::Exp(-2*par[2]));
    return jillings*par[3];	
}


void SaveAngular()
{
  double ae,sig_p,mu_p,tau_p;//par[0],par[1],par[2],par[3]
  TString infile_name = "ResolMCnoDrive_FitMPW1p40_N16_rat631_x0100_y0100_z0100.root";
  TFile *f1 = new TFile(infile_name);    
  //TFile *f6 = new TFile("convolved_N16.root");   
  //TH1F *hConvol = (TH1F*)f6->Get("hConvl"); hConvol->Scale(1e4);
  //hConvol->Sumw2();

  TF1 *funcPosResol = new TF1("fitPosResol",posResol,-200,200,5);
  TF1 *funcAngularResol = new TF1("fitAngularResol",angularResol,-0.5,1.0,4);
 
  //From SNO data//FTU
  TH1F *hpx = (TH1F*)f1->Get("hfitX");
  TH1F *hAngular = (TH1F*)f1->Get("hcosThetaSNOP_MultiWater");
  TH1F *hAngularCut1 = (TH1F*)f1->Get("hcosThetaCut1SNOP_MultiWater");
  TH1F *hAngularCut2 = (TH1F*)f1->Get("hcosThetaCut2SNOP_MultiWater");
  TH1F *hAngularCut3 = (TH1F*)f1->Get("hcosThetaCut3SNOP_MultiWater");
  TH1F *hAngularCut4 = (TH1F*)f1->Get("hcosThetaCut4SNOP_MultiWater");
  TH1F *hAngularCut5 = (TH1F*)f1->Get("hcosThetaCut5SNOP_MultiWater");
  
///fit Position Resolution Sets initial values and parameter names
   //hpx->SetStats(1);
   //cout<<"integrate "<<hpx->Integral()<<endl;
   //funcPosResol->SetParameter(0,0.55);//funcPosResol->FixParameter(0,0.55);
   //funcPosResol->SetParameter(1,17.0);
   //funcPosResol->SetParameter(2,-0.5);
   //funcPosResol->SetParameter(3,32.0);
   //funcPosResol->SetParameter(4,1e5);//1e5
   //funcPosResol->SetParNames("#alpha_{e}","#sigma_{P}","#mu_{P}","#tau_{P}","scale");
   
   //hpx->Sumw2();// hpx->Rebin(2.4); //hpx->Scale(1./5);
   //hpx->GetXaxis()->SetRangeUser(-300,300);
   
  // hpx->Fit(funcPosResol,"RP");
  //hConvol->Fit(funcPosResol,"RP");
   ///fit Angular Resolution 
   double trial_scale = hAngularCut5->Integral(-0.5,1.0);
   funcAngularResol->FixParameter(0,0.5);
   funcAngularResol->SetParameters(0.5,4.0,18.8,trial_scale);
   funcAngularResol->SetParNames("#alpha_{M}","#beta_{M}","#beta_{S}","scale");
   ///Rebin histograms 
   hAngular->Rebin(10.0);hAngularCut1->Rebin(10.0);hAngularCut2->Rebin(10.0);hAngularCut3->Rebin(10.0);hAngularCut4->Rebin(10.0);hAngularCut5->Rebin(10.0);    
   hAngular->Sumw2();hAngularCut1->Sumw2();hAngularCut2->Sumw2();hAngularCut3->Sumw2();hAngularCut4->Sumw2();hAngularCut5->Sumw2();
// Fit histogram in range defined by function
   //hAngular->Fit(funcAngularResol,"P");
   //hAngularCut1->Fit(funcAngularResol,"P");
   //hAngularCut2->Fit(funcAngularResol,"P");
   //hAngularCut3->Fit(funcAngularResol,"P");
   //hAngularCut4->Fit(funcAngularResol,"P");
   hAngularCut5->Fit(funcAngularResol,"PRq");//q for quiet
   
   double ndf = funcAngularResol->GetNDF();
   double chi2 = funcAngularResol->GetChisquare();    
  
   TF1 *inte_funcAngularResol = new TF1("inte_fitAngularResol",inte_angularResol,-1.0,1.0,4);
   //fitting parameters for Integral: par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale
   double aM1, bM1, bS1, scale1;
   aM1 = funcAngularResol->GetParameter(0);
   bM1 = funcAngularResol->GetParameter(1);
   bS1 = funcAngularResol->GetParameter(2);
   scale1 = funcAngularResol->GetParameter(3);
   inte_funcAngularResol->SetParameters(aM1,bM1,bS1,scale1);
  //check the integral function is right
  // cout<<inte_funcAngularResol->Eval(1.0)-inte_funcAngularResol->Eval(-1.0)<<" "<<funcAngularResol->Integral(-1.0,1.0)<<endl;
  
   double integral =  inte_funcAngularResol->Eval(1.0)-inte_funcAngularResol->Eval(-1.0); 
   double percent1 = 0.5;double percent2 = 0.8;double percent3 = 0.95;
   double F1 = inte_funcAngularResol->Eval(1.0);
   //Solve function: F(x)-F(1) + percent1*integral = 0
   //par[0]=aM, par[1]=bM, par[2]=bS, par[3]=scale, par[4]=F1, par[5]=percent1, par[6] = integral
   TF1 SolvFunc("solve_func",fSolve,-1.0,1.0,7);
   SolvFunc.SetParameters(aM1,bM1,bS1,scale1,F1,percent1,integral);
   ROOT::Math::WrappedTF1 wf1(SolvFunc); 
   // Create the Integrator
   ROOT::Math::BrentRootFinder brf;
   // Set parameters of the method
   brf.SetFunction( wf1, -1.0, 1.0);
   brf.Solve();
   double cosTheta50 = brf.Root();
  // cout << "cosTheta50% "<<brf.Root() << endl;
   
   SolvFunc.SetParameters(aM1,bM1,bS1,scale1,F1,percent2,integral);
   ROOT::Math::WrappedTF1 wf2(SolvFunc); 
   // Create the Integrator
   ROOT::Math::BrentRootFinder brf2;
   // Set parameters of the method
   brf2.SetFunction( wf2, -1.0, 1.0);
   brf2.Solve();
   double cosTheta80 = brf2.Root();
   //cout << "cosTheta80% "<<brf2.Root() << endl;

   SolvFunc.SetParameters(aM1,bM1,bS1,scale1,F1,percent3,integral);
   ROOT::Math::WrappedTF1 wf3(SolvFunc); 
   // Create the Integrator
   ROOT::Math::BrentRootFinder brf3;
   // Set parameters of the method
   brf3.SetFunction( wf3, -1.0, 1.0);
   brf3.Solve();
   double cosTheta95 = brf3.Root();
  // cout << "cosTheta95% "<<brf3.Root() << endl;
   TString savefilename = "SaveFit"+infile_name;
   TFile *savefile = new TFile(savefilename,"recreate");
   std::cout<<infile_name<<" "<<funcAngularResol->GetParameter(1)<<" "<<funcAngularResol->GetParError(1)<<" "<<funcAngularResol->GetParameter(2)<<" "<<funcAngularResol->GetParError(2)<<" "<<chi2<<"/"<<ndf<<" "<<cosTheta50<<" "<<cosTheta80<<" "<<cosTheta95<<std::endl;
   savefile->cd();
   hAngularCut5->Write();
   funcAngularResol->Write();
///The Exact Boulay Function    
   //TF1 *funcBoulayPosResol = new TF1("fBoulayPosResol",boulayPosResol,-200,200,5);
   //funcBoulayPosResol->SetParameters(0.55,17.23,-0.5414,32.25,1.0e5);
   //funcBoulayPosResol->SetParNames("#alpha_{e}","#sigma_{P}","#mu_{P}","#tau_{P}","scale");
   //cout<<funcBoulayPosResol->Integral(-200,200)<<endl;
   
   //TF1 *funcBoulayAngularResol = new TF1("fBoulayAngularResol",boulayAngularResol,-0.5,1.0,4);
   //funcBoulayAngularResol->SetParameters(0.5,4.0,18.8,200);
   //cout<<funcBoulayAngularResol->Integral(-1,1)<<endl;
   //funcBoulayAngularResol->SetParNames("#alpha_{M}","#beta_{M}","#beta_{S}","scale");
   //gPad->Update();
/// 
   //TPaveStats *st = (TPaveStats*)hpx->FindObject("stats");
   //st->SetOptStat(1111);
   //hpx->SetStats(1).
   
   //gStyle->SetOptStat(1111);
//   gStyle->SetOptFit(1);
   //hpx->SetStats(1);
   //hpx->Write();
   //hConvol->Write();
   //hAngular->SetStats(1);
 //  hAngular->Write();hAngularCut1->Write();hAngularCut2->Write();hAngularCut3->Write();hAngularCut4->Write();     
 //  funcBoulayPosResol->Write(); funcBoulayAngularResol->Write();

}
