// This code loops over all of the events in a given root file and makes plots
// It is run as a function analysisSolar.cc+(path/inputFileName.root, path/outputFileName.root)
#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DS/UniversalTime.hh>
#include <RAT/SunUtil.hh>
#include <RAT/DataCleaningUtility.hh>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h> 
#include <string>

void AnalysisN16mpw()
{
  const double ITRval = 0.55;const double ITR_PMT_Hi = 5.0; const double ITR_PMT_Lo = -2.5;
  double driveP0 = 1, driveP1 = 0;//driveP0 = 0.997143, driveP1 = -64.2801;
  const double c_light = 299.792458;
  double grVelocity = c_light/1.38486;//2.17554021555098529e+02 ;
  double nhitCut = 25;
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
  ifstream in0, in1;
  Double_t xtrue,ytrue,rtrue,ztrue,rt;
  ofstream out,outputbadfiles;
  ostringstream oss;
  in0.open("fileList.dat");
  
  // Define all of the histograms and ntuples that are to be made  
  //MPW no drive correction, with drive correction and ITR+beta14 cuts
  TH1D *hITR = new TH1D("hITR","trig+FECD, ITR",100,0,1.0);
  TH1D *hBeta14 = new TH1D("hBeta14","trig+FECD, beta14",200,-1.0,4.0);
  TH2D *hITR_vs_Beta14 = new TH2D("hITR_vs_Beta14","trig+FECD, itr vs beta14",100,0,1.0,200,-1.0,4.0);
  
  TH1D* radius = new TH1D( "radius", "Radius of all triggered events", 100, 0., 15000.); 
  TH1D* nHit = new TH1D( "nHit", "nHit of all triggered events", 1000, 0.5, 1000.5);   
  TH1D* nHitAV = new TH1D( "nHitAV", "nHit of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5);  
  TH2D* nHit_toSun = new TH2D("nHit_toSun", "nHit vs direction to the Sun of all triggered events", 1000, 0.5, 1000.5, 100, -1., 1.);
  TH2D* nHit_toSunAV = new TH2D("nHit_toSunAV", "nHit vs direction to the Sun of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5, 100, -1., 1.); 
  //TH3D* nHit_toSun_radius = new TH3D("nHit_toSun_radius", "nHit vs direction to the Sun vs radius of all triggered events", 1000, 0.5, 1000.5, 100, -1., 1., 100., 0., 15000.);

  TNtuple* candidateEventInfo = new TNtuple("candidateEventInfo", "Information for interesting events, without drive correction", "runNumber:eventGTID:nHit:radius_mm:cosThetaToSun:x_pos:y_pos:z_pos:x_dir:y_dir:z_dir");
  
  TTree *tree= new TTree("T","Simulation");
  
  double posX, posY, posZ, posRad, time, dirX, dirY, dirZ, nhits, sunDirX, sunDirY, sunDirZ, cosThetaToSun, ITR, beta14, itr, iso;
  double posTheta, posPhi;//for berkeley blob  
  double prompt_cut; 
  UInt_t day, sec, runNumber, eventGTID, fecd;
  double nsecs;
  std::vector<double>* vtRes;
  
  tree->Branch("posX",&posX,"posX/D");
  tree->Branch("posY",&posY,"posY/D");
  tree->Branch("posZ",&posZ,"posZ/D");
  tree->Branch("posTheta",&posTheta,"posTheta/D");
  tree->Branch("posPhi",&posPhi,"posPhi/D");
  tree->Branch("time",&time,"time/D");
  tree->Branch("timeRes",&vtRes);
  tree->Branch("posRad",&posRad,"posRad/D");
  tree->Branch("dirX",&dirX,"dirX/D");
  tree->Branch("dirY",&dirY,"dirY/D");
  tree->Branch("dirZ",&dirZ,"dirZ/D");
  tree->Branch("nhits",&nhits,"nhits/D");
//  tree->Branch("sunDirX",&sunDirX,"sunDirX/D");
//  tree->Branch("sunDirY",&sunDirY,"sunDirY/D");
//  tree->Branch("sunDirZ",&sunDirZ,"sunDirZ/D");
//  tree->Branch("cosThetaToSun",&cosThetaToSun,"cosThetaToSun/D");
  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("eventGTID",&eventGTID,"eventGTID/i");
  tree->Branch("day",&day,"day/i");
  tree->Branch("sec",&sec,"sec/i");
  tree->Branch("fecd",&fecd,"fecd/i");
  tree->Branch("nsecs",&nsecs,"nsecs/D");
  tree->Branch("itr",&itr,"itr/D");
  tree->Branch("beta14",&iso,"iso/D");
  tree->Branch("prompt_cut",&prompt_cut,"prompt_cut/D");

  char filenames[1500];
  Int_t countValid = 0; 
  Int_t i=0;
  Int_t TriggerType ;
  Int_t nentries=0;
  ULong64_t bitword ;
  while(in0>>filenames)
  {
   cout<<" filenames "<<filenames<<endl ;
    // Load the RAT file
   //TFile *ftemp = new TFile(filenames,"read");
   //check whether file is Zombie or not contains Tree
   //bool checkCorrupt = ftemp->IsZombie();
   //bool checkContain = ftemp->GetListOfKeys()->Contains("T");
   //if(checkCorrupt || !checkContain){cout<<"File "<<filenames<<" is corrupted, passing it!"<<endl;}
   //else
   {
    RAT::DU::DSReader dsReader(filenames);
    
    // To later get the run info
    const RAT::DS::Run& run = dsReader.GetRun();
    
    size_t numEvents = dsReader.GetEntryCount();
   
    cout<<numEvents<<endl; 
    // Loop over all of the events
    for(size_t iEntry = 0; iEntry < numEvents; iEntry++ )
    {
       //posX=-10000, posY=-10000, posZ=-10000, posRad=-10000, time=99999, ITR=-2, beta14=-10, dirX=-2.0, dirY=-2.0, dirZ=-2.0, nhits=0, sunDirX=0.0,sunDirY=0.0, sunDirZ=0.0, cosThetaToSun=-2.0, itr=-2, iso=-10, prompt_cut = -1;
       //posTheta = -10;  posPhi = -10;
       //for drive correction, only for position !!!
       //posXcor=-10000, posYcor=-10000, posZcor=-10000, posRadCor=-10000;
       //day=0, sec=0, runNumber=0, eventGTID=0;
       //nsecs=0 ;	
       const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
       // const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event 
       
       size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
       
       ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

       if(iEntry%5000 == 0) {cout<<iEntry<<endl;}
	
    	// Looping over triggered events in each ds
	if (iEVCount) {
    	 for(size_t iEV=0; iEV<iEVCount; iEV++) //iEVCount
         {
	     // Get the event
	     const RAT::DS::EV& ev = rDS.GetEV(iEV);
	     const RAT::DS::CalPMTs& calpmts = ev.GetCalPMTs();
	     runNumber = run.GetRunID();
	     eventGTID = ev.GetGTID();
	     RAT::DS::DataQCFlags dcflags = ev.GetDataCleaningFlags();
             if( RAT::EventIsClean( ev, dcAnalysisWord ) )
	     {	  
 		  // Fill the nHit histogram, all events should have nHit
		  Float_t nHitOfEvent = ev.GetNhits();
		  nHit->Fill(nHitOfEvent);
		  nhits= nHitOfEvent ;
		  TriggerType = ev.GetTrigType() ;
		  //cout<<" TriggerType "<<TriggerType<<endl ;
		  vtRes->clear(); 
		  // Now get the radius, if applicable
		  Double_t radiusOfEvent = -1., radiusOfEventCor = -1.;
		  TVector3 pos_fit, pos_cor;
		  TVector3 u_fit;
		  const string fitName = "MultiPathProcessor";

                  for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
                  {
                    fecd = calpmts.GetFECDPMT(ipmt).GetID();
                  }

		  if(((TriggerType&2)==2))// || ((TriggerType&18)==18))
		  {
		    if( ev.FitResultExists(fitName) ) 
		    {  // It needs to exist 
                     RAT::DS::FitVertex fVertex =ev.GetFitResult("MultiPathProcessor").GetVertex(0);
                     pos_fit = fVertex.GetPosition(); 
		     u_fit = fVertex.GetDirection();		
		     //NOTE: drive correction here
		     pos_fit = driveP0*pos_fit + driveP1*u_fit;
		     if(pos_fit.Mag()<9000) 
		     {  //fit valid: inside cavity, specified by MPW fitter		    
                        countValid++;
                        RAT::DS::FitResult fResult = ev.GetFitResult(fitName);     // Get fit result
                        time = fVertex.GetTime();
                        radiusOfEvent = pos_fit.Mag();
                        posX = pos_fit.X();posY = pos_fit.Y();posZ = pos_fit.Z();
                        posTheta = pos_fit.CosTheta(); posPhi = pos_fit.Phi();
                        posRad = radiusOfEvent;
                        radius->Fill(radiusOfEvent);
                        dirX = u_fit.X();dirY = u_fit.Y();dirZ = u_fit.Z();		

                        //!!! apply ITR and iso cuts here, use pos_fit or pos_cor 
                        double countPMT_ITR = 0 ;double p1 = 0.0; double p4 = 0.0;
                        for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                        {
                          TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                          double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
                          double tRes = (calpmts.GetPMT(ipmt)).GetTime()-fVertex.GetTime()-(pmtpos-pos_fit).Mag()/grVelocity;
                          vtRes->push_back(tRes);
                          if(tRes<ITR_PMT_Hi && tRes>ITR_PMT_Lo) countPMT_ITR++;
                          //for beta14 cuts
                          //const TVector3 pmtDir1 = (pmtpos - pos_cor).Unit();
                            
		          //do promt time cut
		          prompt_cut = abs((pos_fit-pmtpos).Mag() - c_light/1.40*(hitTime-fVertex.GetTime())); 
			  //cout<<prompt_cut<<" "<<(pos_fit-pmtpos).Mag()<<" "<<hitTime-fVertex.GetTime()<<" "<<c_light/1.40*(hitTime-fVertex.GetTime())<<endl;
		          if(prompt_cut<=2000)
 		          {
  		            const TVector3 pmtDir1 = (pmtpos - pos_fit).Unit();
                            for(unsigned int ipmt2 = ipmt+1; ipmt2 <calpmts.GetCount(); ipmt2++ )
                            {
                               //const TVector3 pmtDir2 = (pmtInfo.GetPosition(calpmts.GetPMT(ipmt2).GetID()) - pos_cor).Unit();
                               const TVector3 pmtDir2 = (pmtInfo.GetPosition(calpmts.GetPMT(ipmt2).GetID()) - pos_fit).Unit();
                               const double cosTheta12 = pmtDir1.Dot( pmtDir2 );
                               p1 += cosTheta12; // Legendre Polynomial order 1
                               const double cosThetaSquared = cosTheta12 * cosTheta12;
                               p4 += ( 35.0 * cosThetaSquared * cosThetaSquared - 30.0 * cosThetaSquared + 3.0 ) / 8.0; // Legendre Polynomial order 4
                            }
                          }
			 }
                        p1 = 2.0 * p1 / static_cast<double>(calpmts.GetCount() * (calpmts.GetCount() - 1 ) );
                        p4 = 2.0 * p4 / static_cast<double>(calpmts.GetCount() * (calpmts.GetCount() - 1 ) );
                        const double beta14 = p1 + 4 * p4;
                        iso = beta14;
                        double ITR = countPMT_ITR/double(calpmts.GetCount()); 
                        itr = ITR;	
                        hITR->Fill(ITR);hBeta14->Fill(beta14);hITR_vs_Beta14->Fill(ITR,beta14);
 
                        tree->Fill();
                       } // fit Region Valid 
       
                       // Fill the nHit histogram if the event was inside the AV
                       if (radiusOfEvent >= 0. && radiusOfEvent <= 6000.) {
                          nHitAV->Fill(nHitOfEvent);
                        }
                         
                     // **** Here you can put code to retrieve the ITR and Beta14 to select on those along with a the statement (radiusOfEvent <= 5500.)   
                         
                     // **** You can change the below ntuple to suit, adding in the universal time, it has everything else already
                  candidateEventInfo->Fill(run.GetRunID(), ev.GetGTID(), nHitOfEvent, radiusOfEvent, cosThetaToSun, pos_fit.X(),  pos_fit.Y(), pos_fit.Z(), u_fit.X(),  u_fit.Y(),  u_fit.Z());  
		 }//fit Name valid 		 
	    } // trig cuts
	   } //data cleaning
	 }//loop triggered events
    }//triggered events
   }//loop events 
   }//check file is zombie or not contains Tree
  } // end of the while loop for reading files
  
  
  // Give everything an axis label, as should always be
  radius->SetXTitle( "Radius [mm]" );radius->SetYTitle( "Counts" );
  nHit->SetXTitle( "nHit" );nHit->SetYTitle( "Counts" );  
  nHitAV->SetXTitle( "nHit" );nHitAV->SetYTitle( "Counts" );  
  nHit_toSun->SetXTitle( "nHit" );nHit_toSun->SetYTitle( "cos(#theta) to Sun" );
  nHit_toSun->SetZTitle( "Counts" );
  nHit_toSunAV->SetXTitle( "nHit" );nHit_toSunAV->SetYTitle( "cos(#theta) to Sun" );nHit_toSunAV->SetZTitle( "Counts" );
  //nHit_toSun_radius->SetXTitle( "nHit" );nHit_toSun_radius->SetYTitle( "cos(#theta) to Sun" );nHit_toSun_radius->SetZTitle( "Radius [mm]" );
  
  cout<<countValid<<endl;
  // Write the histograms to fileName
  TString newname = "MPWProcess_";
  TString processname = newname+filenames;
  TFile *file=new TFile(processname,"RECREATE");
  file->cd();
  tree->Write();
  hITR->Write();hBeta14->Write();hITR_vs_Beta14->Write();
  radius->Write();nHit->Write();nHitAV->Write();nHit_toSun->Write();nHit_toSunAV->Write();//nHit_toSun_radius->Write();
  candidateEventInfo->Write();
  
}
