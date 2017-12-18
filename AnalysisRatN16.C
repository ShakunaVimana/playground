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
#include <vector>
void AnalysisRatN16()
{
  
  ifstream in0, in1;
  Double_t xtrue,ytrue,rtrue,ztrue,rt;
  ofstream out,outputbadfiles;
  ostringstream oss;
  in0.open("ratfile.dat");
  double grVelocity = 2.17554021555098529e+02;
  const RAT::DU::PMTInfo& pmtInfo = RAT::DU::Utility::Get()->GetPMTInfo();
  
  // Define all of the histograms and ntuples that are to be made  
  TH1D* radius = new TH1D( "radius", "Radius of all triggered events", 100, 0., 15000.);
  
  TH1D* nHit = new TH1D( "nHit", "nHit of all triggered events", 1000, 0.5, 1000.5);   
  TH1D* nHitAV = new TH1D( "nHitAV", "nHit of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5); 
  
  TH2D* nHit_toSun = new TH2D("nHit_toSun", "nHit vs direction to the Sun of all triggered events", 1000, 0.5, 1000.5, 100, -1., 1.);
  TH2D* nHit_toSunAV = new TH2D("nHit_toSunAV", "nHit vs direction to the Sun of all triggered events reconstructed in the AV", 1000, 0.5, 1000.5, 100, -1., 1.);
  
  TH3D* nHit_toSun_radius = new TH3D("nHit_toSun_radius", "nHit vs direction to the Sun vs radius of all triggered events", 1000, 0.5, 1000.5, 100, -1., 1., 100., 0., 15000.);
  
  TNtuple* candidateEventInfo = new TNtuple("candidateEventInfo", "Information for interesting events", "runNumber:eventGTID:nHit:radius_mm:cosThetaToSun:x_pos:y_pos:z_pos:x_dir:y_dir:z_dir");
  
  TTree *tree= new TTree("T","SolarNu");
  
  double posX, posY, posZ, posRad, time, energy, ITR, beta14, dirX, dirY, dirZ, nhits, sunDirX, sunDirY, sunDirZ, cosThetaToSun, itr, iso, thij;
  double posTheta, posPhi;
  UInt_t day, sec, runNumber, eventGTID, fecd;
  double nsecs;
  std::vector<double>* vtRes;
  
  tree->Branch("posX",&posX,"posX/D");
  tree->Branch("posY",&posY,"posY/D");
  tree->Branch("posZ",&posZ,"posZ/D");
  tree->Branch("posTheta",&posTheta,"posTheta/D");
  tree->Branch("posPhi",&posPhi,"posPhi/D");
  tree->Branch("time",&time,"time/D");
//  tree->Branch("timeRes",&tRes,"timeRes/D");
  tree->Branch("timeRes",&vtRes);
  tree->Branch("posRad",&posRad,"posRad/D");
  tree->Branch("energy",&energy,"energy/D");
  tree->Branch("dirX",&dirX,"dirX/D");
  tree->Branch("dirY",&dirY,"dirY/D");
  tree->Branch("dirZ",&dirZ,"dirZ/D");
  tree->Branch("nhits",&nhits,"nhits/D");
  tree->Branch("sunDirX",&sunDirX,"sunDirX/D");
  tree->Branch("sunDirY",&sunDirY,"sunDirY/D");
  tree->Branch("sunDirZ",&sunDirZ,"sunDirZ/D");
  tree->Branch("cosThetaToSun",&cosThetaToSun,"cosThetaToSun/D");
  tree->Branch("runNumber",&runNumber,"runNumber/i");
  tree->Branch("eventGTID",&eventGTID,"eventGTID/i");
  tree->Branch("day",&day,"day/i");
  tree->Branch("sec",&sec,"sec/i");
  tree->Branch("fecd",&fecd,"fecd/i");
  tree->Branch("nsecs",&nsecs,"nsecs/D");
  tree->Branch("itr",&itr,"itr/D");
  tree->Branch("beta14",&iso,"iso/D");
  tree->Branch("theta_ij",&thij,"thij/D");
  char filenames[1500];
  Int_t i=0;
  Int_t TriggerType ;
  Int_t nentries=0;
  ULong64_t bitword ;
  while(in0>>filenames){
   cout<<" filenames "<<filenames<<endl ;
    // Load the RAT file
    RAT::DU::DSReader dsReader(filenames);
    
    // To later get the run info
    const RAT::DS::Run& run = dsReader.GetRun();
    
    size_t numEvents = dsReader.GetEntryCount();
    
    // Loop over all of the events
    for(size_t iEntry = 0; iEntry < numEvents; iEntry++ )
    {
	//posX=-10000, posY=-10000, posZ=-10000, posRad=-10000, time=99999, energy=0, ITR=-2, beta14=-10, dirX=-2.0, dirY=-2.0, dirZ=-2.0, nhits=0, sunDirX=0.0,sunDirY=0.0, sunDirZ=0.0, cosThetaToSun=-2.0, itr=-2, iso=-10, tRes=500 ;
	//posTheta = -10;  posPhi = -10;
        //day=0, sec=0, runNumber=0, eventGTID=0;
        //nsecs=0 ;	
	const RAT::DS::Entry& rDS = dsReader.GetEntry( iEntry ); // ds for each event
   	// const RAT::DS::MC& mc = rDS.GetMC(); // MC branch for each event 
	
	size_t iEVCount = rDS.GetEVCount(); // Number of triggered events in this MC event
	
	ULong64_t dcAnalysisWord = RAT::GetDataCleaningWord( "analysis_mask" );

	if (iEntry%5000 == 0) {
    	  cout<<iEntry<<endl;
	}
	
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
		  //cout<< " iEntry "<<iEntry<<" GTID "<<eventGTID<<" Applied bits "<< dcflags.GetApplied( dcflags.GetLatestPass() ).ToString() <<" bits flagged "<<dcflags.GetFlags( dcflags.GetLatestPass() ).ToString()<<endl ;
		  
		  //if bitword==dcflags.GetFlags( dcflags.GetLatestPass() ).ToString())
		  // Fill the nHit histogram, all events should have nHit
		  Float_t nHitOfEvent = ev.GetNhitsCleaned();
		  //nHit->Fill(nHitOfEvent);
		  nhits= nHitOfEvent ;
		  
		  TriggerType = ev.GetTrigType() ;
		  //cout<<" TriggerType "<<TriggerType<<endl ;
		  
                  vtRes->clear();
		  // Now get the radius, if applicable
		  Double_t radiusOfEvent = -1.;
		  TVector3 pos_fit;
		  bool directionFlag = false;
		  TVector3 u_fit;
		  const string fitName = "waterFitter";

                  for(unsigned int ipmt=0;ipmt<calpmts.GetFECDCount();ipmt++)//do FECD cuts
                  {
                    fecd = calpmts.GetFECDPMT(ipmt).GetID();
                  }

		  if(((TriggerType&2)==2))// || ((TriggerType&18)==18))
		  {
		    if( ev.FitResultExists(fitName) ) {  // It needs to exist 
		      if( ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() ) { // Needs to contain position
                        if( ev.GetFitResult(fitName).GetVertex(0).ValidPosition() ) { // Needs a valid position
                           RAT::DS::FitResult fResult = ev.GetFitResult(fitName);     // Get fit result
                           if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
                             RAT::DS::FitVertex fVertex = fResult.GetVertex(0);          // Get first fit vertex
                             pos_fit = fVertex.GetPosition();
                             time = fVertex.GetTime();
                             radiusOfEvent = pos_fit.Mag();
                             posX = pos_fit.X();
                             posY = pos_fit.Y();
                             posZ = pos_fit.Z();
                             posTheta =  pos_fit.CosTheta(); posPhi =  pos_fit.Phi();
                             posRad = radiusOfEvent ;
                             radius->Fill(radiusOfEvent);
                           }
                         } // end of loop if valid position
                       }// end of loop if contains position
			
			if( ev.GetFitResult(fitName).GetVertex(0).ContainsDirection() ) { // Needs to contain direction
			  if( ev.GetFitResult(fitName).GetVertex(0).ValidDirection() ) { // Needs a valid direction
			    
			    RAT::DS::FitResult fResult = ev.GetFitResult(fitName);     // Get fit result
			    if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
			      RAT::DS::FitVertex fVertex = fResult.GetVertex(0);          // Get first fit vertex
			      u_fit = fVertex.GetDirection();
			      dirX = u_fit.X();
			      dirY = u_fit.Y();
			      dirZ = u_fit.Z();
			      directionFlag = true;
			    }
			  }
			}
			
			if( ev.GetFitResult(fitName).GetVertex(0).ContainsEnergy() ) {
			  RAT::DS::FitResult fResult = ev.GetFitResult(fitName);
			  if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
			    energy= fResult.GetVertex(0).GetEnergy(); }
			}

                        //Calculate time residue
		        if( ev.GetFitResult(fitName).GetVertex(0).ContainsPosition() ) { // Needs to contain position
                          if( ev.GetFitResult(fitName).GetVertex(0).ValidPosition() ) {	
                            RAT::DS::FitResult fResult = ev.GetFitResult(fitName);     // Get fit result
                            if(fResult.GetValid()) { //cout<< " Fit is Valid "<<endl ;
                              RAT::DS::FitVertex fVertex = fResult.GetVertex(0); 
                              for(unsigned int ipmt=0;ipmt<calpmts.GetCount();ipmt++)
                              {
                                 TVector3 pmtpos = pmtInfo.GetPosition(calpmts.GetPMT(ipmt).GetID());
                                 double hitTime =(calpmts.GetPMT(ipmt)).GetTime();
                                 vtRes->push_back((calpmts.GetPMT(ipmt)).GetTime()-fVertex.GetTime()-(pmtpos-pos_fit).Mag()/grVelocity);
                              }
                            }
                          }
                        }

			if( !ev.ClassifierResultExists("ITR:waterFitter") ) continue;
			//try{ ( ev.GetClassifierResult( "ITR:waterFitter" ).GetValid()== false ) ; continue; }
			//catch(RAT::DS::ClassifierResult::NoClassifierResultError& e){cout << "No 'ITR' classifier invoked in the macro." << endl; return;}
			//try { itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );}
			//catch(RAT::DS::ClassifierResult::NoClassificationError& e){ cout<< " ITR classifier has failed "<< endl ; return ;}
			if ( !ev.GetClassifierResult( "ITR:waterFitter" ).GetValid() ) continue;
			itr = ev.GetClassifierResult( "ITR:waterFitter" ).GetClassification( "ITR" );		    
			
			if( !ev.ClassifierResultExists("isotropy:waterFitter") ) continue;
			//try{ (ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid()== false ) ; continue;}
			//catch(RAT::DS::ClassifierResult::NoClassifierResultError& e){cout << "No 'Beta14' classifier invoked in the macro." << endl; return;}
			//try { iso = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "snobeta14" );}
			//catch(RAT::DS::ClassifierResult::NoClassificationError& e){ cout<< " Beta14 classifier has failed "<< endl ; return ;}
			if ( !ev.GetClassifierResult( "isotropy:waterFitter" ).GetValid() ) continue;
			iso = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "snobeta14" );
                        thij = ev.GetClassifierResult( "isotropy:waterFitter" ).GetClassification( "thetaij" );
			//cout <<" itr "<<itr<<" iso "<<iso<<endl;
			
		        } // end of fit results
		      
                        tree->Fill(); 

		      // Fill the nHit histogram if the event was inside the AV
		      if (radiusOfEvent >= 0. && radiusOfEvent <= 6000.) {
			nHitAV->Fill(nHitOfEvent);
		      }
		    } //trigger cuts
		} // end of if loop for triggered word  
	    } // End loop over triggered events
	}
      } // End ds loop  
  } // end of the while loop
  // Give everything an axis label, as should always be
  radius->SetXTitle( "Radius [mm]" );
  radius->SetYTitle( "Counts" );
  nHit->SetXTitle( "nHit" );
  nHit->SetYTitle( "Counts" );  
  nHitAV->SetXTitle( "nHit" );
  nHitAV->SetYTitle( "Counts" );  
  nHit_toSun->SetXTitle( "nHit" );
  nHit_toSun->SetYTitle( "cos(#theta) to Sun" );
  nHit_toSun->SetZTitle( "Counts" );
  nHit_toSunAV->SetXTitle( "nHit" );
  nHit_toSunAV->SetYTitle( "cos(#theta) to Sun" );
  nHit_toSunAV->SetZTitle( "Counts" );
  nHit_toSun_radius->SetXTitle( "nHit" );
  nHit_toSun_radius->SetYTitle( "cos(#theta) to Sun" );
  nHit_toSun_radius->SetZTitle( "Radius [mm]" );
  
  // Write the histograms to fileName
  TString newname = "RatProcess_";
  TString processname = newname+filenames; 
  TFile *file=new TFile(processname,"RECREATE");
  file->cd();
  tree->Write();
  radius->Write();
  nHit->Write();
  nHitAV->Write();
  nHit_toSun->Write();
  nHit_toSunAV->Write();
  nHit_toSun_radius->Write();
  candidateEventInfo->Write();
  
}
