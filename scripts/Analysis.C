#include <iostream>
#include <algorithm>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1F.h"
#include "TMath.h"
#include "TClonesArray.h"

#include "TreeTrack.hh"

void Analysis()
{
	// TreeTrack library
	gROOT -> LoadMacro("TreeTrack.cc++");

	// Load data
	TString fileName = "file_name_is_here";
	TString dataName = "../" + fileName + ".root";
	std::cout << "\033[1;33m[Notice] " << dataName.Data() << " is going to be processed.\033[0m" << std::endl;
	TFile* data = new TFile(dataName.Data());
	if ( !data -> IsOpen() )
	{
		std::cout << "\033[1;31m[Error] Fail to open file\033[0m" << std::endl;
		return;
	}
	TTree* tree = (TTree*) data -> Get("data");
	if ( !tree )
	{
		std::cout << "\033[1;31m[Error] Fail to get tree\033[0m" << std::endl;
		return;
	}
	Double_t theta;
	Double_t phi;
	Double_t eDep[216];
	tree -> SetBranchAddress("theta", &theta);
	tree -> SetBranchAddress("phi"  , &phi  );
	TString detID;
	for ( Int_t i = 0; i < 216; i++ )
	{
		detID.Form("DetID%d", i);
		tree -> SetBranchAddress(detID.Data(), &eDep[i]);
	}

	// Prepare tree output
	Double_t beame     = 0;
	Int_t    target    = 0;
	Int_t    polar     = 0;
	Int_t    runnb     = 0;
	UInt_t   timeH     = 0;
	Double_t timeCycle = 0;
	Int_t    trigger   = 0;
	Int_t    nTracks   = 0;
	TClonesArray* trackArray = new TClonesArray("TreeTrack", 100);

	TString outputName = fileName + "_tree.root";
	TFile* outputFile = new TFile(outputName.Data(), "RECREATE");
	TTree* outputTree = new TTree("data", "data");
	outputTree -> Branch("beame"     , &beame         );
	outputTree -> Branch("target"    , &target        );
	outputTree -> Branch("polar"     , &polar         );
	outputTree -> Branch("runnb"     , &runnb         );
	outputTree -> Branch("time_h"    , &timeH         );
	outputTree -> Branch("time_cycle", &timeCycle     );
	outputTree -> Branch("trigger"   , &trigger       );
	outputTree -> Branch("Tracks"    , &trackArray    );
	outputTree -> Branch("N"         , &nTracks, "N/I");

	// Hit banks
	std::map<Int_t, Double_t> FWC1HitBank;
	std::map<Int_t, Double_t> FWC2HitBank;
	std::map<Int_t, Double_t> FTHHitBank;
	std::map<Int_t, Double_t> FRHHitBank[5];

	// Track reconstruction materials
	Double_t ADC[8];
	Int_t    El[8];
	TH1F* hFailReas = new TH1F("hFailReas", "Fail Reason"     , 11, -0.5, 10.5); 
	TH1F* hNTracks  = new TH1F("hNTracks" , "Number of Tracks", 11, -0.5, 10.5);
	std::map<Int_t, Double_t>::iterator itBank0;
	std::map<Int_t, Double_t>::iterator itBank1;
	std::map<Int_t, Double_t>::iterator itBank2;
	Int_t El1, El2, iFPC;
	TreeTrack* tTrack;

	// Loop
	Int_t nEntries = tree -> GetEntries();
	std::cout << "\033[1;32m[Notice] " << nEntries << " entries will be processed.\033[0m" << std::endl;
	for ( Int_t i = 0; i < nEntries; i++ )
	{
		if ( (i+1) % 10000 == 0 ) std::cout << i+1 << "-th entry processing..." << std::endl;

		tree -> GetEntry(i);

		// Initialize
		FWC1HitBank.clear();
		FWC2HitBank.clear();
		FTHHitBank.clear();
		for ( Int_t j = 0; j < 5; j++ ) FRHHitBank[j].clear();
		for ( Int_t j = 0; j < 8; j++ )
		{
			ADC[j] = -1;
			El[j]  = -1;
		}
		nTracks = 0;

		// Fill banks
		for ( Int_t j = 0; j < 216; j++ )
		{
			if ( eDep[j] <= 0 ) continue;

			if      ( j <  24 ) FWC1HitBank  .insert(std::pair<Int_t, Double_t>(j +   1, eDep[j])); // FWC1
			else if ( j <  48 ) FWC2HitBank  .insert(std::pair<Int_t, Double_t>(j -  23, eDep[j])); // FWC2
			else if ( j <  96 ) FTHHitBank   .insert(std::pair<Int_t, Double_t>(j -  47, eDep[j])); // FTH
			else if ( j < 120 ) FRHHitBank[0].insert(std::pair<Int_t, Double_t>(j -  95, eDep[j])); // FRH1
			else if ( j < 144 ) FRHHitBank[1].insert(std::pair<Int_t, Double_t>(j - 119, eDep[j])); // FRH2
			else if ( j < 168 ) FRHHitBank[2].insert(std::pair<Int_t, Double_t>(j - 143, eDep[j])); // FRH3
			else if ( j < 192 ) FRHHitBank[3].insert(std::pair<Int_t, Double_t>(j - 167, eDep[j])); // FRH4
			else                FRHHitBank[4].insert(std::pair<Int_t, Double_t>(j - 191, eDep[j])); // FRH5
		}

		// Loop over FTH hits
		if ( FTHHitBank.size() == 0 ) hFailReas -> Fill(0);
		for ( itBank0 = FTHHitBank.begin(); itBank0 != FTHHitBank.end(); itBank0++ )
		{
			El[2]  = itBank0 -> first;
			ADC[2] = itBank0 -> second;

			El1 = -1;
			El2 = -1;

			// Track check for FWC1
			if ( El[2] % 2 == 0 )
			{
				// Even
				if ( El[2] == 48 ) El2 = 1;
				else El2 = El[2]/2 + 1;
				El1 = El[2] / 2;
			}
			else
			{
				// Odd
				if ( El[2] == 1 ) El2 = 24;
				else El2 = (El[2]+1)/2 - 1;
				El1 = (El[2]+1) / 2;
			}

			itBank1 = FWC1HitBank.find(El1);
			itBank2 = FWC1HitBank.find(El2);
			if      ( itBank1 != FWC1HitBank.end() )
			{
				El[0]  = FWC1HitBank.find(El1) -> first;
				ADC[0] = FWC1HitBank.find(El1) -> second;
			}
			else if ( itBank2 != FWC1HitBank.end() )
			{
				El[0]  = FWC1HitBank.find(El2) -> first;
				ADC[0] = FWC1HitBank.find(El2) -> second;
			}
			else
			{
				hFailReas -> Fill(1);
				continue;
			}

			// Track check for FWC2
			if ( El[2] % 2 == 0 )
			{
				// Even
				if ( El[2] == 48 ) El1 = 1;
				else El1 = (El[2]+2) / 2;
				El2 = El[2] / 2;
			}
			else
			{
				// Odd
				if ( El[2] == 47 ) El2 = 1;
				else El2 = (El[2]+3) / 2;
				El1 = (El[2]+1) / 2;
			}

			itBank1 = FWC2HitBank.find(El1);
			itBank2 = FWC2HitBank.find(El2);
			if      ( itBank1 != FWC2HitBank.end() )
			{
				El[1]  = FWC2HitBank.find(El1) -> first;
				ADC[1] = FWC2HitBank.find(El1) -> second;
			}
			else if ( itBank2 != FWC2HitBank.end() )
			{
				El[1]  = FWC2HitBank.find(El2) -> first;
				ADC[1] = FWC2HitBank.find(El2) -> second;
			}
			else
			{
				hFailReas -> Fill(2);
				continue;
			}

			// Track check for FRH
			if ( El[2] % 2 == 0 )
			{
				// Even
				if ( El[2] == 48 ) El2 = 1;
				else El2 = El[2]/2 + 1;
				El1 = El[2] / 2;
			}
			else
			{
				// Odd
				if ( El[2] == 1 ) El2 = 24;
				else El2 = (El[2]+1)/2 - 1;
				El1 = (El[2]+1) / 2;
			}

			for ( Int_t k = 0; k < 5; k++ )
			{
				itBank1 = FRHHitBank[k].find(El1);
				itBank2 = FRHHitBank[k].find(El2);
				if      ( itBank1 != FRHHitBank[k].end() )
				{
					El[k+3]  = FRHHitBank[k].find(El1) -> first;
					ADC[k+3] = FRHHitBank[k].find(El1) -> second;
				}
				else if ( itBank2 != FRHHitBank[k].end() )
				{
					El[k+3]  = FRHHitBank[k].find(El2) -> first;
					ADC[k+3] = FRHHitBank[k].find(El2) -> second;
				}
				else
				{
					if ( k == 0 )
					{
						// Require plane 1 of FRH
						hFailReas -> Fill(3);
						continue;
					}
					else break;
				}
			}

			// Track check for FPC
			if   ( phi < 0. ) iFPC = TMath::Floor(48. * (phi+TMath::TwoPi()) / TMath::TwoPi()) + 1;
			else              iFPC = TMath::Floor(48. *  phi                 / TMath::TwoPi()) + 1;
			if ( iFPC != El[2] )
			{
				hFailReas -> Fill(4);
				continue;
			}

			tTrack = new ((*trackArray)[nTracks]) TreeTrack();
			tTrack -> SetTheta(theta);
			tTrack -> SetPhi(phi);
			tTrack -> SetElements(El);
//			std::cout << ADC[0] << ", " << ADC[1] << ", " << ADC[2] << ", " << ADC[3] << ", ";
//			std::cout << ADC[4] << ", " << ADC[5] << ", " << ADC[6] << ", " << ADC[7] << std::endl;
			tTrack -> SetADCs(ADC);
			nTracks++;
		}

		if ( nTracks )
		{
			hNTracks -> Fill(nTracks);
			outputTree -> Fill();
		}
	}

	outputTree -> Write();
	hFailReas  -> Write();
	hNTracks   -> Write();
	outputFile -> Close();
}
