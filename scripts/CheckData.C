#include <iostream>

#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TRandom.h"

#include "TreeTrack.hh"

void CheckData(Double_t resol = 0.)
{
	// Global
	gROOT  -> Reset();
	gROOT  -> LoadMacro("TreeTrack.cc++");
	gStyle -> SetOptStat(0);

	// Read file
	TString fileName = "file_name_is_here";
	TString dataName = fileName + ".root";
	std::cout << "\033[1;33m[Notice] " << dataName.Data() << " is going to be processed.\033[0m" << std::endl;
	TFile* data = new TFile(dataName.Data());
	if ( !data -> IsOpen() )
	{
		std::cout << "\033[1;31m[Error] Cannot open file " << dataName.Data();
		std::cout << "! Check file name or path.\033[0m" << std::endl;

		return;
	}

	// Define tree
	TTree* tree = (TTree*) data -> Get("data");

	// Set branch addresses
	TClonesArray* tracks = new TClonesArray("TreeTrack", 100);
	Int_t    pol;
	Int_t    runNo;
	UInt_t   time;
	Double_t timeC;
	Int_t    trg;
	Int_t    nTrack;
	tree -> SetBranchAddress("polar"     , &pol   );
	tree -> SetBranchAddress("runnb"     , &runNo );
	tree -> SetBranchAddress("time_h"    , &time  );
	tree -> SetBranchAddress("time_cycle", &timeC );
	tree -> SetBranchAddress("trigger"   , &trg   );
	tree -> SetBranchAddress("N"         , &nTrack);
	tree -> SetBranchAddress("Tracks"    , &tracks);

	// Setting up spectra
	const Int_t nL        =   8; // Number of layers
	const Int_t nBinADC   = 600; // Number of bins for ADC axis
	const Int_t nBinTheta = 100; // Number of bins for angle axis
	const TString  lName[nL+1]  = {"FWC1", "FWC2", "FTH", "FRH1", "FRH2", "FRH3", "FRH4", "FRH5", "All"};
	const Double_t lRange[nL+1] = {    6.,     6.,   10.,   200.,   200.,   200.,   200.,   200.,  500.};
	TH1D* hADC[nL+1];
	TH1D* hADC_El[nL][48];
	TH2D* hADC_ADC[nL-1];
	TH2D* hADC_Theta[nL+1];
	TH1D* hTheta  = new TH1D("Theta", "Hit Distribution", nBinTheta, 0, 20);
	TH1D* hXS     = new TH1D("XS"   , "Cross Section"   , nBinTheta, 0, 20);
	TString tmp;
	for ( Int_t l = 0; l < nL + 1; l++ )
	{
		tmp.Form("ADC_%s", lName[l].Data());
		hADC[l] = new TH1D(tmp, tmp, nBinADC, 0, lRange[l]);
		hADC[l] -> SetFillColor(kGreen);
		for ( Int_t e = 0; e < 48; e++ )
		{
			if ( l == nL           ) break;
			if ( l != 2 && e == 24 ) break;
			tmp.Form("ADC_%s_El%i", lName[l].Data(), e + 1);
			hADC_El[l][e] = new TH1D(tmp, tmp, nBinADC, 0, lRange[l]);
			hADC_El[l][e] -> SetFillColor(kGreen);
		}
		
		tmp.Form("ADC_%svsTheta", lName[l].Data());
		hADC_Theta[l] = new TH2D(tmp, tmp, nBinTheta, 0, 20, nBinADC, 0, lRange[l]);
		hADC_Theta[l] -> GetXaxis() -> SetTitle("#theta_{Lab} [deg]");
		hADC_Theta[l] -> GetYaxis() -> SetTitle(lName[l]);
		hADC_Theta[l] -> SetOption("colz");
		if ( l < nL - 1 )
		{
			tmp.Form("ADC_%svs%s", lName[l].Data(), lName[l+1].Data() );
			hADC_ADC[l] = new TH2D(tmp, tmp, nBinADC, 0, lRange[l+1], nBinADC, 0, lRange[l]);
			hADC_ADC[l] -> GetXaxis() -> SetTitle(lName[l+1]);
			hADC_ADC[l] -> GetYaxis() -> SetTitle(lName[l]);
			hADC_ADC[l] -> GetYaxis() -> SetTitleOffset(1.5);
			hADC_ADC[l] -> SetOption("colz");
		}
	}

	// Setting up asymmetry
	TH1D* hNLeft  = new TH1D("hNLeft" , "Number of Hits (Left)" , nBinTheta, 0, 20);
	TH1D* hNRight = new TH1D("hNRight", "Number of Hits (Right)", nBinTheta, 0, 20);
	TH1D* hAsym   = new TH1D("hAsym"  , "Asymmetry"             , nBinTheta, 0, 20);
	Double_t phi;
	Double_t theta;

	// Setting up 2D hit map
	TH2D* hHitMap = new TH2D("hHitMap", "2D Hit Map", 1000, -0.5, 0.5, 1000, -0.5, 0.5);

	// Style
	hNLeft  -> SetLineColor(kBlack);
	hNRight -> SetLineColor(kBlack);
	hAsym   -> SetLineColor(kBlue );
	hTheta  -> SetLineColor(kBlack);
	hXS     -> SetLineColor(kBlack);
	hNLeft  -> GetXaxis() -> SetTitle("#theta_{Lab} [deg]");
	hNRight -> GetXaxis() -> SetTitle("#theta_{Lab} [deg]");
	hAsym   -> GetXaxis() -> SetTitle("#theta_{Lab} [deg]");
	hTheta  -> GetXaxis() -> SetTitle("#theta_{Lab} [deg]");
	hXS     -> GetXaxis() -> SetTitle("#theta_{Lab} [deg]");
	hNLeft  -> SetOption("E");
	hNRight -> SetOption("E");
	hAsym   -> SetOption("E");
	hHitMap -> SetOption("colz");

	// Looping over entries
	Int_t nEntries = tree -> GetEntries();
	std::cout << "\033[1;32m[Notice] " << nEntries << " entries will be processed.\033[0m" << std::endl;

	TreeTrack* tr;
	for ( Int_t i = 0; i < nEntries; i++ )
	{
		tree -> GetEntry(i);

		// Ignore empty event
		if( nTrack < 1 ) continue;

		// Loop over tracks
		for ( Int_t j = 0; j < nTrack; j++ )
		{
			tr = (TreeTrack*) tracks -> At(j);
			if ( !tr )
			{
				std::cout << "[Warning] At j = " << j
				          << ", there is no object although there should be one." << std::endl;
			}

			// Loop over hodoscope layers
			Double_t sumADC  = 0;
			Double_t nowADC  = 0;
			Double_t nextADC = 0;
			Int_t    nowEl   = 0;
			for ( Int_t l = 0; l < nL; l++ )
			{
				nowADC = tr -> fAdc[l];
				nowEl  = tr -> fEl[l];
				nowADC += gRandom -> Gaus(0., resol*nowADC/TMath::Sqrt(8.*TMath::Log(2.)));
				if ( l < nL - 1 )
				{
					nextADC = tr -> fAdc[l+1];
					nextADC += gRandom -> Gaus(0., resol*nextADC/TMath::Sqrt(8.*TMath::Log(2.)));
				}
				sumADC += nowADC;
				hADC[l] -> Fill(nowADC);
				if ( nowEl != -1 ) hADC_El[l][nowEl - 1] -> Fill(nowADC);
				if ( l < nL - 1 ) hADC_ADC[l] -> Fill(nextADC, nowADC);
				hADC_Theta[l] -> Fill(tr -> fTheta*TMath::RadToDeg(), nowADC);
			}
			hADC[nL]       -> Fill(sumADC);
			hADC_Theta[nL] -> Fill(tr -> fTheta*TMath::RadToDeg(), sumADC);
//
			// Fill counts
			phi   = tr -> fPhi   * TMath::RadToDeg();
			theta = tr -> fTheta * TMath::RadToDeg();
			if      ( (phi >  -45.0) && (phi <  45.0) ) hNLeft  -> Fill(theta);
			else if ( (phi < -135.0) || (phi > 135.0) ) hNRight -> Fill(theta);
//
			// Fill theta distibution
			hTheta -> Fill(theta);
//
			// Fill 2D hit map
			theta *= TMath::DegToRad();
			phi   *= TMath::DegToRad();
			hHitMap -> Fill(1.1765*TMath::Tan(theta)*TMath::Cos(phi), 1.1765*TMath::Tan(theta)*TMath::Sin(phi));
		}

		if ( (i+1) % 100000 == 0 ) std::cout << "[Notice] "<< i+1 << "-th event processing..." << std::endl;
	}

	// Caculate Asym
	Double_t countL, countR;
	Double_t asym, sigma;
	for ( Int_t i = 1; i <= hAsym -> GetNbinsX(); i++ )
	{
		countL = hNLeft  -> GetBinContent(i);
		countR = hNRight -> GetBinContent(i);
		if ( countL + countR <= 0 ) continue;

		asym  = (countL - countR) / (countL + countR);
		sigma = TMath::Sqrt((4*countL*countR) / TMath::Power(countL+countR, 3));

		hAsym -> SetBinContent(i, asym);
		hAsym -> SetBinError(i, sigma);
	}

	// Calculate cross section
	Double_t nHits, eNHits;
	Double_t phiDistCor;
	for ( Int_t i = 1; i <= hXS -> GetNbinsX(); i++ )
	{
		 nHits = hTheta -> GetBinContent(i);
		eNHits = hTheta -> GetBinError(i);

		if ( nHits == 0 ) continue;

		phiDistCor = TMath::Sin(hTheta -> GetBinCenter(hXS -> GetNbinsX()) * TMath::DegToRad())
		           / TMath::Sin(hTheta -> GetBinCenter(i)                  * TMath::DegToRad());

		hXS -> SetBinContent(i, phiDistCor * nHits);
		hXS -> SetBinError(i, phiDistCor * eNHits);
	}

	// Draw ADC
	TCanvas* cADC = new TCanvas("cADC", "ADC", 900, 900);
	cADC -> Divide(3, 3);
	for ( Int_t l = 0; l < nL + 1; l++ )
	{
		cADC -> cd(l+1);
		gPad -> SetLogz(1);
		hADC[l] -> Draw();
	}

	// Draw ADC_El
	TCanvas* cADCEl[nL];
	for ( Int_t l = 0; l < nL; l++ )
	{
		tmp.Form("ADC_El_%s", lName[l].Data());
		cADCEl[l] = new TCanvas(tmp.Data(), tmp.Data(), 900, 900);
		if ( l != 2 )
		{
			cADCEl[l] -> Divide(5, 5);
			for ( int e = 0; e < 24; e++ )
			{
				cADCEl[l] -> cd(e + 1) -> SetGrid();
				hADC_El[l][e] -> Draw();
			}
		}
		else
		{
			cADCEl[l] -> Divide(7, 7);
			for ( int e = 0; e < 48; e++ )
			{
				cADCEl[l] -> cd(e + 1) -> SetGrid();
				hADC_El[l][e] -> Draw();
			}
		}
	}

	// Draw ADC vs ADC
	TCanvas* cADC2D = new TCanvas("cADC2D", "ADC 2D", 900, 900);
	cADC2D -> Divide(3, 3);
	for ( Int_t l = 1; l < nL; l++ )
	{
		cADC2D -> cd(l);
		gPad -> SetLogz(1);
		hADC_ADC[l-1] -> Draw();
	}

	// Draw ADC vs Theta
	TCanvas* cADC1D = new TCanvas("cADC1D", "ADC theta", 900, 900);
	cADC1D -> Divide(3, 3);
	for ( Int_t l = 0; l < nL + 1; l++ )
	{
		cADC1D -> cd(l+1);
		gPad -> SetLogz(1);
		hADC_Theta[l] -> Draw();
	}

	// Draw counts
	TCanvas* cCount = new TCanvas("cCount", "Number of Hit", 1600, 800);
	cCount -> Divide(2, 1);
	cCount -> cd(1) -> SetGrid();
//	cCount -> cd(1) -> SetLogy();
	hNLeft -> Draw();
	cCount -> cd(2) -> SetGrid();
//	cCount -> cd(4) -> SetLogy();
	hNRight -> Draw();

	// Draw theta distribution and XS
	TCanvas* cTheta = new TCanvas("cTheta", "Theta Distribution", 1600, 800);
	cTheta -> Divide(2, 1);
	cTheta -> cd(1) -> SetGrid();
	hTheta -> Draw();
	cTheta -> cd(2) -> SetGrid();
	cTheta -> cd(2) -> SetLogy();
	hXS    -> Draw();

	// Draw Asym, CR and AP
	TCanvas* cA = new TCanvas("cA", "Asym", 900, 900);
	cA -> SetGrid();
	hAsym -> Draw();

	// Draw 2D hit map
	TCanvas* cHM = new TCanvas("cHM", "2D Hit Map", 900, 900);
	cHM -> SetGrid();
	hHitMap -> Draw("colz");

	// Save
	TString saveName;
	TString saveNameStart;
	TString saveNameEnd;
	saveName      = fileName + "_Draw.pdf";
	saveNameStart = saveName + "(";
	saveNameEnd   = saveName + ")";
	cADC      -> Print(saveNameStart.Data(), "pdf");
	cADCEl[0] -> Print(saveName.Data()     , "pdf");
	cADCEl[1] -> Print(saveName.Data()     , "pdf");
	cADCEl[2] -> Print(saveName.Data()     , "pdf");
	cADCEl[3] -> Print(saveName.Data()     , "pdf");
	cADCEl[4] -> Print(saveName.Data()     , "pdf");
	cADCEl[5] -> Print(saveName.Data()     , "pdf");
	cADCEl[6] -> Print(saveName.Data()     , "pdf");
	cADCEl[7] -> Print(saveName.Data()     , "pdf");
	cADC2D    -> Print(saveName.Data()     , "pdf");
	cADC1D    -> Print(saveName.Data()     , "pdf");
	cCount    -> Print(saveName.Data()     , "pdf");
	cTheta    -> Print(saveName.Data()     , "pdf");
	cA        -> Print(saveName.Data()     , "pdf");
	cHM       -> Print(saveNameEnd.Data()  , "pdf");
	TString outputName = fileName + "_Draw.root";
	TFile* result = new TFile(outputName.Data(), "RECREATE");
	cADC      -> Write();
	cADCEl[0] -> Write();
	cADCEl[1] -> Write();
	cADCEl[2] -> Write();
	cADCEl[3] -> Write();
	cADCEl[4] -> Write();
	cADCEl[5] -> Write();
	cADCEl[6] -> Write();
	cADCEl[7] -> Write();
	cADC2D    -> Write();
	cADC1D    -> Write();
	cCount    -> Write();
	cTheta    -> Write();
	cA        -> Write();
	cHM       -> Write();
	for ( Int_t l = 0; l < nL + 1; l++ )
	{
		hADC[l] -> Write();
		for ( Int_t e = 0; e < 48; e++ )
		{
			if ( l == nL           ) break;
			if ( l != 2 && e == 24 ) break;
			hADC_El[l][e] -> Write();
		}
		hADC_Theta[l] -> Write();
		if ( l < nL - 1 ) hADC_ADC[l] -> Write();
	}
	hTheta  -> Write();
	hNLeft  -> Write(); 
	hNRight -> Write(); 
	hAsym   -> Write(); 
	hHitMap -> Write();
	result  -> Close();

	std::cout << "\033[1;33m[Notice] Run for " << dataName.Data()
	          << " finished successfuly.\033[0m" << std::endl;

	// Delete
	delete result;
	delete cHM;
	delete cA;
	delete cTheta;
	delete cCount;
	delete cADC1D;
	delete cADC2D;
	for ( Int_t l = 0; l < nL; l++ ) delete cADCEl[l];
	delete cADC;
	delete hHitMap;
	delete hAsym;
	delete hNRight;
	delete hNLeft;
	for ( Int_t l = 0; l < nL + 1; l++ )
	{
		for ( Int_t e = 0; e < 48; e++ )
		{
			if ( l == nL           ) break;
			if ( l != 2 && e == 24 ) break;
			delete hADC_El[l][e];
		}
		delete hADC_Theta[l];
		if ( l < nL - 1 ) delete hADC_ADC[l];
	}
	delete hXS;
	delete hTheta;
	delete tracks;
	delete data;
}
