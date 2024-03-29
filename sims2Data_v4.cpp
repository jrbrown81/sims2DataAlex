#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <TROOT.h>
#include <TFile.h>
#include <TLine.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TSystem.h"
#include <TLatex.h>
#include <TApplication.h>
#include <math.h>
#include <TSystemDirectory.h>
#include <TChain.h>
#include <TLegend.h>
#include <fstream>
#include <TGraph.h>
#include <list>
#include <vector>
#include <limits>
#include <TRandom3.h>
#include <dirent.h>
#include <TStopwatch.h>
#include <TMath.h>
#include <TTree.h>
#include <TString.h>

using std::setw;
using namespace std;
ifstream in_file;
ofstream logfile;

const Float_t G_FWHM = 2.35482;
const Float_t Epair = 4.6; // eV
const Float_t EpsilonR = 11;
const Double_t Epsilon0 = 8.854187817e-12;
const Double_t k = 1.3806488e-23;
const Float_t dd = 0.01; // m
const Float_t mu_e = 0.055; // m2*V/s
const Float_t HV = 1400;
const Float_t T = 303; //K
const Int_t nPtsH = 12;
const Float_t pixEdges[nPtsH]={-4.5,-3.6,-2.8,-2.0,-1.2,-0.4,0.4,1.2,2.0,2.8,3.6,4.5};
const TString pathRoot_in_calib = "/home/userfs/j/jb536/physhad/Medical/jb536/setup";
const Int_t nPixXY = 11;
const Int_t nAMs = 2;
const Float_t screen_width_def = 1600;
const Float_t screen_height_def = 900;

TString pathRoot = "./";

const Float_t Zoffset = 80;
const Int_t nPoints = 100;
// 500 points is used for the final best results, but it's rather slows down the program. Use 100 for good intermediate results
//const Int_t nPoints = 500;
Float_t factor2ConvertAnodeTime2Distance[2] = {95,127};
const Int_t prntOutF = 200000;

// Parameters to replicate experimental configuration
Float_t triggerThresh = 50; // keV
Bool_t errorFlags = 0;
Float_t maxDistDZ = 4.7; // mm, radius
Float_t fwhm[2] = {3.8,5.3}; //energy resolution of detectors (%)
Float_t dX_shift_AM0 = 0;
Float_t dY_shift_AM0 = 1.6;
Float_t dX_shift_AM1 = -1.6;
Float_t dY_shift_AM1 = -1.6;
Double_t D, N;

TH2F *plane[2];
TH2F *plane2[2];
TRandom3 *rand3;
TH2I *pixelPattern[2];
TH1F *usedPixelsH[2];
TString spectrumName = "Na-22, 1400V";
TH1F *dPhiAngle_Ewin;
TH1F *dPhiAngleNorm_Ewin;
TH1F *dPhiAngle1_Ewin;
TH1F *dPhiAngle1norm_Ewin;
const Float_t Ewindow4PhiAnalyis_min = 480;
const Float_t Ewindow4PhiAnalyis_max = 540;
const Int_t nBins_dPhi = 38;
const Float_t ThetaWindowFor4PhiAnalyis_min = 94;
const Float_t ThetaWindowFor4PhiAnalyis_max = 102;
const Float_t relativePhiAngle = 90;

TChain *fChain; 
TStopwatch localTimer;
Float_t totalTimeElapced = 0;
TString num_str;

Int_t    GetEntry(Long64_t entry);
Long64_t LoadTree(Long64_t entry);
Int_t findAM(Float_t z);
Int_t findPixelID(Int_t am, Float_t x, Float_t y);
Double_t sigmaTotal2_t_scaled(Double_t *, Double_t *);
Bool_t getPixel2Dcoordinates(const Int_t, const Int_t, Int_t&, Int_t&);
Float_t getPhiAngleDeg(const Float_t, const Float_t);

int main()
{
	rand3 = new TRandom3();
	rand3->SetSeed(0);
	//gStyle->SetOptStat(0);
	TH1::AddDirectory(0);
	gErrorIgnoreLevel = 1001;
	gStyle->SetPalette(1);
	D = k*T/TMath::Qe()*mu_e;

	// Declaration of leaf types
	Float_t         number;
	Float_t         theta1;
	Float_t         phi1;
	Float_t         E1_a;
	Float_t         X1_a;
	Float_t         Y1_a;
	Float_t         Z1_a;
	Float_t         E1_b;
	Float_t         X1_b;
	Float_t         Y1_b;
	Float_t         Z1_b;
	Float_t         theta2;
	Float_t         phi2;
	Float_t         E2_a;
	Float_t         X2_a;
	Float_t         Y2_a;
	Float_t         Z2_a;
	Float_t         E2_b;
	Float_t         X2_b;
	Float_t         Y2_b;
	Float_t         Z2_b;

	// List of branches
	TBranch        *b_number;   //!
	TBranch        *b_theta1;   //!
	TBranch        *b_phi1;   //!
	TBranch        *b_E1_a;   //!
	TBranch        *b_X1_a;   //!
	TBranch        *b_Y1_a;   //!
	TBranch        *b_Z1_a;   //!
	TBranch        *b_E1_b;   //!
	TBranch        *b_X1_b;   //!
	TBranch        *b_Y1_b;   //!
	TBranch        *b_Z1_b;   //!
	TBranch        *b_theta2;   //!
	TBranch        *b_phi2;   //!
	TBranch        *b_E2_a;   //!
	TBranch        *b_X2_a;   //!
	TBranch        *b_Y2_a;   //!
	TBranch        *b_Z2_a;   //!
	TBranch        *b_E2_b;   //!
	TBranch        *b_X2_b;   //!
	TBranch        *b_Y2_b;   //!
	TBranch        *b_Z2_b;   //!

	fChain = new TChain("tree");
	TString fname2do;
	Int_t Nfiles = 0;
	TSystemDirectory dir("./","./");
	TList *files = dir.GetListOfFiles();
	if (files)
	{
		TSystemFile *file;
		TIter next(files);
		cout << "Counting valid root files in " << pathRoot << " ... ";
		while ((file=(TSystemFile*)next()))
		{
			fname2do = file->GetName();
			if (file->IsDirectory()) files->Remove(file);
			if (!fname2do.EndsWith(".root")) files->Remove(file);
			if (!fname2do.Contains("ill")) files->Remove(file);
		}
		next.Reset();
		cout << "List of files to read:" << endl;
		while ((file=(TSystemFile*)next()))
		{
			fname2do = file->GetName();
			Nfiles++;
		}
		next.Reset();
		if (Nfiles < 1)
		{
			cerr << "No valid ROOT files to read. Exiting." << endl;
			return kFALSE;
		}
		cout << "Total of " << Nfiles << " files found." << endl;
		while ((file=(TSystemFile*)next()))
		{
			fname2do = pathRoot + file->GetName();
			cout << "Reading file:" << endl;
			cout << fname2do << endl;
			fChain->Add(fname2do);
		}
	}

	fChain->SetBranchAddress("number", &number, &b_number);
	fChain->SetBranchAddress("theta1", &theta1, &b_theta1);
	fChain->SetBranchAddress("phi1", &phi1, &b_phi1);
	fChain->SetBranchAddress("E1_a", &E1_a, &b_E1_a);
	fChain->SetBranchAddress("X1_a", &X1_a, &b_X1_a);
	fChain->SetBranchAddress("Y1_a", &Y1_a, &b_Y1_a);
	fChain->SetBranchAddress("Z1_a", &Z1_a, &b_Z1_a);
	fChain->SetBranchAddress("E1_b", &E1_b, &b_E1_b);
	fChain->SetBranchAddress("X1_b", &X1_b, &b_X1_b);
	fChain->SetBranchAddress("Y1_b", &Y1_b, &b_Y1_b);
	fChain->SetBranchAddress("Z1_b", &Z1_b, &b_Z1_b);
	fChain->SetBranchAddress("theta2", &theta2, &b_theta2);
	fChain->SetBranchAddress("phi2", &phi2, &b_phi2);
	fChain->SetBranchAddress("E2_a", &E2_a, &b_E2_a);
	fChain->SetBranchAddress("X2_a", &X2_a, &b_X2_a);
	fChain->SetBranchAddress("Y2_a", &Y2_a, &b_Y2_a);
	fChain->SetBranchAddress("Z2_a", &Z2_a, &b_Z2_a);
	fChain->SetBranchAddress("E2_b", &E2_b, &b_E2_b);
	fChain->SetBranchAddress("X2_b", &X2_b, &b_X2_b);
	fChain->SetBranchAddress("Y2_b", &Y2_b, &b_Y2_b);
	fChain->SetBranchAddress("Z2_b", &Z2_b, &b_Z2_b);
   
	usedPixelsH[0] = new TH1F("usedPixelsH_AM0","Used pixels, AM0",968,0,968);
	usedPixelsH[0]->GetXaxis()->SetTitleOffset(1.2);
	usedPixelsH[0]->GetXaxis()->SetTitle("Pixel number, AM0");

	usedPixelsH[1] = new TH1F("usedPixelsH_AM1","Used pixels, AM1",968,0,968);
	usedPixelsH[1]->GetXaxis()->SetTitleOffset(1.2);
	usedPixelsH[1]->GetXaxis()->SetTitle("Pixel number, AM1");

	for (Int_t im = 0; im < nAMs; im++)
	{
		plane[im] = new TH2F(Form("plane_AM%d",im),Form("My Image, AM%d",im),nPtsH-1,pixEdges,nPtsH-1,pixEdges);
		plane[im]->GetXaxis()->SetTitle("Pixel number");
		plane[im]->GetXaxis()->SetTitleOffset(1.2);
		plane[im]->GetYaxis()->SetTitle("Pixel number");
		plane[im]->GetYaxis()->SetTitleOffset(1.4);
		
		plane2[im] = new TH2F(Form("plane2_AM%d",im),Form("Jamie's Image, AM%d",im),nPtsH-1,pixEdges,nPtsH-1,pixEdges);
		plane2[im]->GetXaxis()->SetTitle("Pixel number");
		plane2[im]->GetXaxis()->SetTitleOffset(1.2);
		plane2[im]->GetYaxis()->SetTitle("Pixel number");
		plane2[im]->GetYaxis()->SetTitleOffset(1.4);
		
		pixelPattern[im] = new TH2I(Form("pixelPattern_AM%d",im),"",nPtsH-1,pixEdges,nPtsH-1,pixEdges);
		pixelPattern[im]->GetXaxis()->SetTitle("Pixel number");
		pixelPattern[im]->GetXaxis()->SetTitleOffset(1.2);
		pixelPattern[im]->GetYaxis()->SetTitle("Pixel number");
		pixelPattern[im]->GetYaxis()->SetTitleOffset(1.4);
		
		TString fname_pxl_cfg = Form("%s/LISA_QET_pixel_map_config_AM%d_GM1.txt",pathRoot_in_calib.Data(),im);
		cout << "fname_pxl_cfg  = " << fname_pxl_cfg << endl;
		if (gSystem->AccessPathName(fname_pxl_cfg))
		{
			cerr << "ERROR: Pixel config file " << fname_pxl_cfg << " doesn't exist. Exiting." << endl;
			return 0;
		}
		in_file.open(fname_pxl_cfg, ios::in);
		for (Int_t i = nPixXY; i >= 1 && !in_file.eof(); i--)
		{
			for (Int_t j = 0; j < nPixXY && !in_file.eof(); j++)
			{
				Int_t pix = -1;
				in_file >> pix;
				if (pix < 1 || pix > 968) continue;
				pixelPattern[im]->SetBinContent(j+1,i,pix);
			}
		}
		in_file.close();
		in_file.clear();
	}
	
	Long_t nEvents2Analyse = fChain->GetEntries();
	TString nev_str = Form("%ld",nEvents2Analyse);
	for (Int_t is = nev_str.Length(); is > 1; is--)
		if ((is-nev_str.Length()+2)%4 == 0) nev_str.Insert(is-1,",",1);
	cout << "nEvents2Analyse = " << nev_str << endl;
	
	cout << "Trigger threshold: " << triggerThresh << " keV" << endl;
	cout << "FWHM for energy smearing: " << fwhm[0] << "% and " << fwhm[1] << "%" << endl;
	if(errorFlags) cout << endl << "Error warnings inhibited" << endl;
    	
	dPhiAngle_Ewin = new TH1F("dPhiAngle_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngle_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	dPhiAngle_Ewin->Sumw2();
	
	dPhiAngleNorm_Ewin = new TH1F("dPhiAngleNorm_Ewin",Form("%s, #Delta#varphi within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi,-180,180);
	dPhiAngleNorm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngleNorm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	dPhiAngleNorm_Ewin->Sumw2();

	dPhiAngle1_Ewin = new TH1F("dPhiAngle1_Ewin",Form("%s, |#Delta#varphi| within  %.0f keV < E < %.0f keV",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi/2,0,180);
	dPhiAngle1_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	
	dPhiAngle1norm_Ewin = new TH1F("dPhiAngle1norm_Ewin",Form("%s, |#Delta#varphi| within  %.0f keV < E < %.0f keV, normalised at #Delta#varphi = 0",spectrumName.Data(),Ewindow4PhiAnalyis_min,Ewindow4PhiAnalyis_max),nBins_dPhi/2,0,180);
	dPhiAngle1norm_Ewin->GetXaxis()->SetTitleOffset(1.2);
	dPhiAngle1norm_Ewin->GetXaxis()->SetTitle("#Delta#varphi, degrees");		
	dPhiAngle1norm_Ewin->Sumw2();
  
	// My temporary data definations
	Int_t eventCounter=0;
	Int_t amID[2]={-1,-1};
	Int_t pixel[2][2]={{-1,-1},{-1,-1}};
	Int_t pixel2[2][2]={{-1,-1},{-1,-1}};
	Float_t energy[2][2]={{0,0},{0,0}};
	Int_t pos[2][2]={{0,0},{0,0}};
	//Int_t nTrigPixels[2]={0,0};
	Float_t xvec[2][2]={{0,0},{0,0}};
	Float_t yvec[2][2]={{0,0},{0,0}};
	Float_t zvec[2][2]={{0,0},{0,0}};

	// Definitions for new trees for output of exp. data format
	//Long_t event = -1;
	Int_t nAMs = 2;
	Long_t *event_e = new Long_t[nAMs];
	Int_t *nTrigPixels_e = new Int_t[nAMs];
	Int_t *AM_e = new Int_t[nAMs];
	Int_t *GM_e = new Int_t[nAMs];
	Int_t *nAMsInEvent_e = new Int_t[nAMs];
	Int_t *pixel_e = new Int_t[nAMs];
	Float_t *E = new Float_t[nAMs];
	Float_t *E_neg = new Float_t[nAMs];
	Int_t *time_e = new Int_t[nAMs];
	Int_t *triggerFlag_e = new Int_t[nAMs];
	Int_t *nloop_e = new Int_t[nAMs];
	Int_t *timeDetect_e = new Int_t[nAMs];
	Int_t *timeDetectPos_e = new Int_t[nAMs];
	Float_t *Temp_e = new Float_t[nAMs];
	Int_t *pos_e = new Int_t[nAMs];
	Long_t event_eev;
	Int_t *AM_flag_eev = new Int_t[nAMs];
	Int_t *cath_flag_eev = new Int_t[nAMs];
	Int_t nAMsInEvent_eev;
	Int_t *nTrigPixels_eev = new Int_t[nAMs];
	Int_t *Nvar2 = new Int_t[nAMs];
	Float_t *E0_e = new Float_t[nAMs];
	Float_t *E0_x = new Float_t[nAMs];
	Float_t *E0_y = new Float_t[nAMs];
	Float_t *E0_z = new Float_t[nAMs];
	Int_t *E0_pixel = new Int_t[nAMs];
	Float_t *E1_e = new Float_t[nAMs];
	Float_t *E1_x = new Float_t[nAMs];
	Float_t *E1_y = new Float_t[nAMs];
	Float_t *E1_z = new Float_t[nAMs];
	Int_t *E1_pixel = new Int_t[nAMs];
	Float_t *theta_sim = new Float_t[nAMs];
	Float_t *phi_sim = new Float_t[nAMs];

	TString rfileOut = Form("sortedData_Ecal_simulated_%dpt.root",nPoints);
	TFile *f2 = new TFile(rfileOut,"recreate");
	// Event tree
	TTree *tree2_ev = new TTree("tree2_events","Event tree of E calib data");
	TTree **tree2 = new TTree*[nAMs];
	tree2_ev->Branch("event",&event_eev,"event/L");
	tree2_ev->Branch("nAMsInEvent",&nAMsInEvent_eev,"nAMsInEvent/I");
	tree2_ev->Branch("nTrigPixels",nTrigPixels_eev,Form("nTrigPixels[%d]/I",nAMs));
	tree2_ev->Branch("AM_flag",AM_flag_eev,Form("AM_flag[%d]/I",nAMs));
	tree2_ev->Branch("cath_flag",cath_flag_eev,Form("cath_flag[%d]/I",nAMs));
	// AM trees
	for (Int_t im = 0; im < nAMs; im++)
	{
		tree2[im] = new TTree(Form("tree2_AM%d",im), Form("E calib data tree AM%d",im));
		tree2[im]->Branch("event",&event_e[im],"event/L");
		tree2[im]->Branch("nTrigPixels",&nTrigPixels_e[im],"nTrigPixels/I");
		tree2[im]->Branch("AM",&AM_e[im],"AM/I");
		tree2[im]->Branch("GM",&GM_e[im],"GM/I");
		tree2[im]->Branch("nAMsInEvent",&nAMsInEvent_e[im],"nAMsInEvent/I");
		tree2[im]->Branch("pixel",&pixel_e[im],"pixel/I");
		tree2[im]->Branch("E",&E[im],"E/F");
		tree2[im]->Branch("E_neg",&E_neg[im],"E_neg/F");
		tree2[im]->Branch("pos",&pos_e[im],"pos/I");
		tree2[im]->Branch("time",&time_e[im],"time/I");
		tree2[im]->Branch("triggerFlag",&triggerFlag_e[im],"triggerFlag/I");
		tree2[im]->Branch("nloop",&nloop_e[im],"nloop/I");
		tree2[im]->Branch("timeDetect",&timeDetect_e[im],"timeDetect/I");
		tree2[im]->Branch("timeDetectPos",&timeDetectPos_e[im],"timeDetectPos/I");
		tree2[im]->Branch("Temp",&Temp_e[im],"Temp/F");
		tree2[im]->Branch("E0_e",&E0_e[im],"E0_e/F");
		tree2[im]->Branch("E0_x",&E0_x[im],"E0_x/F");
		tree2[im]->Branch("E0_y",&E0_y[im],"E0_y/F");
		tree2[im]->Branch("E0_z",&E0_z[im],"E0_z/F");
		tree2[im]->Branch("E0_pixel",&E0_pixel[im],"E0_pixel/I");
		tree2[im]->Branch("E1_e",&E1_e[im],"E1_e/F");
		tree2[im]->Branch("E1_x",&E1_x[im],"E1_x/F");
		tree2[im]->Branch("E1_y",&E1_y[im],"E1_y/F");
		tree2[im]->Branch("E1_z",&E1_z[im],"E1_z/F");
		tree2[im]->Branch("E1_pixel",&E1_pixel[im],"E1_pixel/I");
		tree2[im]->Branch("theta_sim",&theta_sim[im],"theta_sim/F");
		tree2[im]->Branch("phi_sim",&phi_sim[im],"phi_sim/F");
		Nvar2[im] = tree2[im]->GetNbranches();
	}

	TF1 *sigTotal_t_scaled = new TF1("sigTotal_t_scaled",sigmaTotal2_t_scaled,0,dd*dd/mu_e/HV*1e6);
	
	// Main loop over input file
	Long64_t nbytes = 0, nb = 0;
	Long_t savedEvents = 0;
	localTimer.Start();
	//for (Long64_t jentry = 0; jentry < nEvents2Analyse; jentry++)
	for (Long64_t jentry = 0; jentry < 2000000; jentry++)
	{
		nb = fChain->GetEntry(jentry);
		nbytes += nb;
		eventCounter++;
		if (jentry%prntOutF == 0)
		{
			num_str = Form("%lld",jentry);
			for (Int_t is = num_str.Length(); is > 1; is--)
			{
				if ((is-num_str.Length()+2)%4 == 0) num_str.Insert(is-1,",",1);
			}
			cout << Form("%s, Nev = %s, RealTime: %.3fs, CPUTime: %.3fs",num_str.Data(),nev_str.Data(),localTimer.RealTime(), localTimer.CpuTime()) << endl;
			totalTimeElapced += localTimer.RealTime();
			localTimer.ResetRealTime();
			localTimer.ResetCpuTime();
			localTimer.Start();
		}
      
		//cout << event_eev << " " << nAMsInEvent_eev << " " << nTrigPixels_eev[0] << " " << nTrigPixels_eev[1] << " " << AM_flag_eev[0] << " " << AM_flag_eev[1] << " " << cath_flag_eev[0] << " " << cath_flag_eev[1]  << endl;

		// Find which AMs and pixels hit
		// Find AMs
		amID[0]=findAM(Z1_a);
		amID[1]=findAM(Z2_a);
		if(amID[0]==-1 || amID[1]==-1) 
		{
			if(errorFlags) cout << "Error: Gamma ray not in a detector!" << endl;
			continue;
		}
		if(amID[0]==amID[1])
		{
			if(errorFlags) cout << "Error: Same AM hit by both gammas!" << endl;
			continue;
		}
		for(int im = 0; im < nAMs; im++) plane[amID[im]]->Reset();
		for(int im = 0; im < nAMs; im++) plane2[amID[im]]->Reset();
		Int_t nTrigPixels[2]={0,0};

		//if (amID[0] == 0) cout << Form("AM%d: [ pixel %d, X1_a = %.2f, Y1_a = %.2f, E1_a = %.1f ] ; [ pixel %d, X1_b = %.2f, Y1_b = %.2f, E1_b = %.1f ]", amID[0], pixel[amID[0]][0], X1_a, Y1_a, E1_a, pixel[amID[0]][1], X1_b, Y1_b, E1_b) << endl;
		//if (amID[0] == 1) cout << Form("AM%d: [ pixel %d, X1_a = %.2f, Y1_a = %.2f, E1_a = %.1f ] ; [ pixel %d, X1_b = %.2f, Y1_b = %.2f, E1_b = %.1f ]", amID[0], pixel[amID[0]][0], -X1_a, Y1_a, E1_a, pixel[amID[0]][1], -X1_b, Y1_b, E1_b) << endl;
		//if (amID[0] == 0) cout << Form("AM%d: [ pixel %d, X2_a = %.2f, Y2_a = %.2f, E2_a = %.1f ] ; [ pixel %d, X2_b = %.2f, Y2_b = %.2f, E2_b = %.1f ]", amID[1], pixel[amID[1]][0], X2_a, Y2_a, E2_a, pixel[amID[1]][1], X2_b, Y2_b, E2_b) << endl;
		//if (amID[0] == 1) cout << Form("AM%d: [ pixel %d, X2_a = %.2f, Y2_a = %.2f, E2_a = %.1f ] ; [ pixel %d, X2_b = %.2f, Y2_b = %.2f, E2_b = %.1f ]", amID[1], pixel[amID[1]][0], -X2_a, Y2_a, E2_a, pixel[amID[1]][1], -X2_b, Y2_b, E2_b) << endl;
		
	    energy[amID[0]][0] = rand3->Gaus(E1_a,fwhm[amID[0]]*E1_a/G_FWHM/100);
		energy[amID[0]][1] = rand3->Gaus(E1_b,fwhm[amID[0]]*E1_b/G_FWHM/100);
		energy[amID[1]][0] = rand3->Gaus(E2_a,fwhm[amID[1]]*E2_a/G_FWHM/100);
		energy[amID[1]][1] = rand3->Gaus(E2_b,fwhm[amID[1]]*E2_b/G_FWHM/100);

		if (amID[0] == 0) plane2[amID[0]]->Fill(X1_a,Y1_a,energy[amID[0]][0]);
		if (amID[0] == 0) plane2[amID[0]]->Fill(X1_b,Y1_b,energy[amID[0]][1]);
		if (amID[1] == 0) plane2[amID[1]]->Fill(X2_a,Y2_a,energy[amID[1]][0]);
		if (amID[1] == 0) plane2[amID[1]]->Fill(X2_b,Y2_b,energy[amID[1]][1]);
		if (amID[0] == 1) plane2[amID[0]]->Fill(-X1_a,Y1_a,energy[amID[0]][0]);
		if (amID[0] == 1) plane2[amID[0]]->Fill(-X1_b,Y1_b,energy[amID[0]][1]);
		if (amID[1] == 1) plane2[amID[1]]->Fill(-X2_a,Y2_a,energy[amID[1]][0]);
		if (amID[1] == 1) plane2[amID[1]]->Fill(-X2_b,Y2_b,energy[amID[1]][1]);

		xvec[amID[0]][0] = X1_a;
		xvec[amID[0]][1] = X1_b;
		xvec[amID[1]][0] = X2_a;
		xvec[amID[1]][1] = X2_b;
		if (amID[0] == 1)
		{
			xvec[amID[0]][0] = -X1_a;
			xvec[amID[0]][1] = -X1_b;			
		}
		if (amID[1] == 1)
		{
			xvec[amID[1]][0] = -X2_a;
			xvec[amID[1]][1] = -X2_b;			
		}
		xvec[0][0] += dX_shift_AM0;
		xvec[0][1] += dX_shift_AM0;
		xvec[1][0] += dX_shift_AM1;
		xvec[1][1] += dX_shift_AM1;
		if(TMath::Abs(xvec[0][0]) > 4.5 || TMath::Abs(xvec[0][1]) > 4.5) continue;
		if(TMath::Abs(xvec[1][0]) > 4.5 || TMath::Abs(xvec[1][1]) > 4.5) continue;
		
		yvec[amID[0]][0] = Y1_a;
		yvec[amID[0]][1] = Y1_b;
		yvec[amID[1]][0] = Y2_a;
		yvec[amID[1]][1] = Y2_b;
		yvec[0][0] += dY_shift_AM0;
		yvec[0][1] += dY_shift_AM0;
		yvec[1][0] += dY_shift_AM1;
		yvec[1][1] += dY_shift_AM1;
		if(TMath::Abs(yvec[0][0]) > 4.5 || TMath::Abs(yvec[0][1]) > 4.5) continue;
		if(TMath::Abs(yvec[1][0]) > 4.5 || TMath::Abs(yvec[1][1]) > 4.5) continue;
		
		zvec[amID[0]][0] = Z1_a;
		zvec[amID[0]][1] = Z1_b;
		zvec[amID[1]][0] = Z2_a;
		zvec[amID[1]][1] = Z2_b;
				
	    pixel[amID[0]][0] = findPixelID(amID[0],xvec[amID[0]][0],yvec[amID[0]][0]);
	    pixel[amID[0]][1] = findPixelID(amID[0],xvec[amID[0]][1],yvec[amID[0]][1]);
	    pixel[amID[1]][0] = findPixelID(amID[1],xvec[amID[1]][0],yvec[amID[1]][0]);
	    pixel[amID[1]][1] = findPixelID(amID[1],xvec[amID[1]][1],yvec[amID[1]][1]);

		if(TMath::Hypot(xvec[0][0],yvec[0][0]) > maxDistDZ || TMath::Hypot(xvec[1][0],yvec[1][0]) > maxDistDZ) continue;
		
		for(int im = 0; im < nAMs; im++)
		{
			Double_t locZ = TMath::Abs(zvec[amID[im]][0]) - Zoffset;
			locZ = dd*1000 - locZ;
			Double_t driftT = factor2ConvertAnodeTime2Distance[amID[im]]*locZ/1000; // in us
			N = energy[amID[im]][0]*1000/Epair;
			Double_t RR1 = sigTotal_t_scaled->Eval(driftT)/1000;
			//cout << "driftT = " << driftT << ", RR1 = " << RR1 << ", E = " << energy[amID[im]][0] << ", locZ = " << locZ << ", amID = " << amID[im] << endl;
			locZ = TMath::Abs(zvec[amID[im]][1]) - Zoffset;
			locZ = dd*1000 - locZ;
			driftT = factor2ConvertAnodeTime2Distance[amID[im]]*locZ/1000; // in us
			N = energy[amID[im]][1]*1000/Epair;
			Double_t RR2 = sigTotal_t_scaled->Eval(driftT)/1000;
			//cout << "driftT = " << driftT << ", RR2 = " << RR2 << ", E = " << energy[amID[im]][1] << ", locZ = " << locZ << ", amID = " << amID[im] << endl;
			//cout << "-----------------------------------------" << endl;
			for (Int_t i = 0; i < nPoints; i++)
			{
				//Double_t xx = xvec[amID[im]][0];
				//Double_t yy = yvec[amID[im]][0];
				
				Double_t xx = 1e6, yy = 1e6, zz = 1e6;
				while (TMath::Sqrt((xx-xvec[amID[im]][0])*(xx-xvec[amID[im]][0])+(yy-yvec[amID[im]][0])*(yy-yvec[amID[im]][0])+zz*zz) > 3*RR1)
				{
					xx = rand3->Gaus(xvec[amID[im]][0],RR1);
					yy = rand3->Gaus(yvec[amID[im]][0],RR1);
					zz = rand3->Gaus(0,RR1);
				}
				plane[amID[im]]->Fill(xx,yy,energy[amID[im]][0]/nPoints);
				//cout << Form("AM%d: [ pixel %d, xx = %.2f, yy = %.2f, E = %.1f ] ;", amID[im], findPixelID(amID[im],xx,yy), xx, yy, energy[amID[im]][0]);
				
				//xx = xvec[amID[im]][1];
				//yy = yvec[amID[im]][1];
				xx = 1e6;
				yy = 1e6;
				zz = 1e6;
				while (TMath::Sqrt((xx-xvec[amID[im]][1])*(xx-xvec[amID[im]][1])+(yy-yvec[amID[im]][1])*(yy-yvec[amID[im]][1])+zz*zz) > 3*RR2)
				{
					xx = rand3->Gaus(xvec[amID[im]][1],RR2);
					yy = rand3->Gaus(yvec[amID[im]][1],RR2);
					zz = rand3->Gaus(0,RR2);
				}
				plane[amID[im]]->Fill(xx,yy,energy[amID[im]][1]/nPoints);
				//cout << Form(" [ pixel %d, xx = %.2f, yy = %.2f, E = %.1f ]", findPixelID(amID[im],xx,yy), xx, yy, energy[amID[im]][1]) << endl;
			}
			
			for (Int_t ii = 1; ii <= nPixXY; ii++)
			{
				for (Int_t jj = 1; jj <= nPixXY; jj++)
				{
					if (plane[amID[im]]->GetBinContent(ii,jj) >= triggerThresh) nTrigPixels[amID[im]]++;
				}
			}
		}
		
		if(nTrigPixels[0] == 0 || nTrigPixels[1] == 0) continue;
	
		if (energy[amID[0]][0] + energy[amID[0]][1] > Ewindow4PhiAnalyis_min && energy[amID[0]][0] + energy[amID[0]][1] < Ewindow4PhiAnalyis_max
				&& energy[amID[1]][0] + energy[amID[1]][1] > Ewindow4PhiAnalyis_min && energy[amID[1]][0] + energy[amID[1]][1] < Ewindow4PhiAnalyis_max)
		{
			Float_t theta1_loc = theta1*TMath::RadToDeg();
			Float_t theta2_loc = theta2*TMath::RadToDeg();
			if (ThetaWindowFor4PhiAnalyis_min > 90 && ThetaWindowFor4PhiAnalyis_max > 90)
			{
				if (theta1_loc < 90) theta1_loc = 180 - theta1_loc;
				if (theta2_loc < 90) theta2_loc = 180 - theta2_loc;
			}
			if (theta1_loc > ThetaWindowFor4PhiAnalyis_min && theta1_loc < ThetaWindowFor4PhiAnalyis_max
					&& theta2_loc > ThetaWindowFor4PhiAnalyis_min && theta2_loc < ThetaWindowFor4PhiAnalyis_max)
			{
				if (TMath::Hypot(xvec[0][0]-xvec[0][1],yvec[0][0]-yvec[0][1]) > 0.8*TMath::Sqrt(2) && TMath::Hypot(xvec[1][0]-xvec[1][1],yvec[1][0]-yvec[1][1]) > 0.8*TMath::Sqrt(2))
				{
					Float_t xc2 = xvec[0][1];
					Float_t yc2 = yvec[0][1];
					Float_t xc1 = xvec[0][0];
					Float_t yc1 = yvec[0][0];
					Float_t phiAng0 = getPhiAngleDeg(xc2-xc1,yc2-yc1);	
					xc2 = xvec[1][1];
					yc2 = yvec[1][1];
					xc1 = xvec[1][0];
					yc1 = yvec[1][0];
					Float_t phiAng1mirrored = getPhiAngleDeg(xc1-xc2,yc2-yc1);
					Float_t phiAng1mirrored_shifted = phiAng1mirrored + relativePhiAngle;
					if (phiAng1mirrored_shifted > 360) phiAng1mirrored_shifted -= 360;
					Float_t dPh = phiAng0-phiAng1mirrored_shifted;
					if (dPh > 180) dPh = 360 - dPh;
					if (dPh < -180) dPh = -360 - dPh;
					dPhiAngle_Ewin->Fill(dPh);
					dPhiAngleNorm_Ewin->Fill(dPh);
					dPhiAngle1_Ewin->Fill(TMath::Abs(dPh));
					dPhiAngle1norm_Ewin->Fill(TMath::Abs(dPh));
				}
			}
		}
	
		// Fill data for AM trees
		for(Int_t im = 0; im < nAMs; im++)
		{
			for (Int_t ii = 1; ii <= nPixXY; ii++)
			{
				for (Int_t jj = 1; jj <= nPixXY; jj++)
				{
					if (plane[amID[im]]->GetBinContent(ii,jj) < 5) continue;
					event_e[amID[im]] = eventCounter-1;
					nTrigPixels_e[amID[im]] = nTrigPixels[amID[im]];
					AM_e[amID[im]] = amID[im];
					GM_e[amID[im]] = 1;
					nAMsInEvent_e[amID[im]] = 2;
					pixel_e[amID[im]] = pixelPattern[amID[im]]->GetBinContent(ii,jj);
					E[amID[im]] = plane[amID[im]]->GetBinContent(ii,jj);
					E_neg[amID[im]] = 0;
					time_e[amID[im]] = 1000;
					if(E[amID[im]] >= triggerThresh) triggerFlag_e[amID[im]] = 1;
					else triggerFlag_e[amID[im]] = 0;
					nloop_e[amID[im]] = 0;
					timeDetect_e[amID[im]] = 1000;
					timeDetectPos_e[amID[im]] = 1000;
					Temp_e[amID[im]] = 0;
					pos_e[amID[im]] = 0;
					E0_e[amID[im]] = energy[amID[im]][0];
					E0_x[amID[im]] = xvec[amID[im]][0];
					E0_y[amID[im]] = yvec[amID[im]][0];
					E0_z[amID[im]] = zvec[amID[im]][0];
					E0_pixel[amID[im]] = pixel[amID[im]][0];
					E1_e[amID[im]] = energy[amID[im]][1];
					E1_x[amID[im]] = xvec[amID[im]][1];
					E1_y[amID[im]] = yvec[amID[im]][1];
					E1_z[amID[im]] = zvec[amID[im]][1];
					E1_pixel[amID[im]] = pixel[amID[im]][1];
					if (amID[im] == 0)
					{
						theta_sim[amID[im]] = theta1;
						phi_sim[amID[im]] = phi1;
					}
					if (amID[im] == 1)	
					{
						theta_sim[amID[im]] = theta2;
						phi_sim[amID[im]] = phi2;
					}
					tree2[amID[im]]->Fill();
					//cout << "AM = " << AM_e[amID[im]] << ", pixel_e = " << pixel_e[amID[im]] << ", pixel = " << pixel[amID[im]][0] << ", E = " << E[amID[im]] << endl;
				}
			}
			
			//plane[amID[im]]->Draw("colz");
			//plane[amID[im]]->Draw("same,text");
			//c0->SaveAs(Form("event%ld_image_AM%d.png",event_e[amID[im]],amID[im]));
			//plane2[amID[im]]->Draw("colz");
			//plane2[amID[im]]->Draw("same,text");
			//c0->SaveAs(Form("event%ld_image2_AM%d.png",event_e[amID[im]],amID[im]));
		}
		//cout << "------------------------" << endl;
		// Fill data for event tree
		event_eev = eventCounter-1;
		nAMsInEvent_eev = 2;
		for(Int_t im = 0; im < nAMs; im++)
		{
			nTrigPixels_eev[amID[im]] = nTrigPixels_e[amID[im]];
			AM_flag_eev[amID[im]] = 1;
			cath_flag_eev[amID[im]] = 0;
		}
		savedEvents++;
		tree2_ev->Fill();
	}
	// End of loop over input file
	cout << "savedEvents = " << savedEvents << endl;
	
 	if (totalTimeElapced < 60) cout << Form("Total running time is %.2f sec",totalTimeElapced) << endl;
	if (totalTimeElapced >= 60 && totalTimeElapced < 3600) cout << Form("Total running time is %dm:%.0f sec",Int_t(totalTimeElapced/60),totalTimeElapced - Int_t(totalTimeElapced/60)*60) << endl;
	if (totalTimeElapced >= 3600) cout << Form("Total running time is %dh:%dm:%.0fs",Int_t(totalTimeElapced/3600),Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60),
		totalTimeElapced - Int_t(totalTimeElapced/3600)*3600 - Int_t((totalTimeElapced - Int_t(totalTimeElapced/3600)*3600)/60)*60) << endl;
  
	// Write trees and close files
	f2->cd();
	tree2_ev->Write();
	for (Int_t im = 0; im < nAMs; im++) tree2[im]->Write();
	f2->Close();
	cout << "Output written to file: " << rfileOut << endl;
	
	UInt_t screen_width = screen_width_def;
	UInt_t screen_height = screen_height_def;
	Float_t size_factor = 0.95;
	Float_t size_ratio = 1.333;
	TCanvas *c0 = new TCanvas("c0","",10,10,Int_t(screen_height*size_factor*size_ratio),Int_t(screen_height*size_factor));

	/*
	const Float_t screen_width_def = 1600;
	const Float_t screen_height_def = 900;
	UInt_t screen_width = screen_width_def;
	UInt_t screen_height = screen_height_def;
	Float_t size_factor = 0.95;
	Float_t size_ratio = 1.333;
	TCanvas *c0 = new TCanvas("c0","",10,10,Int_t(screen_height*size_factor*size_ratio),Int_t(screen_height*size_factor));
	*/
	usedPixelsH[0]->Draw();
	c0->SaveAs("usedPixelsH_AM0.png");
	usedPixelsH[1]->Draw();
	c0->SaveAs("usedPixelsH_AM1.png");
	pixelPattern[0]->Draw("text");
	c0->SaveAs("pixelPattern_AM0.png");
	pixelPattern[1]->Draw("text");
	c0->SaveAs("pixelPattern_AM1.png");
	plane[0]->Draw("colz");
	c0->SaveAs("image_AM0.png");
	plane[1]->Draw("colz");
	c0->SaveAs("image_AM1.png");
	
	dPhiAngle1norm_Ewin->SetMinimum(0);
	Float_t nf1 = (dPhiAngle1norm_Ewin->GetBinContent(1) + dPhiAngle1norm_Ewin->GetBinContent(dPhiAngle1norm_Ewin->GetNbinsX()))/2;
	for (Int_t ib = 1; ib <= dPhiAngle1norm_Ewin->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiAngle1norm_Ewin->GetBinContent(ib)/nf1;
		dPhiAngle1norm_Ewin->SetBinContent(ib,bc);
	}
	dPhiAngle1norm_Ewin->Draw("hist");
	c0->SaveAs(Form("%s.png",TString(dPhiAngle1norm_Ewin->GetName()).Data()));
	
	dPhiAngleNorm_Ewin->SetMinimum(0);
	Int_t b1 = dPhiAngleNorm_Ewin->GetXaxis()->FindBin(-0.01);
	Int_t b2 = dPhiAngleNorm_Ewin->GetXaxis()->FindBin(0.01);
	Float_t nf = (dPhiAngleNorm_Ewin->GetBinContent(b1) + dPhiAngleNorm_Ewin->GetBinContent(b2))/2;
	for (Int_t ib = 1; ib <= dPhiAngleNorm_Ewin->GetNbinsX(); ib++)
	{
		Float_t bc = dPhiAngleNorm_Ewin->GetBinContent(ib)/nf;
		dPhiAngleNorm_Ewin->SetBinContent(ib,bc);
	}
	dPhiAngleNorm_Ewin->Draw("hist");
	c0->SaveAs(Form("%s.png",TString(dPhiAngleNorm_Ewin->GetName()).Data()));
	
	dPhiAngle_Ewin->SetMinimum(0);
	dPhiAngle_Ewin->Draw("hist");
	c0->SaveAs(Form("%s.png",TString(dPhiAngle_Ewin->GetName()).Data()));
	
	dPhiAngle1_Ewin->SetMinimum(0);
	dPhiAngle1_Ewin->Draw("hist");
	c0->SaveAs(Form("%s.png",TString(dPhiAngle1_Ewin->GetName()).Data()));
	
	delete c0;
	
	TFile *hfile2 = new TFile("histos_phi.root","recreate");
	hfile2->cd();
	dPhiAngle1norm_Ewin->Write();
	dPhiAngleNorm_Ewin->Write();
	dPhiAngle_Ewin->Write();
	dPhiAngle1_Ewin->Write();
	hfile2->Close();
}

Int_t findAM(Float_t z)
{
   Int_t am = -1;
   if(z >= Zoffset) am = 0;
   else if(z <= -Zoffset) am = 1;
   
   return am;
}

Int_t findPixelID(Int_t am, Float_t x, Float_t y)
{
   Int_t pixelID=-1;
   Int_t row=-1, col=-1;
   
   // AM 1 if rotated 180 degrees so invert x
   if(am==1) x=-x;
   
   // added 100um to extremes
   Float_t pixEdges[12]={-4.5,-3.6,-2.8,-2.0,-1.2,-0.4,0.4,1.2,2.0,2.8,3.6,4.5};
   for(Int_t i=0; i<12; i++) {
      if(x>pixEdges[i] && x<=pixEdges[i+1]) col=i+1;
      if(y>pixEdges[i] && y<=pixEdges[i+1]) row=11-i;
   }
//   if(row!=-1 && col!=-1) pixelID=col+(row-1)*11;
   if(row!=-1 && col!=-1) pixelID=col+(row-1)*44;
   
//   if(pixelID!=-1){
   if(am==0) return pixelID;
   else if(am==1 && pixelID!=-1) return pixelID+22;
   else return pixelID;
}

Double_t sigmaTotal2_t_scaled(Double_t *x, Double_t *par)
{
	Double_t sigma2 = 2*D*x[0]*1e-6;
	Double_t c = -sigma2;
	Double_t d = -sigma2*TMath::Qe()*TMath::Qe()*N/20/k/T/EpsilonR/Epsilon0/TMath::Pi()/TMath::Sqrt(5.);
	
	return TMath::Hypot(TMath::Sqrt(2*D*x[0]*1e-6)*1e6*G_FWHM/2,TMath::Power(3.*mu_e*N*TMath::Qe()*x[0]*1e-6/4/TMath::Pi()/EpsilonR/Epsilon0,1./3)*1e6);
}

Bool_t getPixel2Dcoordinates(const Int_t pixx, const Int_t am, Int_t &x1, Int_t &y1)
{
	for (Int_t ix = 1; ix <= pixelPattern[am]->GetNbinsX(); ix++)
	{
		for (Int_t iy = 1; iy <= pixelPattern[am]->GetNbinsY(); iy++)
		{
			if (pixelPattern[am]->GetBinContent(ix,iy) == pixx)
			{
				x1 = ix;
				y1 = iy;
				return kTRUE;
			}
		}
	}
	return kFALSE;
}

Float_t getPhiAngleDeg(const Float_t x, const Float_t y)
{
	Float_t ang = -9999;
	if (x >= 0 && y >= 0) ang = TMath::ATan(y/x)*TMath::RadToDeg();
	if (x < 0 && y >= 0) ang = 180 - TMath::ATan(y/fabs(x))*TMath::RadToDeg();
	if (x < 0 && y < 0) ang = 180 + TMath::ATan(fabs(y/x))*TMath::RadToDeg();
	if (x >= 0 && y < 0) ang = 360 - TMath::ATan(fabs(y)/x)*TMath::RadToDeg();
	return ang;
}
