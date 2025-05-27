#ifndef AnalysisHelper_local
#define AnalysisHelper_local

#include "TSystem.h"
#include "TROOT.h"
#include "Rtypes.h"

#include "TMathBase.h"
#include "TMath.h"
#include "TMatrixDSym.h"

#include "TFile.h"
#include "TDirectoryFile.h"
#include "THashList.h"

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraphAsymmErrors.h"

#include "TString.h"

#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TLine.h"

#include "TColor.h"
#include "TStyle.h"

class AnalysisHelper
{

public:
	enum ETriggers
	{
		kMB = BIT(0),				  // Minimum bias trigger in PbPb 2010-11
		kINT1 = BIT(0),				  // V0A | V0C | SPD minimum bias trigger
		kINT7 = BIT(1),				  // V0AND minimum bias trigger
		kMUON = BIT(2),				  // Single muon trigger in pp2010-11, INT1 suite
		kHighMult = BIT(3),			  // High-multiplicity SPD trigger
		kHighMultSPD = BIT(3),		  // High-multiplicity SPD trigger
		kEMC1 = BIT(4),				  // EMCAL trigger in pp2011, INT1 suite
		kCINT5 = BIT(5),			  // V0OR minimum bias trigger
		kINT5 = BIT(5),				  // V0OR minimum bias trigger
		kCMUS5 = BIT(6),			  // Single muon trigger, INT5 suite
		kMUSPB = BIT(6),			  // Single muon trigger in PbPb 2011
		kINT7inMUON = BIT(6),		  // INT7 in MUON or MUFAST cluster
		kMuonSingleHighPt7 = BIT(7),  // Single muon high-pt, INT7 suite
		kMUSH7 = BIT(7),			  // Single muon high-pt, INT7 suite
		kMUSHPB = BIT(7),			  // Single muon high-pt in PbPb 2011
		kMuonLikeLowPt7 = BIT(8),	  // Like-sign dimuon low-pt, INT7 suite
		kMUL7 = BIT(8),				  // Like-sign dimuon low-pt, INT7 suite
		kMuonLikePB = BIT(8),		  // Like-sign dimuon low-pt in PbPb 2011
		kMuonUnlikeLowPt7 = BIT(9),	  // Unlike-sign dimuon low-pt, INT7 suite
		kMUU7 = BIT(9),				  // Unlike-sign dimuon low-pt, INT7 suite
		kMuonUnlikePB = BIT(9),		  // Unlike-sign dimuon low-pt in PbPb 2011
		kEMC7 = BIT(10),			  // EMCAL/DCAL L0 trigger, INT7 suite
		kEMC8 = BIT(10),			  // EMCAL/DCAL L0 trigger, INT8 suite
		kMUS7 = BIT(11),			  // Single muon low-pt, INT7 suite
		kMuonSingleLowPt7 = BIT(11),  // Single muon low-pt, INT7 suite
		kPHI1 = BIT(12),			  // PHOS L0 trigger in pp2011, INT1 suite
		kPHI7 = BIT(13),			  // PHOS trigger, INT7 suite
		kPHI8 = BIT(13),			  // PHOS trigger, INT8 suite
		kPHOSPb = BIT(13),			  // PHOS trigger in PbPb 2011
		kEMCEJE = BIT(14),			  // EMCAL/DCAL L1 jet trigger
		kEMCEGA = BIT(15),			  // EMCAL/DCAL L1 gamma trigger
		kHighMultV0 = BIT(16),		  // High-multiplicity V0 trigger
		kCentral = BIT(16),			  // Central trigger in PbPb 2011
		kSemiCentral = BIT(17),		  // Semicentral trigger in PbPb 2011
		kDG = BIT(18),				  // Double gap diffractive
		kDG5 = BIT(18),				  // Double gap diffractive
		kZED = BIT(19),				  // ZDC electromagnetic dissociation
		kSPI7 = BIT(20),			  // Power interaction trigger
		kSPI = BIT(20),				  // Power interaction trigger
		kINT8 = BIT(21),			  // 0TVX trigger
		kMuonSingleLowPt8 = BIT(22),  // Single muon low-pt, INT8 suite
		kMuonSingleHighPt8 = BIT(23), // Single muon high-pt, INT8 suite
		kMuonLikeLowPt8 = BIT(24),	  // Like-sign dimuon low-pt, INT8 suite
		kMuonUnlikeLowPt8 = BIT(25),  // Unlike-sign dimuon low-pt, INT8 suite
		kMuonUnlikeLowPt0 = BIT(26),  // Unlike-sign dimuon low-pt, no additional L0 requirement
		kUserDefined = BIT(27),		  // Set when custom trigger classes are set in AliPhysicsSelection
		kTRD = BIT(28),				  // TRD trigger
		kMuonCalo = BIT(29),		  // Muon-calo triggers
		kCaloOnly = BIT(29),		  // MB, EMCAL and PHOS triggers in CALO or CALOFAST cluster
		// Bits 30 and above are reserved for FLAGS
		kFastOnly = BIT(30),						  // The fast cluster fired. This bit is set in to addition another trigger bit, e.g. kMB
		kAny = 0xffffffff,							  // to accept any defined trigger
		kAnyINT = kMB | kINT7 | kINT5 | kINT8 | kSPI7 // to accept any interaction (aka minimum bias) trigger
	};

	AnalysisHelper(const char *filename, TString CollisionSystem, TString Centrality)
	{

		// TString test = "20240215/cen0_10/CPV25/Disp25/200MeV/CoreE/30ns";

		fInputFile = new TFile(Form("%s.root", filename));
		fCollisionSystem = CollisionSystem;
		fCentRange = Centrality;

		if (fCollisionSystem.Contains("PbPb"))
		{
			if (!fCentRange.CompareTo("0-10%"))
				fCentColor = kRed + 1;
			else if (!fCentRange.CompareTo("0-20%"))
				fCentColor = kRed + 1;
			else if (!fCentRange.CompareTo("0-40%"))
				fCentColor = kMagenta + 2;
			else if (!fCentRange.CompareTo("0-5%"))
				fCentColor = kRed + 1;
			else if (!fCentRange.CompareTo("5-10%"))
				fCentColor = 807;
			else if (!fCentRange.CompareTo("10-20%"))
				fCentColor = 800;
			else if (!fCentRange.CompareTo("20-30%"))
				fCentColor = kSpring + 5;
			else if (!fCentRange.CompareTo("10-30%"))
				fCentColor = kYellow - 3;
			else if (!fCentRange.CompareTo("20-40%") || !fCentRange.CompareTo("30-40%"))
				fCentColor = kGreen + 2;
			else if (!fCentRange.CompareTo("20-50%"))
				fCentColor = kCyan + 2;
			else if (!fCentRange.CompareTo("30-50%") || !fCentRange.CompareTo("40-50%"))
				fCentColor = kTeal + 4;
			else if (!fCentRange.CompareTo("40-60%") || !fCentRange.CompareTo("50-60%"))
				fCentColor = kCyan + 2;
			else if (!fCentRange.CompareTo("60-70%"))
				fCentColor = kAzure + 5;
			else if (!fCentRange.CompareTo("70-80%"))
				fCentColor = kAzure - 1;
			else if (!fCentRange.CompareTo("60-80%") || !fCentRange.CompareTo("60-92%"))
				fCentColor = kBlue + 1;
			else if (!fCentRange.CompareTo("40-80%") || !fCentRange.CompareTo("0-94%"))
				fCentColor = kCyan + 2;
			else if (!fCentRange.CompareTo("50-90%"))
				fCentColor = kViolet + 2;
			else if (!fCentRange.CompareTo("80-90%"))
				fCentColor = kViolet + 4;
			else if (!fCentRange.CompareTo("30-100%"))
				fCentColor = kViolet + 2;
		}

		SetCentrality(fCentRange);

		fIsRawData = kFALSE;
	};
	AnalysisHelper(){};
	// AnalysisHelper(TFile *file, Int_t CollisionSystem, Int_t CentBin);

	void SetInputFile(const char *filename)
	{
		fInputFile = TFile::Open(Form("%s.root", filename));
	};

	void SetCentrality(TString cent){
		fCentLowerLimit = TString(cent(0, cent.First('-'))).Atoi();
		fCentUpperLimit = TString(cent(cent.First('-') + 1, cent.First('%'))).Atoi();
	};

	void SetRawDataFlag(Bool_t flag = kFALSE) { fIsRawData = flag; };

	void SetRunConditions(const char *name = "Pi0EtaToGammaGamma",
						  UInt_t trigger = AnalysisHelper::kINT7,
						  const TString CollisionSystem = "PbPb",
						  const Bool_t isMC = kFALSE,
						  const Int_t L1input = -1, // L1H,L1M,L1L
						  const Int_t L0input = -1, // L0
						  const Float_t CenMin = 0.,
						  const Float_t CenMax = 90.,
						  const Int_t NMixed = 10,
						  const Bool_t FlowTask = kFALSE,
						  const Int_t harmonics = -1,
						  const Int_t FlowMethod = -1,
						  const Int_t QnDetector = -1,
						  const Bool_t useCoreE = kFALSE,
						  const Bool_t useCoreDisp = kFALSE,
						  const Double_t NsigmaCPV = 2.5,
						  const Double_t NsigmaDisp = 2.5,
						  const Bool_t usePHOSTender = kTRUE,
						  const Bool_t TOFcorrection = kTRUE,
						  const Bool_t Trgcorrection = kFALSE,
						  const Bool_t NonLinStudy = kFALSE,
						  const Double_t bs = 100.,	  // bunch space in ns.
						  const Double_t distBC = -1, // minimum distance to bad channel.
						  const Double_t Emin = 0.2,  // minimum energy for photon selection in GeV
						  const Bool_t isJJMC = kFALSE,
						  const TString MCtype = "MBMC",
						  const Bool_t ForceActiveTRU = kFALSE,
						  const Bool_t ApplyTOFTrigger = kFALSE,
						  const char *sub_name = "0_10" // options for subwagon names
	)

	{

		fIsRawData = kTRUE;

		TString TriggerName = "";
		if (trigger == (UInt_t)AnalysisHelper::kAny)
			TriggerName = "kAny";
		else if (trigger == (UInt_t)AnalysisHelper::kINT7)
			TriggerName = "kINT7";
		else if (trigger & AnalysisHelper::kPHI7)
			TriggerName = "kPHI7";

		if (trigger & AnalysisHelper::kPHI7)
		{
			if (L1input > 0)
			{
				if (L1input == 7)
					TriggerName = TriggerName + "_" + "L1H";
				else if (L1input == 6)
					TriggerName = TriggerName + "_" + "L1M";
				else if (L1input == 5)
					TriggerName = TriggerName + "_" + "L1L";
			}
			else if (L0input > 0)
				TriggerName = TriggerName + "_" + "L0";
			else
			{
				printf("PHOS trigger analysis requires at least 1 trigger input (L0 or L1[H,M,L]).");
				return;
			}
		}

		Int_t systemID = -1;
		if (CollisionSystem == "pp")
			systemID = 0;
		else if (CollisionSystem == "PbPb")
			systemID = 1;
		else if (CollisionSystem == "pPb" || CollisionSystem == "Pbp")
			systemID = 2;

		TString PIDname = "";
		if (NsigmaCPV > 0)
			PIDname += Form("_CPV%d", (Int_t)(NsigmaCPV * 10));
		if (NsigmaDisp > 0)
		{
			if (useCoreDisp)
				PIDname += Form("_CoreDisp%d", (Int_t)(NsigmaDisp * 10));
			else
				PIDname += Form("_FullDisp%d", (Int_t)(NsigmaDisp * 10));
		}
		if (useCoreE)
			PIDname += "_CoreE";
		else
			PIDname += "_FullE";

		TString taskname = "";

		taskname = Form("%s_%s_%s_Cen%d_%d%s_BS%dns_DBC%dcell_Emin%dMeV",
						name, CollisionSystem.Data(), TriggerName.Data(),
						(Int_t)CenMin, (Int_t)CenMax, PIDname.Data(),
						(Int_t)bs, (Int_t)(distBC), (Int_t)(Emin * 1e+3));

		if ((trigger == (UInt_t)AnalysisHelper::kPHI7) && ApplyTOFTrigger)
			taskname += "_TOFTrigger";
		if (ForceActiveTRU)
			taskname += "_ForceActiveTRU";

		taskname += Form("_%s", sub_name);

		fHashListName = Form("hist_%s", taskname.Data());
		// printf("Could not find THashList %s\n", fHashListName.Data());
	};

	void SetCollisionSystem(TString collsys) { fCollisionSystem = collsys; };

	Int_t GetCentUpperLimit() { return fCentUpperLimit; };
	Int_t GetCentLowerLimit() { return fCentLowerLimit; };

	void SetTH1Style(TH1 *hist,
					 Style_t markerStyle,
					 Color_t markerColor,
					 Color_t lineColor)
	{

		hist->SetMarkerStyle(markerStyle);
		hist->SetMarkerColor(markerColor);
		hist->SetLineColor(lineColor);
	};

	void SetCentStyle(TH1 *hist)
	{
		hist->SetMarkerStyle(20);
		hist->SetMarkerSize(1.0);
		hist->SetMarkerColor(GetCentralityColor());
		hist->SetLineColor(GetCentralityColor());
	};

	// void SetTH1Binning(TH1 *hist)
	// {
	// 	hist->Rebin(fNbinsPt - 1, hist->GetName(), fPtAxis);
	// };

	TObject *Get(const char *name)
	{

		TObject *obj = nullptr;

		if (fCollisionSystem.Contains("PbPb"))
		{
			if (fIsRawData)
			{
				TDirectoryFile *dir = (TDirectoryFile *)fInputFile->Get("PWGGA_PHOSTasks_PHOSRun2");
				THashList *list = (THashList *)dir->Get(fHashListName);
				obj = list->FindObject(name);
			}
			else
			{
				obj = fInputFile->Get(name);
			}
		}

		return obj;
	};

	Color_t GetCentralityColor() { return fCentColor; };
	Style_t GetCentMarkerStyle() { return fCentMarker; };

private:
	TFile *fInputFile;
	Bool_t fIsRawData;
	TString fHashListName;

	TString fCollisionSystem;
	TString fCentRange;
	Int_t fCentBin;
	Int_t fCentUpperLimit;
	Int_t fCentLowerLimit;

	// Int_t				fNbinsPt;
	// Double_t            *fPtAxis;

	Color_t fCentColor;
	Style_t fCentMarker;
};

#endif


TCanvas *DefaultCanvas(){

		// TCanvas *canv = new TCanvas("DefaultCanvas", "DefaultCanvas", 800,600);
		TCanvas *canv = new TCanvas();
		
		// canv->SetGrid();
		canv->SetTicks();
		canv->SetTopMargin(0.02);
		canv->SetBottomMargin(0.1);
		canv->SetRightMargin(0.02);
		canv->SetLeftMargin(0.1);
		canv->Draw();
		canv->cd();

		return canv;

}

TCanvas *DrawRatio(TH1F *hist1, TH1F *hist2, TH1F *hRatio, TLegend *l){

	TCanvas *canv = new TCanvas();	

	TPad *pad1 = new TPad("pad1", "pad1", 0., 0.0125, 1., 1.);
	// pad1->SetGrid();
	pad1->SetTicks();
	pad1->SetBottomMargin(0.3);
	pad1->SetRightMargin(0.03);
	pad1->SetTopMargin(0.01);
	pad1->SetLeftMargin(0.1);
	pad1->Draw();
	pad1->cd();

	hist1->Draw("E0");
	hist2->Draw("same E0");
	l->Draw("same");

	TPad *pad2 = new TPad("pad2", "pad2", 0., 0.0045, 1., 0.3);
	// pad2->SetGrid();
	pad2->SetTicks();
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->SetRightMargin(0.03);
	pad2->SetLeftMargin(0.1);
	pad2->Draw();
	pad2->cd();

	hRatio->GetYaxis()->SetNdivisions(504);
	hRatio->GetYaxis()->SetLabelSize(0.04*3);
	hRatio->GetYaxis()->SetTitleSize(0.04*3);
	hRatio->GetYaxis()->SetTitleOffset(0);
	hRatio->GetXaxis()->SetTickLength(0.04*3);
	hRatio->GetXaxis()->SetLabelSize(0.04*3);
	hRatio->GetXaxis()->SetTitleSize(0.04*3);


	hRatio->Draw("E0");

	return canv;

}

TCanvas *DrawRatio(TH1F *hist1, TH1F *hist2, TString RatioTitle = "MC/Data"){

	TCanvas *canv = new TCanvas();	

	TPad *pad1 = new TPad("pad1", "pad1", 0., 0.0125, 1., 1.);
	// pad1->SetGrid();
	pad1->SetTicks();
	pad1->SetBottomMargin(0.3);
	pad1->SetRightMargin(0.03);
	pad1->SetTopMargin(0.01);
	pad1->SetLeftMargin(0.1);
	pad1->Draw();
	pad1->cd();

	hist1->Draw("E0");
	hist2->Draw("same E0");

	TPad *pad2 = new TPad("pad2", "pad2", 0., 0.0045, 1., 0.3);
	// pad2->SetGrid();
	pad2->SetTicks();
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->SetRightMargin(0.03);
	pad2->SetLeftMargin(0.1);
	pad2->Draw();
	pad2->cd();

	TH1F *hRatio = (TH1F*)hist1->Clone(Form("%s_Ratio", hist1->GetName()));
	hRatio->Divide(hRatio, hist2, 1., 1., "B");
	hRatio->SetTitle("");
	hRatio->GetXaxis()->SetTitle(Form("%s", hist1->GetXaxis()->GetTitle()));
	hRatio->GetYaxis()->SetTitle(RatioTitle);


	hRatio->GetYaxis()->SetNdivisions(504);
	hRatio->GetYaxis()->SetLabelSize(0.04*3);
	hRatio->GetYaxis()->SetTitleSize(0.04*3);
	hRatio->GetYaxis()->SetTitleOffset(0);
	hRatio->GetXaxis()->SetTickLength(0.04*3);
	hRatio->GetXaxis()->SetLabelSize(0.04*3);
	hRatio->GetXaxis()->SetTitleSize(0.04*3);


	hRatio->Draw("E0");

	return canv;

}

TCanvas *DrawRatio(TH1F *hist1, TF1 *fit1, TString RatioTitle = "Data/fit"){

	TCanvas *canv = new TCanvas();	

	TPad *pad1 = new TPad("pad1", "pad1", 0., 0.0125, 1., 1.);
	pad1->SetGrid();
	pad1->SetTicks();
	pad1->SetBottomMargin(0.3);
	pad1->SetTopMargin(0.02);
	pad1->SetRightMargin(0.03);
	pad1->SetLeftMargin(0.1);
	pad1->Draw();
	pad1->cd();

	hist1->Draw("E0");
	fit1->Draw("same");

	TPad *pad2 = new TPad("pad2", "pad2", 0., 0.0045, 1., 0.3);
	pad2->SetGrid();
	pad2->SetTicks();
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->SetRightMargin(0.03);
	pad2->SetLeftMargin(0.1);
	pad2->Draw();
	pad2->cd();

	TH1F *hRatio = (TH1F*)hist1->Clone(Form("%s_Ratio", hist1->GetName()));
	hRatio->Divide(fit1);
	hRatio->SetTitle("");
	hRatio->GetXaxis()->SetTitle(Form("%s", hist1->GetXaxis()->GetTitle()));
	hRatio->GetYaxis()->SetTitle(RatioTitle);


	hRatio->GetYaxis()->SetNdivisions(504);
	hRatio->GetYaxis()->SetLabelSize(0.04*3);
	hRatio->GetYaxis()->SetTitleSize(0.04*3);
	hRatio->GetYaxis()->SetTitleOffset(0);
	hRatio->GetXaxis()->SetTickLength(0.04*3);
	hRatio->GetXaxis()->SetLabelSize(0.04*3);
	hRatio->GetXaxis()->SetTitleSize(0.04*3);


	hRatio->Draw("E0");

	return canv;

}

void SetHistStyle(TH1* hist,
                  Style_t markerStyle,
                  Size_t markerSize,
                  Color_t markerColor,
                  Color_t lineColor){

	hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerSize(markerSize);
    hist->SetMarkerColor(markerColor);
    hist->SetLineColor(lineColor);
}

void SetHistStyle(TH1* hist,
                  Style_t markerStyle,
                  Color_t Color){

	hist->SetMarkerStyle(markerStyle);
    hist->SetMarkerColor(Color);
    hist->SetLineColor(Color);
}

void SetDrawTH1Style(TCanvas *canv,
					 TH1 *histo,
					 TLegend *l){

	Double_t textsizeLabels, textsizeFac;
	Double_t textSizeLabelsPixel  = 40*4/5;

    Float_t xTitleOffset    = 0.85;
    Float_t yTitleOffset    = 0.85;
    Int_t xNDivisions       = 510;
    Int_t yNDivisions       = 510;
    Font_t textFontLabel    = 42;
    Font_t textFontTitle    = 52;

	if (canv->XtoPixel(canv->GetX2()) < canv->YtoPixel(canv->GetY1())){
        textsizeLabels       = (Double_t)textSizeLabelsPixel/canv->XtoPixel(canv->GetX2());
        textsizeFac          = (Double_t)1./canv->XtoPixel(canv->GetX2());
    } else {
        textsizeLabels       = (Double_t)textSizeLabelsPixel/canv->YtoPixel(canv->GetY1());
        textsizeFac          = (Double_t)1./canv->YtoPixel(canv->GetY1());
    }

    histo->GetYaxis()->SetLabelFont(textFontLabel);
    histo->GetXaxis()->SetLabelFont(textFontLabel);
    histo->GetYaxis()->SetTitleFont(textFontLabel);
    histo->GetXaxis()->SetTitleFont(textFontLabel);

    histo->GetXaxis()->SetLabelSize(0.6*textsizeLabels);
    histo->GetXaxis()->SetTitleSize(0.75*textsizeLabels);
    histo->GetXaxis()->SetTitleOffset(xTitleOffset);
    histo->GetXaxis()->SetNdivisions(xNDivisions,kTRUE);

    histo->GetYaxis()->SetDecimals();
    histo->GetYaxis()->SetLabelSize(0.6*textsizeLabels);
    histo->GetYaxis()->SetTitleSize(0.75*textsizeLabels);
    histo->GetYaxis()->SetTitleOffset(yTitleOffset);
    histo->GetYaxis()->SetNdivisions(yNDivisions,kTRUE);

    l->SetTextFont(textFontLabel);
	l->SetTextSize(0.75*textsizeLabels);


// 									TH2* histo,
//                                     TString XTitle,
//                                     TString YTitle,
//                                     Size_t xLableSize,
//                                     Size_t xTitleSize,
//                                     Size_t yLableSize,
//                                     Size_t yTitleSize,
//                                     Float_t xTitleOffset    = 1,
//                                     Float_t yTitleOffset    = 1,
//                                     Int_t xNDivisions       = 510,
//                                     Int_t yNDivisions       = 510,
//                                     Font_t textFontLabel    = 42,
//                                     Font_t textFontTitle    = 62

// SetStyleHistoTH2ForGraphs(histo2DPi0M02Dummy, "#it{#sigma}_{long}^{2} (cm^{2})","counts (a.u.)",0.85*textsizeLabelsM02, textsizeLabelsM02,
//             0.85*textsizeLabelsM02, textsizeLabelsM02,0.88, 1.15); //

}

// struct TH1DrawArray{

// 	Int_t stacksize = 10;
// 	std::vector<TH1F*> hists1;
// 	std::vector<TH1F*> hists2;
// 	std::vector<TH1F*> ratios;
// 	TCanvas *canv;

// 	hists1.resize(stacksize);
// 	hists2.resize(stacksize);
// 	ratios.resize(stacksize);

// 	void Init(Int_t nhists = 10, std::vector<TH1F*> InputHists1, std::vector<TH1F*> InputHists2){

// 		stacksize = nhists;

// 		hists1.resize(stacksize);
// 		hists2.resize(stacksize);
// 		ratios.resize(stacksize);

// 		InputHists1.resize(stacksize);
// 		InputHists2.resize(stacksize);

// 		for (Int_t i = 0; i < stacksize; i++){
// 			if (InputHists1.at(i))
// 				hists1.at(i) = InputHists1.at(i);
// 			if (InputHists2.at(i))
// 				hists2.at(i) = InputHists2.at(i);
// 			if (hists1.at(i) && hists2.at(i)){
// 				TH1F *hRatio = (TH1F*)hists1.at(i)->Clone(Form("hRatio_%d", i));
// 				hRatio->Divide(hRatio, hists2.at(i), 1., 1., "B");
// 				ratios.at(i) = hRatio;
// 			}
// 		}

// 		DefaultCanvas();
// 	}

// 	void DefaultCanvas(){

// 		canv->SetGrid();
// 		canv->SetTicks();

// 		canv->cd();
// 	}


// 	void MakeRatio(Int_t ihist){
// 		if ()
// 		if (!hists1.at(ihist) && !hists2.at(ihist)){
// 			printf("No element %d in the array!!!\n", ihist);
// 			return;
// 		} else{
// 			TH1F *hRatio = (TH1F*)hists1.at(ihist)->Clone(Form("hRatio_%d", ihist));
// 			hRatio->Divide(hRatio, hists2.at(ihist), 1., 1., "B");
// 			ratios.at(ihist) = hRatio;
// 		}
// 	}

// }