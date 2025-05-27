/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Analysis task for neutral meson and inclusive/direct photon measurements (AODs only)
// + supproting QA studies
// Author: Vladislav Kuskov
// since march 2024

#include "TF1.h"
#include "TH3.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TProfile.h"

#include "TMath.h"

#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliGenEventHeader.h"

#include "AliPHOSEventCuts.h"
#include "AliPHOSClusterCuts.h"

#include "AliMultSelection.h"
#include "AliOADBContainer.h"

#include "THashList.h"

#include "AliAnalysisPHOSNeutralMesonsAndPhotons.h"

ClassImp(AliAnalysisPHOSNeutralMesonsAndPhotons);

//_____________________________________________________________________________
AliAnalysisPHOSNeutralMesonsAndPhotons::AliAnalysisPHOSNeutralMesonsAndPhotons(const char* name) : 
    AliAnalysisTaskSE(name),
    fOutputContainer(nullptr),
    fEvent(nullptr),
    fMCEvent(nullptr),
    fRunNumber(0),
    fMCStack(nullptr),
    fCaloPhotonsPHOS(nullptr),
    fCaloPhotonsMix(nullptr),
    fCaloPhotonsPHOSList(nullptr),
    fUtils(nullptr),
    fPHOSGeo(nullptr),
    fAODCells(nullptr),
    fPHOSClusterCuts(nullptr),
    fCaloTriggerMimicHelper(nullptr),
    fPHOSTrigType(AliCaloTriggerMimicHelper::kPHOSAny),
    fIsMC(kFALSE),
    fIsPHOSTriggerAnalysis(kFALSE),
    fIsMBTriggerAnalysis(kTRUE),
    fUseCoreEnergy(kTRUE),
    fDoNonlinCorr(kTRUE),
    fUsePHOSClusterCuts(kFALSE),
    fUserNonlinFunc(nullptr),
    fClearGammaWeight(kFALSE),
    fCentArrayPi0(nullptr),
    fCentArrayEta(nullptr),
    fCentArrayGamma(nullptr),
    fCentArrayK0S(nullptr),
    fDoClustQA(kFALSE),
    fDoCellQA(kFALSE),
    fNMixEvents(10),
    fCentralityEstimator("V0M"),
    fMultSelection(nullptr),
    fCentrality(-1),
    fCentralityMin(0.),
    fCentralityMax(10.),
    fNCenBin(5),
    fEminCut(0.2),
    fDispSigma(2.5),
    fCPVSigma(2.5),
    fTOFCut(30.),
    fMinBiasLowEClustCut(kTRUE),
    fHistInfo(nullptr),
    fHistSelectEvent(nullptr),
    fHistVertexZ(nullptr),
    fHistCentMain(nullptr),
    fHistCentV0MvsCL0(nullptr),
    fHistCentV0MvsCL1(nullptr),
    fHistCentV0MvsV0C(nullptr),
    fHistCentCL0vsCL1(nullptr),
    fHistCentZNAvsZNC(nullptr),
    fHistTOFClust(nullptr),
    fHistCaloPhotonTOFvsE(nullptr),
    fHistCaloPhotonTOFCut(nullptr),
    fHistClustTOFvsDDL(nullptr),
    fHistClustTOFvsDDLEnCut(nullptr),
    fHistClustMultVsCentrality(nullptr),
    fHistM02vsPt(nullptr),
    fHistM20vsPt(nullptr),
    fHistClustFullE(nullptr),
    fHistClustCoreE(nullptr),
    fHistNonlinTest(nullptr),
    fHistMggTOFCutEffBase(nullptr),
    fHistMixMggTOFCutEffBase(nullptr),
    fHistMggTOFCutEffProbe(nullptr),
    fHistMixMggTOFCutEffProbe(nullptr),
    fHistTruePi0MggVsRecPt(nullptr),
    fHistTruePi0MggVsRecPtFromK0s(nullptr),
    fHistTruePi0MggVsRecPtFromK0L(nullptr),
    fHistTruePi0MggVsRecPtFromLambda(nullptr),
    fHistTruePi0MggVsRecPtFromPandN(nullptr),
    fHistTruePi0MggVsRecPtFromPions(nullptr),
    fHistTruePi0MggVsRecPtFromKaons(nullptr),
    fHistTruePi0MggVsTruePt(nullptr),
    fHistTruePi0MotherIDvsRecPt(nullptr),
    fHistTruePi0MotherIDvsTruePt(nullptr),
    fHistTruePi0MotherIDvsRecPtW(nullptr),
    fHistTrueEtaMggVsRecPt(nullptr),
    fHistTrueEtaMggVsTruePt(nullptr),
    fHistMCPartIDvsPt(nullptr),
    fHistMCPartIDvsPtW(nullptr),
    fHistPrimPi0InAccPt(nullptr),
    fHistPrimPi0InAccPtW(nullptr),
    fHistPrimPi0BothPhInAccPt(nullptr),
    fHistPrimPi0BothPhInAccPtW(nullptr),
    fHistPrimPi0OnePhInAccPt(nullptr),
    fHistPrimPi0PtvsPhInAccPt(nullptr),
    fHistPrimEtaInAccPt(nullptr),
    fHistPrimEtaBothPhInAccPt(nullptr),
    fHistPrimEtaOnePhInAccPt(nullptr),
    fHistPrimEtaPtvsPhInAccPt(nullptr),
    fHistPrimEtaPt(nullptr),
    fHistPrimEtaPtW(nullptr),
    fHistPrimGammaPt(nullptr),
    fHistPrimGammaPtW(nullptr),
    fHistPrimGammaInAccPt(nullptr),
    fHistPrimGammaInAccPtW(nullptr)
{

  for (Int_t i = 0; i < kVtxBins; i++)
    for (Int_t j = 0; j < kCentBins; j++)
      for (Int_t k = 0; k < kPRBins; k++)
        fPHOSEvents[i][j][k] = nullptr;

  for (Int_t imod = 0; imod < kMods; imod++) {

    fPHOSModEnCorr[imod] = 1.;

    fHistTOFClustMod[imod] = nullptr;
    fHistClustMultVsCentralityMod[imod] = nullptr;
    fHistCaloPhotonTOFvsEMod[imod] = nullptr;
    fHistClustFullEMod[imod] = nullptr;
    fHistClustCoreEMod[imod] = nullptr;
    fHistCaloPhotonPtMod[imod] = nullptr;
    fHistCaloPhotonPtModW[imod] = nullptr;

    fHistMggMod[imod] = nullptr;
    fHistMixMggMod[imod] = nullptr;
    fHistMggPhIDCutMod[imod] = nullptr;
    fHistMixMggPhIDCutMod[imod] = nullptr;

    fHistClustEvsXZMod[imod] = nullptr;
  }

  for (Int_t icut = 0; icut < kPIDCuts; icut++) {
    fHistCaloPhotonPt[icut] = nullptr;
    fHistCaloPhoton2Pt[icut] = nullptr;
    fHistCaloPhotonPtW[icut] = nullptr;
    fHistMCCaloPartIDvsPt[icut] = nullptr;
    fHistMCCaloPartIDvsPtW[icut] = nullptr;

    fHistMgg[icut] = nullptr;
    fHistMgg2[icut] = nullptr;
    fHistMggW[icut] = nullptr;
    fHistMixMgg[icut] = nullptr;
    fHistMixMgg2[icut] = nullptr;
    fHistMixMggW[icut] = nullptr;

    fHistMggCutEff[icut] = nullptr;
    fHistMixMggCutEff[icut] = nullptr;
  }

  for (Int_t imap = 0; imap < 6; imap++)
    fPHOSBadMap[imap] = nullptr;

  for(Int_t i = 0; i < 20; i++){
    fPi0PtWeight[i]   = nullptr;
    fEtaPtWeight[i]   = nullptr;
    fK0SPtWeight[i]   = nullptr;
    fGammaPtWeight[i] = nullptr;

    fPi0PtWeight[i] = new TF1(Form("fPi0PtWeight_%d",i) ,"1.", 0, 200);
    fEtaPtWeight[i] = new TF1(Form("fEtaPtWeight_%d",i) ,"1.", 0, 200);
    fK0SPtWeight[i] = new TF1(Form("fK0SPtWeight_%d",i) ,"1.", 0, 200);
    fGammaPtWeight[i] = new TF1(Form("fGammaPtWeight_%d",i) ,"1.", 0, 200);
  }

  fUserNonlinFunc = new TF1("fUserNonlinFunc", "1.", 0, 200);

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, THashList::Class());
}
//_____________________________________________________________________________
AliAnalysisPHOSNeutralMesonsAndPhotons::AliAnalysisPHOSNeutralMesonsAndPhotons(const AliAnalysisPHOSNeutralMesonsAndPhotons& rh) : 
    AliAnalysisTaskSE(rh.GetName()),
    fOutputContainer(nullptr),
    fEvent(nullptr),
    fMCEvent(nullptr),
    fRunNumber(0),
    fMCStack(nullptr),
    fCaloPhotonsPHOS(nullptr),
    fCaloPhotonsMix(nullptr),
    fCaloPhotonsPHOSList(nullptr),
    fUtils(nullptr),
    fPHOSGeo(nullptr),
    fAODCells(nullptr),
    fPHOSClusterCuts(nullptr),
    fCaloTriggerMimicHelper(nullptr),
    fPHOSTrigType(AliCaloTriggerMimicHelper::kPHOSAny),
    fIsMC(kFALSE),
    fIsPHOSTriggerAnalysis(kFALSE),
    fIsMBTriggerAnalysis(kTRUE),
    fUseCoreEnergy(kTRUE),
    fDoNonlinCorr(kFALSE),
    fDoTOFEffCorr(kFALSE),
    fUsePHOSClusterCuts(kFALSE),
    fUserNonlinFunc(nullptr),
    fClearGammaWeight(kFALSE),
    fCentArrayPi0(nullptr),
    fCentArrayEta(nullptr),
    fCentArrayGamma(nullptr),
    fCentArrayK0S(nullptr),
    fDoClustQA(kFALSE),
    fDoCellQA(kFALSE),
    fNMixEvents(10),
    fCentralityEstimator("V0M"),
    fMultSelection(nullptr),
    fCentrality(-1),
    fCentralityMin(0.),
    fCentralityMax(10.),
    fNCenBin(5),
    fEminCut(0.2),
    fDispSigma(2.5),
    fCPVSigma(2.5),
    fTOFCut(30.),
    fMinBiasLowEClustCut(kTRUE),
    fHistInfo(nullptr),
    fHistSelectEvent(nullptr),
    fHistVertexZ(nullptr),
    fHistCentMain(nullptr),
    fHistCentV0MvsCL0(nullptr),
    fHistCentV0MvsCL1(nullptr),
    fHistCentV0MvsV0C(nullptr),
    fHistCentCL0vsCL1(nullptr),
    fHistCentZNAvsZNC(nullptr),
    fHistTOFClust(nullptr),
    fHistCaloPhotonTOFvsE(nullptr),
    fHistCaloPhotonTOFCut(nullptr),
    fHistClustTOFvsDDL(nullptr),
    fHistClustTOFvsDDLEnCut(nullptr),
    fHistClustMultVsCentrality(nullptr),
    fHistM02vsPt(nullptr),
    fHistM20vsPt(nullptr),
    fHistClustFullE(nullptr),
    fHistClustCoreE(nullptr),
    fHistNonlinTest(nullptr),
    fHistMggTOFCutEffBase(nullptr),
    fHistMixMggTOFCutEffBase(nullptr),
    fHistMggTOFCutEffProbe(nullptr),
    fHistMixMggTOFCutEffProbe(nullptr),
    fHistTruePi0MggVsRecPt(nullptr),
    fHistTruePi0MggVsRecPtFromK0s(nullptr),
    fHistTruePi0MggVsRecPtFromK0L(nullptr),
    fHistTruePi0MggVsRecPtFromLambda(nullptr),
    fHistTruePi0MggVsRecPtFromPandN(nullptr),
    fHistTruePi0MggVsRecPtFromPions(nullptr),
    fHistTruePi0MggVsRecPtFromKaons(nullptr),
    fHistTruePi0MggVsTruePt(nullptr),
    fHistTruePi0MotherIDvsRecPt(nullptr),
    fHistTruePi0MotherIDvsTruePt(nullptr),
    fHistTruePi0MotherIDvsRecPtW(nullptr),
    fHistTrueEtaMggVsRecPt(nullptr),
    fHistTrueEtaMggVsTruePt(nullptr),
    fHistMCPartIDvsPt(nullptr),
    fHistMCPartIDvsPtW(nullptr),
    fHistPrimPi0InAccPt(nullptr),
    fHistPrimPi0InAccPtW(nullptr),
    fHistPrimPi0BothPhInAccPt(nullptr),
    fHistPrimPi0BothPhInAccPtW(nullptr),
    fHistPrimPi0OnePhInAccPt(nullptr),
    fHistPrimPi0PtvsPhInAccPt(nullptr),
    fHistPrimEtaInAccPt(nullptr),
    fHistPrimEtaBothPhInAccPt(nullptr),
    fHistPrimEtaOnePhInAccPt(nullptr),
    fHistPrimEtaPtvsPhInAccPt(nullptr),
    fHistPrimEtaPt(nullptr),
    fHistPrimEtaPtW(nullptr),
    fHistPrimGammaPt(nullptr),
    fHistPrimGammaPtW(nullptr),
    fHistPrimGammaInAccPt(nullptr),
    fHistPrimGammaInAccPtW(nullptr)
{

  if (fOutputContainer)
    delete fOutputContainer;
  fOutputContainer = new THashList();

  for (Int_t i = 0; i < kVtxBins; i++)
    for (Int_t j = 0; j < kCentBins; j++)
      for (Int_t k = 0; k < kPRBins; k++)
        fPHOSEvents[i][j][k] = nullptr; // Container for PHOS photons

  for (Int_t imod = 0; imod < kMods; imod++) {

    fPHOSModEnCorr[imod] = 1.;

    fHistTOFClustMod[imod] = nullptr;
    fHistClustMultVsCentralityMod[imod] = nullptr;
    fHistCaloPhotonTOFvsEMod[imod] = nullptr;
    fHistClustFullEMod[imod] = nullptr;
    fHistClustCoreEMod[imod] = nullptr;
    fHistCaloPhotonPtMod[imod] = nullptr;
    fHistCaloPhotonPtModW[imod] = nullptr;

    fHistMggMod[imod] = nullptr;
    fHistMixMggMod[imod] = nullptr;
    fHistMggPhIDCutMod[imod] = nullptr;
    fHistMixMggPhIDCutMod[imod] = nullptr;

    fHistClustEvsXZMod[imod] = nullptr;
  }

  for (Int_t icut = 0; icut < kPIDCuts; icut++) {
    fHistCaloPhotonPt[icut] = nullptr;
    fHistCaloPhoton2Pt[icut] = nullptr;
    fHistCaloPhotonPtW[icut] = nullptr;
    fHistMCCaloPartIDvsPt[icut] = nullptr;
    fHistMCCaloPartIDvsPtW[icut] = nullptr;

    fHistMgg[icut] = nullptr;
    fHistMgg2[icut] = nullptr;
    fHistMggW[icut] = nullptr;
    fHistMixMgg[icut] = nullptr;
    fHistMixMgg2[icut] = nullptr;
    fHistMixMggW[icut] = nullptr;

    fHistMggCutEff[icut] = nullptr;
    fHistMixMggCutEff[icut] = nullptr;
  }

  for (Int_t imap = 0; imap < 6; imap++)
    fPHOSBadMap[imap] = nullptr;

  for(Int_t i = 0; i < 20; i++){
    fPi0PtWeight[i]   = nullptr;
    fEtaPtWeight[i]   = nullptr;
    fK0SPtWeight[i]   = nullptr;
    fGammaPtWeight[i] = nullptr;

    fPi0PtWeight[i] = new TF1(Form("fPi0PtWeight_%d",i) ,"1.", 0, 200);
    fEtaPtWeight[i] = new TF1(Form("fEtaPtWeight_%d",i) ,"1.", 0, 200);
    fK0SPtWeight[i] = new TF1(Form("fK0SPtWeight_%d",i) ,"1.", 0, 200);
    fGammaPtWeight[i] = new TF1(Form("fGammaPtWeight_%d",i) ,"1.", 0, 200);
  }

  fUserNonlinFunc = new TF1("fUserNonlinFunc", "1.", 0, 200);

}
//_____________________________________________________________________________
AliAnalysisPHOSNeutralMesonsAndPhotons::~AliAnalysisPHOSNeutralMesonsAndPhotons()
{
  // Destructor
  if (fOutputContainer && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    fOutputContainer->Clear();
    delete fOutputContainer;
    fOutputContainer = nullptr;
  }

  if (fCaloPhotonsPHOS){
    fCaloPhotonsPHOS->Clear();
    delete fCaloPhotonsPHOS;
    fCaloPhotonsPHOS = nullptr;
  }

  for (Int_t i = 0; i < kVtxBins; i++)
    for (Int_t j = 0; j < kCentBins; j++)
      for (Int_t k = 0; k < kPRBins; k++)
        if (fPHOSEvents[i][j][k]) {
          fPHOSEvents[i][j][k]->Clear();
          delete fPHOSEvents[i][j][k];
          fPHOSEvents[i][j][k] = nullptr;
        }

    if(fPHOSClusterCuts){
    delete fPHOSClusterCuts;
    fPHOSClusterCuts = nullptr;
  }

  if (fUserNonlinFunc)
    delete fUserNonlinFunc;
  fUserNonlinFunc = nullptr;

  for(Int_t i = 0; i < 20; i++){
    if(fPi0PtWeight[i]){
      delete fPi0PtWeight[i];
      fPi0PtWeight[i] = nullptr;
    }

    if(fEtaPtWeight[i]){
      delete fEtaPtWeight[i];
      fEtaPtWeight[i] = nullptr;
    }

    if(fK0SPtWeight[i]){
      delete fK0SPtWeight[i];
      fK0SPtWeight[i] = nullptr;
    }

    if(fGammaPtWeight[i]){
      delete fGammaPtWeight[i];
      fGammaPtWeight[i] = nullptr;
    }
  }

  if(fCentArrayPi0){
    delete fCentArrayPi0;
    fCentArrayPi0 = nullptr;
  }
  if(fCentArrayEta){
    delete fCentArrayEta;
    fCentArrayEta = nullptr;
  }
  if(fCentArrayGamma){
    delete fCentArrayGamma;
    fCentArrayGamma = nullptr;
  }
  if(fCentArrayK0S){
    delete fCentArrayK0S;
    fCentArrayK0S = nullptr;
  }
  // No need to delete histograms in array fhHistos[]!
  // They are deleted as content of fOutputContainer
}
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::Terminate(Option_t*)
{
}
//_____________________________________________________________________________
AliAnalysisPHOSNeutralMesonsAndPhotons& AliAnalysisPHOSNeutralMesonsAndPhotons::operator=(const AliAnalysisPHOSNeutralMesonsAndPhotons& ref)
{
  // assignment operator

  this->~AliAnalysisPHOSNeutralMesonsAndPhotons();
  new (this) AliAnalysisPHOSNeutralMesonsAndPhotons(ref);
  return *this;
}
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::UserCreateOutputObjects()
{

  if (fOutputContainer)
    delete fOutputContainer;

  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  const Int_t nPt = 99;
  const Double_t Pt[nPt] = { 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, // 21
                             2.2, 2.4, 2.5, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 4.8, 5.0, 5.5, 5.6, 6.0, 6.4, 6.5, 7.0, 7.2, // 42
                             7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 22.0, 24.0,  // 60
                             25.0, 26.0, 28.0, 30.0, 32.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75, 80, 85, 90, 95, 100,  // 79
                             105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200 };   // 99

  const Int_t nMgg = 180;
  const Double_t MggMax = 0.72;

  const Char_t* CutNames[kPIDCuts] = { "All", "DispCut", "CPVCut", "PhIDCut" };
  const Char_t* CutEffNames[kPIDCuts] = { "Probe", "DispCut", "CPVCut", "PhIDCut" };

  fCaloTriggerMimicHelper = new AliCaloTriggerMimicHelper("PHOSTrigAnalysis", 2, fIsMC);
  fCaloTriggerMimicHelper->SetPHOSTrigger(fPHOSTrigType);
  if (fCaloTriggerMimicHelper) {
    fOutputContainer->Add(fCaloTriggerMimicHelper->GetTriggerMimicHelperHistograms());
  }

  fHistInfo = new TH1F("hInfo", ";Option;Value", 20, 0, 20);
  fHistInfo->GetXaxis()->SetBinLabel(1, "CollSys"); // pp - 1, PbPb - 2, pPb - 3
  fHistInfo->GetXaxis()->SetBinLabel(2, "CollEn");
  fHistInfo->GetXaxis()->SetBinLabel(3, "isMC"); // 1 - true, 2 - false
  fHistInfo->GetXaxis()->SetBinLabel(4, "CentMin");
  fHistInfo->GetXaxis()->SetBinLabel(5, "CentMax");
  fHistInfo->GetXaxis()->SetBinLabel(6, "UseCoreE");      // 1 - true, 2 - false
  fHistInfo->GetXaxis()->SetBinLabel(7, "UseTOFCorr");    // 1 - true, 2 - false
  fHistInfo->GetXaxis()->SetBinLabel(8, "UseNonlinCorr"); // 1 - true, 2 - false
  fHistInfo->GetXaxis()->SetBinLabel(9, "EnMinCut");
  fHistInfo->GetXaxis()->SetBinLabel(10, "TOF");
  fHistInfo->GetXaxis()->SetBinLabel(11, "CPV");
  fHistInfo->GetXaxis()->SetBinLabel(12, "Disp");

  fHistInfo->SetBinContent(1, 2);
  fHistInfo->SetBinContent(2, 5);
  fHistInfo->SetBinContent(3, (fIsMC) ? 1 : 0);
  fHistInfo->SetBinContent(4, (Double_t)fCentralityMin);
  fHistInfo->SetBinContent(5, (Double_t)fCentralityMax);
  fHistInfo->SetBinContent(6, (fUseCoreEnergy) ? 1 : 0);
  fHistInfo->SetBinContent(7, (fDoTOFEffCorr) ? 1 : 0);
  fHistInfo->SetBinContent(8, (fDoNonlinCorr) ? 1 : 0);
  fHistInfo->SetBinContent(9, fEminCut);
  fHistInfo->SetBinContent(10, fTOFCut);
  fHistInfo->SetBinContent(11, fCPVSigma);
  fHistInfo->SetBinContent(12, fDispSigma);

  fOutputContainer->Add(fHistInfo);

  fHistSelectEvent = new TH1F("hEventSummary", ";Step;Number of events", 10, 1, 10);
  fHistSelectEvent->GetXaxis()->SetBinLabel(1, "All Events");
  fHistSelectEvent->GetXaxis()->SetBinLabel(2, "Trig. events");
  fHistSelectEvent->GetXaxis()->SetBinLabel(3, "|Z_{vtx}| < 10 cm");
  fHistSelectEvent->GetXaxis()->SetBinLabel(4, "Pileup");
  fHistSelectEvent->GetXaxis()->SetBinLabel(5, "N_{contr} > 1");
  fHistSelectEvent->GetXaxis()->SetBinLabel(6, "Centrality");

  fOutputContainer->Add(fHistSelectEvent);

  fHistVertexZ = new TH1F("hVertexZ", ";Z, cm;Counts", 100, -50., 50.);
  fOutputContainer->Add(fHistVertexZ);

  fHistCentMain = new TH1F(Form("hCentrality%s", fCentralityEstimator.Data()),
                           Form("%s Estimator;%s %%;Counts", fCentralityEstimator.Data(), fCentralityEstimator.Data()), 100, 0., 100);

  fHistCentV0MvsCL0 = new TH2F("hCentralityV0MvsCL0", "Centrality V0M vs. CL0;V0M;CL0", 100, 0., 100, 100, 0., 100.);
  fHistCentV0MvsCL1 = new TH2F("hCentralityV0MvsCL1", "Centrality V0M vs. CL1;V0M;CL1", 100, 0., 100, 100, 0., 100.);
  fHistCentV0MvsV0C = new TH2F("hCentralityV0AvsV0C", "Centrality V0A vs. V0C;V0A;V0C", 100, 0., 100, 100, 0., 100.);
  fHistCentCL0vsCL1 = new TH2F("hCentralityCL0vsCL1", "Centrality CL0 vs. CL1;CL0;CL1", 100, 0., 100, 100, 0., 100.);
  fHistCentZNAvsZNC = new TH2F("hCentralityV0AvsV0C", "Centrality V0A vs. V0C;V0A;V0C", 100, 0., 100, 100, 0., 100.);

  fOutputContainer->Add(fHistCentMain);
  fOutputContainer->Add(fHistCentV0MvsCL0);
  fOutputContainer->Add(fHistCentV0MvsCL1);
  fOutputContainer->Add(fHistCentV0MvsV0C);
  fOutputContainer->Add(fHistCentCL0vsCL1);
  fOutputContainer->Add(fHistCentZNAvsZNC);

  fHistClustMultVsCentrality = new TH2F(Form("hPHOSClustMultVs%sCentrality", fCentralityEstimator.Data()),
                                        ";Centrality %;PHOS clusters multiplicity", 100, 0., 100., 200, 0., 200.);
  fOutputContainer->Add(fHistClustMultVsCentrality);

  for (Int_t imod = 1; imod < kMods + 1; imod++) {
    fHistClustMultVsCentralityMod[imod - 1] = new TH2F(Form("hPHOSClustMultVs%sCentrality_Mod%d", fCentralityEstimator.Data(), imod),
                                                       ";Centrality %;PHOS clusters multiplicity", 100, 0., 100., 200, 0., 200.);
    fOutputContainer->Add(fHistClustMultVsCentralityMod[imod - 1]);
  }

  fHistTOFClust = new TH2F("hClustTOFvsE", "Clust TOF vs E;E, GeV;TOF, ns", nPt - 1, Pt, 1000, -500, 500);
  fOutputContainer->Add(fHistTOFClust);

  for (Int_t imod = 1; imod < kMods + 1; imod++) {
    fHistTOFClustMod[imod - 1] = new TH2F(Form("hClustTOFvsEnMod%d", imod),
                                          Form("Clust TOF vs E mod. %d;E, GeV;TOF, ns", imod),
                                          nPt - 1, Pt, 1000, -500, 500);
    fOutputContainer->Add(fHistTOFClustMod[imod - 1]);
  }

  if (fDoClustQA) {
    fHistCaloPhotonTOFvsE = new TH2F("hCaloPhotonTOFvsEn", "CaloPhoton TOF vs E;E, GeV;TOF, ns", nPt - 1, Pt, 1000, -500, 500);
    fOutputContainer->Add(fHistCaloPhotonTOFvsE);
    for (Int_t imod = 1; imod < kMods + 1; imod++) {
      fHistCaloPhotonTOFvsEMod[imod - 1] = new TH2F(Form("hCaloPhotonTOFvsEnMod%d", imod),
                                                    "CaloPhoton TOF vs E;E, GeV;TOF, ns", nPt - 1, Pt, 1000, -500, 500);
      fOutputContainer->Add(fHistCaloPhotonTOFvsEMod[imod - 1]);
    }

    fHistCaloPhotonTOFCut = new TH1F("hCaloPhotonTOFCut", "CaloPhoton TOF vs E;E, GeV;TOF, ns", 120, -300, 300);
    fHistClustTOFvsDDL = new TH2F("hClustTOFvsDDL", "TOF_{clust}, all clust;DDL;TOF_{clust}, ns", 20, 1, 21, 600, -300, 300);
    fHistClustTOFvsDDLEnCut = new TH2F("hClustTOFvsDDLEnCut", "TOF_{clust}, E > 1.5 GeV;DDL;TOF_{clust}, ns", 20, 1, 21, 600, -300, 300);

    fOutputContainer->Add(fHistCaloPhotonTOFCut);
    fOutputContainer->Add(fHistClustTOFvsDDL);
    fOutputContainer->Add(fHistClustTOFvsDDLEnCut);
  }

  for (Int_t icut = 0; icut < kPIDCuts; icut++) {
    fHistCaloPhotonPt[icut] = new TH1F(Form("hCaloPhotonPt%s", CutNames[icut]), ";#it{p}_{T}, GeV/c;Counts", nPt - 1, Pt);
    fHistCaloPhotonPtW[icut] = new TH1F(Form("hCaloPhotonPtW%s", CutNames[icut]), ";#it{p}_{T}, GeV/c;Counts", nPt - 1, Pt);
    fOutputContainer->Add(fHistCaloPhotonPt[icut]);
    fOutputContainer->Add(fHistCaloPhotonPtW[icut]);

    fHistCaloPhoton2Pt[icut] = new TH1F(Form("hCaloPhoton2Pt%s", CutNames[icut]), ";#it{p}_{T}, GeV/c;Counts", nPt - 1, Pt);
    fOutputContainer->Add(fHistCaloPhoton2Pt[icut]);

  }

  fHistM02vsPt = new TH2F("hClustM02vsPt", ";#it{p}_{T}, GeV/c;Counts", nPt - 1, Pt, 100, 0., 10.);
  fHistM20vsPt = new TH2F("hClustM20vsPt", ";#it{p}_{T}, GeV/c;Counts", nPt - 1, Pt, 100, 0., 10.);
  for (Int_t imod = 0; imod < kMods; imod++){
    fHistClustEvsXZMod[imod] = new TH2F(Form("hClustEvsXZMod%d",imod+1), Form("Cluster E(X,Z) M%d;X cells;Z cells", imod+1), 
                                        64,0.5,64.5, 56,0.5,56.5);
    fOutputContainer->Add(fHistClustEvsXZMod[imod]);                      
  }

  fOutputContainer->Add(fHistM02vsPt);
  fOutputContainer->Add(fHistM20vsPt);

  fHistClustFullE = new TH1F("hClustFullEn", "Full energy;E, GeV;Counts", nPt - 1, Pt);
  fHistClustCoreE = new TH1F("hClustCoreEn", "Core energy;E, GeV;Counts", nPt - 1, Pt);

  fOutputContainer->Add(fHistClustFullE);
  fOutputContainer->Add(fHistClustCoreE);

  fHistNonlinTest = new TH2F("hNonlinTest", "Clust E nonlineatity test;E, GeV;Nonlin. modificator", nPt - 1, Pt, 100, 0.95, 1.05);
  fOutputContainer->Add(fHistNonlinTest);

  if (fDoClustQA) {

    fHistMggTOFCutEffBase = new TH2F("hMggTOFCutEffBase",
                                     ";M_{#gamma#gamma}, GeV/c^{2}; E_{clust}, GeV", nMgg, 0., MggMax, nPt - 1, Pt);

    fHistMggTOFCutEffProbe = new TH2F("hMggTOFCutEffProbe",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; E_{clust}, GeV", nMgg, 0., MggMax, nPt - 1, Pt);

    fHistMixMggTOFCutEffBase = new TH2F("hMixMggTOFCutEffBase",
                                        ";M_{#gamma#gamma}, GeV/c^{2}; E_{clust}, GeV", nMgg, 0., MggMax, nPt - 1, Pt);

    fHistMixMggTOFCutEffProbe = new TH2F("hMixMggTOFCutEffProbe",
                                         ";M_{#gamma#gamma}, GeV/c^{2}; E_{clust}, GeV", nMgg, 0., MggMax, nPt - 1, Pt);

    fOutputContainer->Add(fHistMggTOFCutEffBase);
    fOutputContainer->Add(fHistMixMggTOFCutEffBase);
    fOutputContainer->Add(fHistMggTOFCutEffProbe);
    fOutputContainer->Add(fHistMixMggTOFCutEffProbe);

    for (Int_t imod = 1; imod < kMods + 1; imod++) {
      fHistClustFullEMod[imod - 1] = new TH1F(Form("hClustFullEnMod%d", imod),
                                              "Full energy;E, GeV;Counts", nPt - 1, Pt);
      fHistClustCoreEMod[imod - 1] = new TH1F(Form("hClustCoreEnMod%d", imod),
                                              "Core energy;E, GeV;Counts", nPt - 1, Pt);
      fHistCaloPhotonPtMod[imod - 1] = new TH1F(Form("hCaloPhPtMod%d", imod),
                                                ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
      fHistCaloPhotonPtModW[imod - 1] = new TH1F(Form("hCaloPhPtModW%d", imod),
                                                ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
      fHistMggMod[imod - 1] = new TH2F(Form("hMggAllMod%d%d", imod, imod),
                                       ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
      fHistMixMggMod[imod - 1] = new TH2F(Form("hMixMggAllMod%d%d", imod, imod),
                                          ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
      fHistMggPhIDCutMod[imod - 1] = new TH2F(Form("hMggPhIDCutMod%d%d", imod, imod),
                                              ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
      fHistMixMggPhIDCutMod[imod - 1] = new TH2F(Form("hMixMggPhIDCutMod%d%d", imod, imod),
                                                 ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);

      for (Int_t icut = 0; icut < kPIDCuts; icut++) {
        fHistMggCutEffMod[imod - 1][icut] = new TH2F(Form("hMgg%sEffMod%d%d", CutEffNames[icut], imod, imod),
                                                     ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}^{#gamma}", nMgg, 0., MggMax, nPt - 1, Pt);
        fHistMixMggCutEffMod[imod - 1][icut] = new TH2F(Form("hMixMgg%sEffMod%d%d", CutEffNames[icut], imod, imod),
                                                        ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}^{#gamma}", nMgg, 0., MggMax, nPt - 1, Pt);
      }
    }

    for (Int_t imod = 1; imod < kMods + 1; imod++)
      fOutputContainer->Add(fHistClustFullEMod[imod - 1]);
    for (Int_t imod = 1; imod < kMods + 1; imod++)
      fOutputContainer->Add(fHistClustCoreEMod[imod - 1]);
    for (Int_t imod = 1; imod < kMods + 1; imod++)
      fOutputContainer->Add(fHistCaloPhotonPtMod[imod - 1]);
    for (Int_t imod = 1; imod < kMods + 1; imod++)
      fOutputContainer->Add(fHistCaloPhotonPtModW[imod - 1]);

    for (Int_t imod = 1; imod < kMods + 1; imod++)
      fOutputContainer->Add(fHistMggMod[imod - 1]);
    for (Int_t imod = 1; imod < kMods + 1; imod++)
      fOutputContainer->Add(fHistMixMggMod[imod - 1]);
    for (Int_t imod = 1; imod < kMods + 1; imod++)
      fOutputContainer->Add(fHistMggPhIDCutMod[imod - 1]);
    for (Int_t imod = 1; imod < kMods + 1; imod++)
      fOutputContainer->Add(fHistMixMggPhIDCutMod[imod - 1]);

    for (Int_t imod = 1; imod < kMods + 1; imod++)
      for (Int_t icut = 0; icut < kPIDCuts; icut++)
        fOutputContainer->Add(fHistMggCutEffMod[imod - 1][icut]);

    for (Int_t imod = 1; imod < kMods + 1; imod++)
      for (Int_t icut = 0; icut < kPIDCuts; icut++)
        fOutputContainer->Add(fHistMixMggCutEffMod[imod - 1][icut]);
  }

  for (Int_t icut = 0; icut < kPIDCuts; icut++) {
    fHistMgg[icut] = new TH2F(Form("hMgg%s", CutNames[icut]),
                              ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistMggW[icut] = new TH2F(Form("hMggW%s", CutNames[icut]),
                              ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistMixMgg[icut] = new TH2F(Form("hMixMgg%s", CutNames[icut]),
                                 ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistMixMggW[icut] = new TH2F(Form("hMixMggW%s", CutNames[icut]),
                                 ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);

    fHistMggCutEff[icut] = new TH2F(Form("hMgg%sEff", CutEffNames[icut]),
                                    ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}^{#gamma}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistMixMggCutEff[icut] = new TH2F(Form("hMixMgg%sEff", CutEffNames[icut]),
                                       ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}^{#gamma}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);

    fHistMgg2[icut] = new TH2F(Form("hMgg2%s", CutNames[icut]),
                              ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);

    fHistMixMgg2[icut] = new TH2F(Form("hMixMgg2%s", CutNames[icut]),
                                 ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);                          
  }

  for (Int_t icut = 0; icut < kPIDCuts; icut++){
    fOutputContainer->Add(fHistMgg[icut]);
    fOutputContainer->Add(fHistMggW[icut]);
  }
  for (Int_t icut = 0; icut < kPIDCuts; icut++){
    fOutputContainer->Add(fHistMgg2[icut]);
    fOutputContainer->Add(fHistMixMgg2[icut]);
  }
  for (Int_t icut = 0; icut < kPIDCuts; icut++){
    fOutputContainer->Add(fHistMixMgg[icut]);
    fOutputContainer->Add(fHistMixMggW[icut]);
  }
  for (Int_t icut = 0; icut < kPIDCuts; icut++)
    fOutputContainer->Add(fHistMggCutEff[icut]);
  for (Int_t icut = 0; icut < kPIDCuts; icut++)
    fOutputContainer->Add(fHistMixMggCutEff[icut]);

  if (fIsMC) {

    fHistTruePi0MggVsRecPt = new TH2F("hTruePi0MggVsRecPt",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistTruePi0MggVsTruePt = new TH2F("hTruePi0MggVsTruePt",
                                       ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistTruePi0MotherIDvsRecPt = new TH2F("hTruePi0MotherIDvsRecPt",
                                       ";#it{p}_{T}, GeV/#it{c};Index", nPt - 1, Pt, 20,1,20);
    fHistTruePi0MotherIDvsRecPtW = new TH2F("hTruePi0MotherIDvsRecPtW",
                                       ";#it{p}_{T}, GeV/#it{c};Index", nPt - 1, Pt, 20,1,20);

    fHistTruePi0MotherIDvsTruePt = new TH2F("hTruePi0MotherIDvsTruePt",
                                       ";#it{p}_{T}, GeV/#it{c};Index", nPt - 1, Pt, 20,1,20);                                   

    fHistTrueEtaMggVsRecPt = new TH2F("hTrueEtaMggVsRecPt",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistTrueEtaMggVsTruePt = new TH2F("hTrueEtaMggVsTruePt",
                                       ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);

    std::map<Int_t, TString> pdgMapSecPi0 = { { 1, "K_{0}^{S}" }, { 2, "#Lambda" }, { 3, "K_{0}^{L}" }, 
                                              { 4, "p^{#pm}" }, { 5, "n(#bar{n})" }, { 6, "#pi^{#pm}" }, 
                                              { 7, "K^{#pm}" }, { 8, "#rho^{0,#pm}" }, { 9, "#Sigma" }, 
                                              { 10, "#Delta" }, { 11, "K^{*}" }, { 12, "#pi^{0}"}, { 15, "Rest" } };

    for (auto& ibin : pdgMapSecPi0) {
      if (ibin.first < 0 || ibin.first > fHistTruePi0MotherIDvsRecPt->GetNbinsY())
        continue;
      fHistTruePi0MotherIDvsRecPt->GetYaxis()->SetBinLabel(ibin.first, ibin.second);
      fHistTruePi0MotherIDvsRecPtW->GetYaxis()->SetBinLabel(ibin.first, ibin.second);
      fHistTruePi0MotherIDvsTruePt->GetYaxis()->SetBinLabel(ibin.first, ibin.second);
    }

    fOutputContainer->Add(fHistTruePi0MggVsRecPt);
    fOutputContainer->Add(fHistTruePi0MggVsTruePt);
    fOutputContainer->Add(fHistTruePi0MotherIDvsRecPt);
    fOutputContainer->Add(fHistTruePi0MotherIDvsRecPtW);
    fOutputContainer->Add(fHistTruePi0MotherIDvsTruePt);
    fOutputContainer->Add(fHistTrueEtaMggVsRecPt);
    fOutputContainer->Add(fHistTrueEtaMggVsTruePt);

    fHistTruePi0MggVsRecPtFromK0s = new TH2F("hTruePi0MggVsRecPtFromK0s",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistTruePi0MggVsRecPtFromK0L = new TH2F("hTruePi0MggVsRecPtFromK0L",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistTruePi0MggVsRecPtFromLambda = new TH2F("hTruePi0MggVsRecPtFromLambda",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistTruePi0MggVsRecPtFromPandN = new TH2F("hTruePi0MggVsRecPtFromPandN",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistTruePi0MggVsRecPtFromPions = new TH2F("hTruePi0MggVsRecPtFromPions",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);
    fHistTruePi0MggVsRecPtFromKaons = new TH2F("hTruePi0MggVsRecPtFromKaons",
                                      ";M_{#gamma#gamma}, GeV/c^{2}; #it{p}_{T}, GeV/#it{c}", nMgg, 0., MggMax, nPt - 1, Pt);

    fOutputContainer->Add(fHistTruePi0MggVsRecPtFromK0s);
    fOutputContainer->Add(fHistTruePi0MggVsRecPtFromK0L);
    fOutputContainer->Add(fHistTruePi0MggVsRecPtFromLambda);
    fOutputContainer->Add(fHistTruePi0MggVsRecPtFromPandN);
    fOutputContainer->Add(fHistTruePi0MggVsRecPtFromPions);
    fOutputContainer->Add(fHistTruePi0MggVsRecPtFromKaons);                            
                                      
    std::map<Int_t, TString> pdgMap = { { 1, "K_{0}^{S}" }, { 2, "#Lambda" }, { 3, "K_{0}^{L}" }, 
                                        { 4, "p^{#pm}" }, { 5, "n(#bar{n})" }, { 6, "#pi^{#pm}" }, 
                                        { 7, "K^{#pm}" }, { 8, "#rho^{0,#pm}" }, { 9, "#Sigma" }, 
                                        { 10, "#Delta" }, { 11, "K^{*}" }, { 15, "Rest" } };

    fHistMCPartIDvsPt = new TH2F("hMCPartIDvsPt", ";#it{p}_{T}, GeV/#it{c};Index", nPt - 1, Pt, 20, 1, 20);
    fHistMCPartIDvsPtW = new TH2F("hMCPartIDvsPtW", ";#it{p}_{T}, GeV/#it{c};Index", nPt - 1, Pt, 20, 1, 20);
    for (auto& ibin : pdgMap) {
      if (ibin.first < 0 || ibin.first > fHistMCPartIDvsPt->GetNbinsY())
        continue;
      fHistMCPartIDvsPt->GetYaxis()->SetBinLabel(ibin.first, ibin.second);
      fHistMCPartIDvsPtW->GetYaxis()->SetBinLabel(ibin.first, ibin.second);
    }
    fOutputContainer->Add(fHistMCPartIDvsPt);
    fOutputContainer->Add(fHistMCPartIDvsPtW);

    std::map<Int_t, TString> pdgCaloMap = { { 1, "#gamma" }, { 2, "e^{#pm}" }, { 3, "#pi^{#pm}" }, 
                                            { 4, "K^{#pm}" }, { 5, "K_{0}^{L}" }, { 6, "p^{+}" }, 
                                            { 7, "p^{-}" }, { 8, "n" }, { 9, "#bar{n}" }, { 10, "Rest" } };

    for (Int_t icut = 0; icut < kPIDCuts; icut++) {
      fHistMCCaloPartIDvsPt[icut] = new TH2F(Form("hMCCaloPartIDvsPt%s", CutNames[icut]),
                                             ";#it{p}_{T}, GeV/#it{c};Index", nPt - 1, Pt, 20, 1, 21);
      fHistMCCaloPartIDvsPtW[icut] = new TH2F(Form("hMCCaloPartIDvsPtW%s", CutNames[icut]),
                                             ";#it{p}_{T}, GeV/#it{c};Index", nPt - 1, Pt, 20, 1, 21);
      for (auto& ibin : pdgCaloMap) {
        if (ibin.first < 0 || ibin.first > fHistMCCaloPartIDvsPt[icut]->GetNbinsY())
          continue;
        fHistMCCaloPartIDvsPt[icut]->GetYaxis()->SetBinLabel(ibin.first, ibin.second);
        fHistMCCaloPartIDvsPtW[icut]->GetYaxis()->SetBinLabel(ibin.first, ibin.second);
      }
      fOutputContainer->Add(fHistMCCaloPartIDvsPt[icut]);
      fOutputContainer->Add(fHistMCCaloPartIDvsPtW[icut]);
    }

    fHistPrimPi0Pt = new TH1F("hPrimPi0Pt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimPi0PtW = new TH1F("hPrimPi0PtW", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimPi0InAccPt = new TH1F("hPrimPi0InAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimPi0InAccPtW = new TH1F("hPrimPi0InAccPtW", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimPi0BothPhInAccPt = new TH1F("hPrimPi0BothPhInAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimPi0BothPhInAccPtW = new TH1F("hPrimPi0BothPhInAccPtW", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimPi0OnePhInAccPt = new TH1F("hPrimPi0OnePhInAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimPi0PtvsPhInAccPt = new TH2F("hPrimPi0PtvsPhInAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt, 20, 0., 1.);

    fHistPrimEtaPt = new TH1F("hPrimEtaPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimEtaPtW = new TH1F("hPrimEtaPtW", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimEtaInAccPt = new TH1F("hPrimEtaInAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimEtaBothPhInAccPt = new TH1F("hPrimEtaBothPhInAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimEtaOnePhInAccPt = new TH1F("hPrimEtaOnePhInAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimEtaPtvsPhInAccPt = new TH2F("hPrimEtaPtvsPhInAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt, 20, 0., 1.);

    fHistPrimGammaPt = new TH1F("hPrimGammaPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimGammaPtW = new TH1F("hPrimGammaPtW", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimGammaInAccPt = new TH1F("hPrimGammaInAccPt", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);
    fHistPrimGammaInAccPtW = new TH1F("hPrimGammaInAccPtW", ";#it{p}_{T}, GeV/#it{c};Counts", nPt - 1, Pt);

    fOutputContainer->Add(fHistPrimPi0Pt);
    fOutputContainer->Add(fHistPrimPi0PtW);
    fOutputContainer->Add(fHistPrimPi0InAccPt);
    fOutputContainer->Add(fHistPrimPi0InAccPtW);
    fOutputContainer->Add(fHistPrimPi0BothPhInAccPt);
    fOutputContainer->Add(fHistPrimPi0BothPhInAccPtW);
    fOutputContainer->Add(fHistPrimPi0OnePhInAccPt);
    fOutputContainer->Add(fHistPrimPi0PtvsPhInAccPt);

    fOutputContainer->Add(fHistPrimEtaPt);
    fOutputContainer->Add(fHistPrimEtaPtW);
    fOutputContainer->Add(fHistPrimEtaInAccPt);
    fOutputContainer->Add(fHistPrimEtaBothPhInAccPt);
    fOutputContainer->Add(fHistPrimEtaOnePhInAccPt);
    fOutputContainer->Add(fHistPrimEtaPtvsPhInAccPt);

    fOutputContainer->Add(fHistPrimGammaPt);
    fOutputContainer->Add(fHistPrimGammaPtW);
    fOutputContainer->Add(fHistPrimGammaInAccPt);
    fOutputContainer->Add(fHistPrimGammaInAccPtW);

  }

  PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::UserExec(Option_t*)
{

  // AliInfo("Event has been started!");

  fEvent = (AliAODEvent*)InputEvent();
  fMCStack = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());

  fHistSelectEvent->Fill(1); // all events
  if (!fEvent) {
    AliInfo("ERROR: Could not retrieve event");
    return;
  }

  if (fIsMC && !fMCStack) {
    AliInfo("Could not retrieve MC stack for MC study");
    return;
  }

  if(!fPHOSClusterCuts){
    AliError("fPHOSClusterCuts is not set! return");
    return;
  }

  fRunNumber = fEvent->GetRunNumber();

  if (!fAODCells)
    fAODCells = fEvent->GetPHOSCells();

  TString trigClasses = fEvent->GetFiredTriggerClasses();
  //printf("%s\n",trigClasses.Data());
  if(!fIsMC 
      && !trigClasses.Contains("-CENT") //accept CENT, CENTNO[TRD|PMD]
      && !trigClasses.Contains("-FAST") //accept FAST, not MUFAST
      && !trigClasses.Contains("-CALO") //accept CALO, CALOFAST
    ){
    AliWarning(Form("Skip event with triggers %s",trigClasses.Data()));
    return;
  }

  UInt_t SelectMask = fInputHandler->IsEventSelected();
  // UInt_t MBTriggerMask = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB;
  // UInt_t MBTriggerMask = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kMB; //test
  UInt_t MBTriggerMask = GetCollisionCandidates();

  Bool_t isMBTriggerSelected = SelectMask & MBTriggerMask;
  Bool_t isPHOSTriggerSelected = SelectMask & AliVEvent::kPHI7;
  Bool_t EventIsSelected = kFALSE;

  if (fIsMBTriggerAnalysis) {
    if (isMBTriggerSelected)
      EventIsSelected = kTRUE;
  }

  if (!EventIsSelected) {
    if (fIsPHOSTriggerAnalysis)
      AliInfo("MB is not triggered.");
    else
      AliInfo("PHI7 is not triggered.");
    return;
  }

  if (fIsMBTriggerAnalysis) { // test
    if (!isMBTriggerSelected) {
      AliInfo("MB is not triggered.");
      return;
    }
  }

  fHistSelectEvent->Fill(2); // triggered events

  const AliAODVertex* vVertex = fEvent->GetPrimaryVertex();
  fVertex[0] = vVertex->GetX();
  fVertex[1] = vVertex->GetY();
  fVertex[2] = vVertex->GetZ();
  fZvtx = (Int_t)((fVertex[2] + 10.) / 1.); // it should be 0-kVtxBins-1.
  if (fZvtx < 0)
    fZvtx = 0; // protection to avoid fZvtx = -1.
  if (fZvtx > kVtxBins)
    fZvtx = kVtxBins - 1; // protection to avoid fZvtx = kVtxBins

  fHistVertexZ->Fill(fVertex[2]);

  if (fVertex[2] > 10.) {
    AliInfo("Event rejected Zvrtx > 10. cm");
    return;
  }

  fHistSelectEvent->Fill(3); // events in Zvrtx < 10. cm

  if (!fUtils)
    fUtils = new AliAnalysisUtils();

  if (fEvent->IsPileupFromSPD()) {
    AliInfo("Event rejected due to pileup");
    return;
  }

  AliPHOSEventCuts* fPHOSEventCuts = new AliPHOSEventCuts("PHOSEventCuts");
  fPHOSEventCuts->SetMCFlag(fIsMC);
  fPHOSEventCuts->SetMaxAbsZvtx(10.);
  fPHOSEventCuts->SetRejectPileup(kTRUE);
  fPHOSEventCuts->SetRejectDAQIncompleteEvent(kTRUE);
  fPHOSEventCuts->SetPileupFinder(AliPHOSEventCuts::kMultiVertexer);

  if(!(fPHOSEventCuts->AcceptEvent(fEvent))){
    AliInfo("event is rejected.");
    return;
  }

  // Pileup in MC
  // if (fIsMC){
  //   AliAODMCHeader *aodMCHeader = (AliAODMCHeader *)fEvent->FindListObject(AliAODMCHeader::StdBranchName());
  //   if (aodMCHeader){
  //     // find cocktail header
  //     Int_t nGenerators = aodMCHeader->GetNCocktailHeaders();
  //     if (nGenerators > 0){
  //       for (Int_t igen = 0; igen < nGenerators; igen++){
  //         AliGenEventHeader *eventHeaderGen = aodMCHeader->GetCocktailHeader(igen);
  //         TString genname = eventHeaderGen->ClassName();
  //         bool isPileUp = AliAnalysisUtils::IsPileupInGeneratedEvent(aodMCHeader, genname);
  //         if (isPileUp)
  //           return;
  //       }
  //     }
  //   }
  // }

  fHistSelectEvent->Fill(4); // events w/o pileup

  Int_t nPrimContributors = vVertex->GetNContributors();
  if (!fIsMC && (nPrimContributors < 1)) {
    return;
  }

  fHistSelectEvent->Fill(5); // events NContributors > 1

  SelectCentrality();
  if (fCentrality > -1)
    if (fCentrality < fCentralityMin || fCentrality > fCentralityMax) {
      AliInfo(Form("Reject this event because centrality %s %f %% is out of the configuration of this task.", fCentralityEstimator.Data(), fCentrality));
      return;
    }

  fHistSelectEvent->Fill(6); // events NContributors > 1

  // centrality bins for mixed events container
  fCentBin = 0;
  if (fCentrality < 5.)
    fCentBin = 0;
  else if (fCentrality < 10.)
    fCentBin = 1;
  else if (fCentrality < 20.)
    fCentBin = 2;
  else if (fCentrality < 30.)
    fCentBin = 3;
  else if (fCentrality < 40.)
    fCentBin = 4;
  else if (fCentrality < 50.)
    fCentBin = 5;
  else if (fCentrality < 80.)
    fCentBin = 6;

  if (!fPHOSGeo) {
    InitPHOSGeometry();
  }

  // reaction plane
  AliEventplane* eventPlane = fEvent->GetEventplane();
  if (!eventPlane) { // Event has no event plane
    AliInfo("Reaction plane has not been found!");
    return;
  }
  // V0A
  const Int_t harmonics = 2;
  double qx = 0., qy = 0.;
  double rpV0A =
    eventPlane->CalculateVZEROEventPlane(fEvent, 8, harmonics, qx, qy);
  // V0C
  double rpV0C =
    eventPlane->CalculateVZEROEventPlane(fEvent, 9, harmonics, qx, qy);

  // Whole V0
  fRP = eventPlane->CalculateVZEROEventPlane(fEvent, 10, harmonics, qx, qy);

  while (rpV0A < 0)
    rpV0A += TMath::TwoPi() / harmonics;
  while (rpV0A > TMath::TwoPi() / harmonics)
    rpV0A -= TMath::TwoPi() / harmonics;

  while (rpV0C < 0)
    rpV0C += TMath::TwoPi() / harmonics;
  while (rpV0C > TMath::TwoPi() / harmonics)
    rpV0C -= TMath::TwoPi() / harmonics;

  while (fRP < 0)
    fRP += TMath::TwoPi() / harmonics;
  while (fRP > TMath::TwoPi() / harmonics)
    fRP -= TMath::TwoPi() / harmonics;

  Int_t irp = Int_t(kPRBins * (fRP) / TMath::Pi());
  if (irp < 0)
    irp = 0;
  if (irp >= kPRBins)
    irp = kPRBins - 1;

  if (!fPHOSEvents[fZvtx][fCentBin][irp])
    fPHOSEvents[fZvtx][fCentBin][irp] = new TList();
  fCaloPhotonsPHOSList = fPHOSEvents[fZvtx][fCentBin][irp];

  ProcessCaloPhotons();
  FillMgg();
  EstimatePIDCutEfficiency();
  if (fDoClustQA)
    EstimateTOFCutEfficiency();

  if (fIsMC)
    ProcessMCParticles();

  if (fCaloPhotonsPHOS->GetEntriesFast() > 0) {
    fCaloPhotonsPHOS->Expand(fCaloPhotonsPHOS->GetEntriesFast());
    fCaloPhotonsPHOSList->AddFirst(fCaloPhotonsPHOS);
    fCaloPhotonsPHOS = 0;
    if (fCaloPhotonsPHOSList->GetSize() > fNMixEvents) {
      TClonesArray* tmp = static_cast<TClonesArray*>(fCaloPhotonsPHOSList->Last());
      fCaloPhotonsPHOSList->RemoveLast();
      delete tmp;
    }
  }

  // AliInfo("Event has been processed!");
  PostData(1, fOutputContainer);

} // end of UserExec()
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::ProcessCaloPhotons()
{

  if (!fCaloPhotonsPHOS)
    fCaloPhotonsPHOS = new TClonesArray("AliCaloPhoton");
  else
    fCaloPhotonsPHOS->Clear();

  Int_t multClust = fEvent->GetNumberOfCaloClusters();
  Int_t inPHOS = 0; // counter of caloClusters
  Int_t PHOSMultMod[4] = { 0, 0, 0, 0 };

  Int_t PHOSMultTOF = 0; // counter of caloClusters
  Int_t PHOSMultTOFMod[4] = { 0, 0, 0, 0 };

  Double_t energy = 0.;
  Double_t coreE = 0.;

  Double_t tof = -1000.;
  Int_t ddl = -1;

  Bool_t CPVBit = kFALSE;
  Bool_t DispBit = kFALSE;
  Bool_t TOFBit = kFALSE;

  for (Int_t iclu = 0; iclu < multClust; iclu++) {

    AliAODCaloCluster* clu = fEvent->GetCaloCluster(iclu);
    if (clu->GetType() != AliVCluster::kPHOSNeutral)
      continue;

    Float_t pos[3];
    clu->GetPosition(pos);

    TVector3 global1(pos);
    Int_t relId[4];
    fPHOSGeo->GlobalPos2RelId(global1, relId);
    Int_t mod = relId[0];
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3];

    if (mod < 1 || mod > 4) {
      AliInfo(Form("Wrong module number %d", mod));
      return;
    }

    coreE = clu->GetCoreEnergy();
    energy = clu->E();
    tof = clu->GetTOF();
    ddl = (5 - mod) * 4 + (cellX - 1) / 16;

    TVector3 local;
    fPHOSGeo->Global2Local(local, global1, mod);

    TLorentzVector p1;
    TLorentzVector p1core;

    clu->GetMomentum(p1, fVertex);
    clu->GetMomentum(p1core, fVertex);
    p1core *= coreE / energy;

    if (fDoNonlinCorr) {
      if (fIsMC){
        p1 *= fUserNonlinFunc->Eval(energy);
        p1core *= fUserNonlinFunc->Eval(coreE);
      }
      else{
        p1 *= fPHOSModEnCorr[mod-1];
        p1core *= fPHOSModEnCorr[mod-1];
      }
    }
    fHistNonlinTest->Fill(p1.E(), p1.E() / energy);

    TLorentzVector p0 = (fUseCoreEnergy) ? p1core : p1;

    if (fMinBiasLowEClustCut){

      if (p1.E() > 2. && clu->GetNCells() < 3)
        continue;

      if (p1.E() > 2. && clu->GetM02() < 0.2) // M02 cut as well
        continue;

      if (p1.E() < 0.1) continue;

      if(!PassSTDCut(clu)) continue;

    }
    else{
      if (clu->GetNCells() < 2)
        continue;
      if (clu->GetM02() < 0.2)
        continue;
    }

    if (p0.E() < fEminCut)
      continue;

    CPVBit = clu->GetEmcCpvDistance() > fCPVSigma;
    DispBit = (fUseCoreEnergy) ? (clu->Chi2() < fDispSigma * fDispSigma) : (clu->GetDispersion() < fDispSigma * fDispSigma);
    // TOFBit = (fIsMC) ? kTRUE : (TMath::Abs(tof * 1e+9) < fTOFCut/2.); // no TOF cut for MC
    TOFBit = TMath::Abs(tof * 1e+9) < fTOFCut/2.;

    if (fDoClustQA) {
      fHistTOFClust->Fill(p0.E(), tof * 1e+9);
      fHistTOFClustMod[mod - 1]->Fill(p0.E(), tof * 1e+9);

      fHistClustTOFvsDDL->Fill(ddl, tof * 1e+9);
      if (p0.E() > 1.5)
        fHistClustTOFvsDDLEnCut->Fill(ddl, tof * 1e+9);

      if (CPVBit && DispBit) {
        fHistCaloPhotonTOFvsE->Fill(p0.E(), tof * 1e+9);
        fHistCaloPhotonTOFvsEMod[mod - 1]->Fill(p0.E(), tof * 1e+9);
      }
    }

    AliCaloPhoton* CaloPhoton = new ((*fCaloPhotonsPHOS)[inPHOS])AliCaloPhoton(p0.Px(), p0.Py(), p0.Pz(), p0.E());

    CaloPhoton->SetModule(mod);
    CaloPhoton->SetDistToBad((Int_t)(1. + clu->GetDistanceToBadChannel() / 2.2));
    CaloPhoton->SetDistToBadfp(clu->GetDistanceToBadChannel() / 2.2); //in unit of cells with floating point. 2.2 cm is crystal size
    CaloPhoton->SetBC(iclu);       // reference to CaloCluster
    CaloPhoton->SetTagInfo(0);     // No pi0 partners found so far
    CaloPhoton->SetTagged(kFALSE); // Reconstructed pairs found
    CaloPhoton->SetEMCx(local.X());
    CaloPhoton->SetEMCz(local.Z());
    CaloPhoton->SetCluster(clu);
    CaloPhoton->SetLambdas(clu->GetM20(), clu->GetM02());
    CaloPhoton->SetNCells(clu->GetNCells());
    CaloPhoton->SetMomV2(&p1core);
    CaloPhoton->SetNsigmaFullDisp(TMath::Sqrt(clu->GetDispersion()));
    CaloPhoton->SetNsigmaCoreDisp(TMath::Sqrt(clu->Chi2()));
    CaloPhoton->SetNsigmaCPV(clu->GetEmcCpvDistance());
    CaloPhoton->SetTOFBit(TOFBit);
    CaloPhoton->SetTime(tof); // in ns
    CaloPhoton->SetWeight(1.);

    // PID cuts
    if (fUsePHOSClusterCuts){
      CaloPhoton->SetCPVBit(fPHOSClusterCuts->IsNeutral(CaloPhoton));
      CaloPhoton->SetDispBit(fPHOSClusterCuts->AcceptDisp(CaloPhoton));
      CaloPhoton->SetPhoton(fPHOSClusterCuts->AcceptPhoton(CaloPhoton));
    }
    else{
      CaloPhoton->SetCPVBit(CPVBit);
      CaloPhoton->SetDispBit(DispBit);
      CaloPhoton->SetPhoton(DispBit && CPVBit);
    }

    if (fIsMC) {

      Bool_t sure = kTRUE;
      Int_t primary = FindClusterPrimary(CaloPhoton, sure);
      CaloPhoton->SetPrimary(primary);

      Int_t clulb = clu->GetLabelAt(0);
      AliAODMCParticle *prim = (AliAODMCParticle*)fMCStack->At(clulb);

      TF1 *f1Pi0Weight   = (TF1*)GetAdditionalPi0PtWeightFunction(fCentrality);
      TF1 *f1EtaWeight   = (TF1*)GetAdditionalEtaPtWeightFunction(fCentrality);
      TF1 *f1GammaWeight = (TF1*)GetAdditionalGammaPtWeightFunction(fCentrality);
      TF1 *f1K0SWeight   = (TF1*)GetAdditionalK0SPtWeightFunction(fCentrality);

      Double_t MCWeight = 1.;
      
      Double_t TruePi0Pt = 0.;
      Double_t TrueEtaPt = 0.;
      //Double_t TrueGammaPt = 0.;
      Double_t TrueK0SPt = 0;
      Double_t TrueL0Pt = 0;

      //weight is always defined as relative X/pi0 ratio, and pi0 weight is absolute number.
      if(IsFrom(primary,TruePi0Pt,111)) MCWeight = f1Pi0Weight->Eval(TruePi0Pt);//pi0
      if(IsFrom(primary,TruePi0Pt,211)) MCWeight = f1Pi0Weight->Eval(TruePi0Pt);//pi+/-

      if(IsFrom(primary,TrueK0SPt,310)) MCWeight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);//K0S
      if(IsFrom(primary,TrueK0SPt,130)) MCWeight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);//K0L
      if(IsFrom(primary,TrueK0SPt,321)) MCWeight = f1K0SWeight->Eval(TrueK0SPt) * f1Pi0Weight->Eval(TrueK0SPt);//K+/-

      if(IsFrom(primary,TrueEtaPt,221)) MCWeight = f1EtaWeight->Eval(TrueEtaPt) * f1Pi0Weight->Eval(TrueEtaPt);//eta

      // if ()

      CaloPhoton->SetWeight(MCWeight);

      } else {
        CaloPhoton->SetPrimary(-1);
        CaloPhoton->SetPrimaryAtVertex(-1);
        CaloPhoton->SetWeight(1.);
    }

    inPHOS++;
    // to do IsTrig()

    if (!CaloPhoton->IsTOFOK())
      continue;

    Double_t Weight = CaloPhoton->GetWeight();

    fHistClustEvsXZMod[mod - 1]->Fill(cellX, cellZ, p0.E());

    PHOSMultMod[mod - 1]++;

    // const Char_t* CutNames[kPIDCuts] = { "All", "DispCut", "CPVCut", "PhIDCut" };

    const Bool_t PhotonCutFlag[kPIDCuts] = { kTRUE, // All
                                             CaloPhoton->IsDispOK(), // DispCut
                                             CaloPhoton->IsCPVOK(), // CPVCut
                                             CaloPhoton->IsPhoton() }; // PhIDCut

    if (fDoClustQA) {
      fHistClustFullEMod[mod - 1]->Fill(p1.E(), Weight);
      fHistClustCoreEMod[mod - 1]->Fill(p1core.E(), Weight);
      if (PhotonCutFlag[3]){
        fHistCaloPhotonPtModW[mod - 1]->Fill(CaloPhoton->Pt(), Weight);
        fHistCaloPhotonPtMod[mod - 1]->Fill(CaloPhoton->Pt());
        fHistCaloPhotonTOFCut->Fill(tof * 1e+9);
      }
    }

    for (Int_t icut = 0; icut < kPIDCuts; icut++) {
      if (!PhotonCutFlag[icut])
        continue;

      fHistCaloPhotonPt[icut]->Fill(CaloPhoton->Pt());
      fHistCaloPhotonPtW[icut]->Fill(CaloPhoton->Pt(), Weight);

      if (!fIsMC)
        continue;

      // Int_t leadlb = clu->GetLabelAt(0);
      Bool_t sure = kTRUE;
      Int_t leadlb = FindClusterPrimary(CaloPhoton, sure);
      Int_t pdg = -1;
      Int_t MCIndex = -1;
      Double_t PrimPt = -1.;

      AliAODMCParticle* particle = nullptr;
      Double_t CaloIndex = 10.5;
      if (leadlb > -1) {
        particle = (AliAODMCParticle*)fMCStack->At(leadlb);
        pdg = particle->GetPdgCode();

        if (pdg == 22) CaloIndex = 1.5; // gamma
        else if (TMath::Abs(pdg) == 11){ // electron
          if (MCMotherIsFound(leadlb, PrimPt, 22)) CaloIndex = 1.5; // gamma
          else if (IsPrimaryMCParticle(particle)) CaloIndex = 2.5; // electron
          else  CaloIndex = 11.5; // primary electrons (???)
        }
        else if (TMath::Abs(pdg) == 211) CaloIndex = 3.5; // pi+-
        else if (TMath::Abs(pdg) == 321) CaloIndex = 4.5; // K+-
        else if (TMath::Abs(pdg) == 130) CaloIndex = 5.5; // K0L
        else if (pdg == 2212) CaloIndex = 6.5; // proton
        else if (pdg == -2212) CaloIndex = 7.5; // antiproton
        else if (pdg == 2112) CaloIndex = 8.5; // neutron
        else if (pdg == -2112) CaloIndex = 9.5; // antineutron

        fHistMCCaloPartIDvsPt[icut]->Fill(CaloPhoton->Pt(), CaloIndex);
        fHistMCCaloPartIDvsPtW[icut]->Fill(CaloPhoton->Pt(), CaloIndex, Weight);
      }
      else{
        fHistMCCaloPartIDvsPt[icut]->Fill(CaloPhoton->Pt(), 10.5);
        fHistMCCaloPartIDvsPtW[icut]->Fill(CaloPhoton->Pt(), 10.5, Weight);
      }
    }

    fHistCaloPhoton2Pt[0]->Fill(CaloPhoton->Pt(), Weight);
    if (CaloPhoton->IsDispOK()) fHistCaloPhoton2Pt[1]->Fill(CaloPhoton->Pt(), Weight);
    if (CaloPhoton->IsCPVOK()) fHistCaloPhoton2Pt[2]->Fill(CaloPhoton->Pt(), Weight);
    if (CaloPhoton->IsPhoton()) fHistCaloPhoton2Pt[3]->Fill(CaloPhoton->Pt(), Weight);

    fHistClustFullE->Fill(p1.E(), Weight);
    fHistClustCoreE->Fill(p1core.E(), Weight);

    fHistM02vsPt->Fill(CaloPhoton->Pt(), clu->GetM02(), Weight);
    fHistM20vsPt->Fill(CaloPhoton->Pt(), clu->GetM20(), Weight);
    
  } // loop over clusters

  for (Int_t imod = 1; imod < 5; imod++)
    fHistClustMultVsCentralityMod[imod - 1]->Fill(fCentrality, PHOSMultMod[imod - 1]);

  fHistClustMultVsCentrality->Fill(fCentrality, inPHOS);

} // end of ProcessCaloPhotons()
//_____________________________________________________________________________
// void AliAnalysisPHOSNeutralMesonsAndPhotons::FillCaloPhotons(){

//   const Int_t nCaloPhotons = fCaloPhotonsPHOS->GetEntriesFast();

//   Double_t pT = 0,energy = 0;
//   Double_t phi = -999, dphi = -999.;
//   Int_t primary = -1;
//   Double_t weight = 1.;
//   Double_t TrueK0SPt = 0;
//   Double_t TrueL0Pt = 0;
//   Double_t TrueS0Pt = 0;

//   for (Int_t iph = 0; iph < nCaloPhotons; iph++){
//     AliCaloPhoton *ph = (AliCaloPhoton*)fCaloPhotonsPHOS->At(iph);
//     if (!ph->IsTOFOK()) continue;

//     weight = 1.;
//     if(fIsMC){
//       primary = ph->GetPrimary();
//       weight = ph->GetWeight();
//     }
//   }
// } // end of FillCaloPhotons()
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::FillMgg()
{

  Int_t primary1 = -1;
  Int_t primary2 = -1;
  Double_t TruePi0Pt = 0;
  Double_t TrueK0SPt = 0;
  Double_t TrueL0Pt = 0;
  Int_t commonID = -1;

  Int_t nCaloPhotons = fCaloPhotonsPHOS->GetEntries();

  for (Int_t i1 = 0; i1 < nCaloPhotons - 1; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1);
    if (!ph1->IsTOFOK())
      continue;

    for (Int_t i2 = i1 + 1; i2 < nCaloPhotons; i2++) {
      AliCaloPhoton* ph2 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i2);
      if (!ph2->IsTOFOK())
        continue;

      Double_t en1 = ph1->E();
      Double_t en2 = ph2->E();

      Double_t pT1 = ph1->Pt();
      Double_t pT2 = ph2->Pt();
      Double_t pTlead = ph1->Pt();
      if (pT2 > pTlead)
        pTlead = pT2;

      TLorentzVector p12 = *ph1 + *ph2;

      Double_t Weight = 1.;

      if (fIsMC){

        commonID = FindCommonParent(primary1,primary2);

        if(commonID > -1) Weight = ph1->GetWeight();
        else Weight = ph1->GetWeight()*ph2->GetWeight();

      }
      // const Char_t* CutNames[kPIDCuts] = { "All", "DispCut", "CPVCut", "PhIDCut" };
      const Bool_t PhotonCutFlag[kPIDCuts] = { kTRUE, // All
                                               ph1->IsDispOK() && ph2->IsDispOK(), // DispCut
                                               ph1->IsCPVOK() && ph2->IsCPVOK(), // CPVCut
                                               ph1->IsPhoton()  && ph2->IsPhoton() }; // PhIDCut

      for (Int_t icut = 0; icut < kPIDCuts; icut++) {
        if (!PhotonCutFlag[icut])
          continue;
        fHistMgg[icut]->Fill(p12.M(), p12.Pt());
        fHistMggW[icut]->Fill(p12.M(), p12.Pt(), Weight);
        if (fDoClustQA) {
          if (icut == 0 && ph1->Module() == ph2->Module())
            fHistMggMod[ph1->Module() - 1]->Fill(p12.M(), p12.Pt(), Weight);
          if (icut == 3 && ph1->Module() == ph2->Module())
            fHistMggPhIDCutMod[ph1->Module() - 1]->Fill(p12.M(), p12.Pt(), Weight);
        }
      }

      fHistMgg2[0]->Fill(p12.M(), p12.Pt());
      if (ph1->IsDispOK() && ph2->IsDispOK()) fHistMgg2[1]->Fill(p12.M(), p12.Pt(), Weight);
      if (ph1->IsCPVOK() && ph2->IsCPVOK()) fHistMgg2[2]->Fill(p12.M(), p12.Pt(), Weight);
      if (ph1->IsPhoton() && ph2->IsPhoton()) fHistMgg2[3]->Fill(p12.M(), p12.Pt(), Weight);

      if (fIsMC && PhotonCutFlag[kPIDCuts-1]) {
        if (ph1->GetPrimary() == ph2->GetPrimary() && ph1->GetPrimary() > 0) {
          AliAODMCParticle* mother = (AliAODMCParticle*)fMCStack->At(ph1->GetPrimary());
          if (mother->GetPdgCode() == 111) {
            fHistTruePi0MggVsRecPt->Fill(p12.M(), p12.Pt(), Weight);
            fHistTruePi0MggVsTruePt->Fill(p12.M(), mother->Pt(), Weight);
            if (ph1->GetPrimaryAtVertex() == ph2->GetPrimaryAtVertex() && (ph1->GetPrimaryAtVertex() > -1)){
              AliAODMCParticle *PrimMother = (AliAODMCParticle*)fMCStack->At(ph1->GetPrimaryAtVertex());
              fHistTruePi0MotherIDvsRecPt->Fill(p12.Pt(), GetPi0SourceIndex(PrimMother->GetPdgCode()));
              fHistTruePi0MotherIDvsRecPtW->Fill(p12.Pt(), GetPi0SourceIndex(PrimMother->GetPdgCode()), Weight);
              if (PrimMother->GetPdgCode() == 310){
                fHistTruePi0MggVsRecPtFromK0s->Fill(p12.M(), p12.Pt(), Weight);
              } else if (PrimMother->GetPdgCode() == 130){
                  fHistTruePi0MggVsRecPtFromK0L->Fill(p12.M(), p12.Pt(), Weight);
              } else if (PrimMother->GetPdgCode() == 3122){
                fHistTruePi0MggVsRecPtFromLambda->Fill(p12.M(), p12.Pt(), Weight);
              } else if (TMath::Abs(PrimMother->GetPdgCode()) == 2112 || TMath::Abs(PrimMother->GetPdgCode()) == 2212){
                fHistTruePi0MggVsRecPtFromPandN->Fill(p12.M(), p12.Pt(), Weight);
              } else if (PrimMother->GetPdgCode() == 321){
                fHistTruePi0MggVsRecPtFromPions->Fill(p12.M(), p12.Pt(), Weight);
              } else if (PrimMother->GetPdgCode() == 211){
                fHistTruePi0MggVsRecPtFromKaons->Fill(p12.M(), p12.Pt(), Weight);
              }
            }
          } else if (mother->GetPdgCode() == 221){
            fHistTrueEtaMggVsRecPt->Fill(p12.M(), p12.Pt(), Weight);
            fHistTrueEtaMggVsTruePt->Fill(p12.M(), mother->Pt(), Weight);
          }
        }
      }

    } // ph2 loop for real event
  }   // ph1 loop for real event

  for (Int_t i1 = 0; i1 < nCaloPhotons; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1);
    if (!ph1->IsTOFOK())
      continue;

    for (Int_t ev = 0; ev < fCaloPhotonsPHOSList->GetSize(); ev++) {
      TClonesArray* mixPHOS = static_cast<TClonesArray*>(fCaloPhotonsPHOSList->At(ev));
      for (Int_t i2 = 0; i2 < mixPHOS->GetEntriesFast(); i2++) {
        AliCaloPhoton* ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if (!ph2->IsTOFOK())
          continue;

        Double_t Weight = ph1->GetWeight()*ph2->GetWeight();

        Double_t en1 = ph1->E();
        Double_t en2 = ph2->E();

        Double_t pT1 = ph1->Pt();
        Double_t pT2 = ph2->Pt();
        Double_t pTlead = ph1->Pt();
        if (pT2 > pTlead)
          pTlead = pT2;

        TLorentzVector p12 = *ph1 + *ph2;
        // const Char_t* CutNames[kPIDCuts] = { "All", "DispCut", "CPVCut", "PhIDCut" };
        const Bool_t PhotonCutFlag[kPIDCuts] = { kTRUE, // All
                                                 ph1->IsDispOK() && ph2->IsDispOK(), // DispCut
                                                 ph1->IsCPVOK() && ph2->IsCPVOK(), // CPVCut
                                                 ph1->IsPhoton() && ph2->IsPhoton() }; // PhIDCut

        for (Int_t icut = 0; icut < kPIDCuts; icut++) {
          if (!PhotonCutFlag[icut])
            continue;
          fHistMixMgg[icut]->Fill(p12.M(), p12.Pt());
          fHistMixMggW[icut]->Fill(p12.M(), p12.Pt(), Weight);
          if (fDoClustQA) {
            if (icut == 0 && ph1->Module() == ph2->Module())
              fHistMixMggMod[ph1->Module() - 1]->Fill(p12.M(), p12.Pt(), Weight);
            if (icut == 3 && ph1->Module() == ph2->Module())
              fHistMixMggPhIDCutMod[ph1->Module() - 1]->Fill(p12.M(), p12.Pt(), Weight);
          }
        }

        fHistMixMgg2[0]->Fill(p12.M(), p12.Pt());
        if (ph1->IsDispOK() && ph2->IsDispOK()) fHistMixMgg2[1]->Fill(p12.M(), p12.Pt(), Weight);
        if (ph1->IsCPVOK() && ph2->IsCPVOK()) fHistMixMgg2[2]->Fill(p12.M(), p12.Pt(), Weight);
        if (ph1->IsPhoton() && ph2->IsPhoton()) fHistMixMgg2[3]->Fill(p12.M(), p12.Pt(), Weight);

      } // ph2 loop for mixed events
    }   // loop over fNMixEvents
  }     // ph1 loop for mixed events

} // end of FillMgg()
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::EstimatePIDCutEfficiency()
{

  TLorentzVector p12;
  Double_t m12 = 0;
  Double_t pT = 0;

  const Int_t nCaloPhotons = fCaloPhotonsPHOS->GetEntriesFast();
  for (Int_t i1 = 0; i1 < nCaloPhotons; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1);
    if (!ph1->IsTOFOK())
      continue;
    // apply tight cut to photon1
    if (ph1->Energy() < 0.5 || ph1->GetNsigmaCPV() < 4 || ph1->GetNsigmaCoreDisp() > 2.5)
      continue;
    for (Int_t i2 = 0; i2 < nCaloPhotons; i2++) {
      if (i2 == i1)
        continue;

      AliCaloPhoton* ph2 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i2);
      if (!ph2->IsTOFOK())
        continue;

      Double_t Weight = ph1->GetWeight()*ph2->GetWeight();

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      pT = ph2->Pt();

      const Bool_t PhotonCutFlag[kPIDCuts] = { kTRUE,  // All
                                               ph2->IsDispOK(),  // DispCut
                                               ph2->IsCPVOK(), // CPVCut
                                               ph2->IsPhoton() }; // PhIDCut

      for (Int_t icut = 0; icut < kPIDCuts; icut++) {
        if (!PhotonCutFlag[icut])
          continue;
        fHistMggCutEff[icut]->Fill(m12, pT, Weight);
        if (fDoClustQA && ph1->Module() == ph2->Module())
          fHistMggCutEffMod[icut][ph1->Module() - 1]->Fill(m12, pT, Weight);
      }

    } // loop over ph1 in real event
  }   // loop over ph2 in real event

  for (Int_t i1 = 0; i1 < nCaloPhotons; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1);
    if (!ph1->IsTOFOK())
      continue;
      
    // apply tight cut to photon1
    if(ph1->Energy() < 0.5 || ph1->GetNsigmaCPV() < 4 || ph1->GetNsigmaCoreDisp() > 2.5) continue;

    for (Int_t ev = 0; ev < fCaloPhotonsPHOSList->GetSize(); ev++) {
      TClonesArray* mixPHOS = static_cast<TClonesArray*>(fCaloPhotonsPHOSList->At(ev));
      for (Int_t i2 = 0; i2 < mixPHOS->GetEntriesFast(); i2++) {
        AliCaloPhoton* ph2 = (AliCaloPhoton*)mixPHOS->At(i2);
        if (!ph2->IsTOFOK())
          continue;

        Double_t Weight = ph1->GetWeight()*ph2->GetWeight();

        p12 = *ph1 + *ph2;
        m12 = p12.M();
        pT = ph2->Pt();

        const Bool_t PhotonCutFlag[kPIDCuts] = { kTRUE, // All
                                                 ph2->IsDispOK(),  // DispCut
                                                 ph2->IsCPVOK(),  // CPVCut
                                                 ph2->IsPhoton() }; // PhIDCut

        for (Int_t icut = 0; icut < kPIDCuts; icut++) {
          if (!PhotonCutFlag[icut])
            continue;
          fHistMixMggCutEff[icut]->Fill(m12, pT, Weight);
          if (fDoClustQA && ph1->Module() == ph2->Module())
            fHistMixMggCutEffMod[icut][ph1->Module() - 1]->Fill(m12, pT, Weight);
        }

      } // ph2 loop for mixed events
    }   // loop over fNMixEvents
  }     // ph1 loop for mixed events

} // end of EstimatePIDCutEfficiency()
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::EstimateTOFCutEfficiency()
{

  TLorentzVector p12;
  Double_t m12 = 0;
  Double_t energy = 0;

  const Int_t nCaloPhotons = fCaloPhotonsPHOS->GetEntriesFast();
  for (Int_t i1 = 0; i1 < nCaloPhotons; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1);
    if (!ph1->IsTOFOK())
      continue;

    for (Int_t i2 = 0; i2 < nCaloPhotons; i2++) {
      if (i2 == i1)
        continue;
      AliCaloPhoton* ph2 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i2);

      p12 = *ph1 + *ph2;
      m12 = p12.M();
      energy = ph2->Energy();

      if (!(ph1->IsCPVOK() && ph2->IsCPVOK() && ph1->IsDispOK() && ph2->IsDispOK()))
        continue; // PhIDCut

      fHistMggTOFCutEffBase->Fill(m12, energy);
      if (ph2->IsTOFOK())
        fHistMggTOFCutEffProbe->Fill(m12, energy);

    } // loop over ph1 in real event
  }   // loop over ph2 in real event

  //Mixed events loop
  for (Int_t i1 = 0; i1 < nCaloPhotons; i1++) {
    AliCaloPhoton* ph1 = (AliCaloPhoton*)fCaloPhotonsPHOS->At(i1);
    if (!ph1->IsTOFOK())
      continue;

    for (Int_t ev = 0; ev < fCaloPhotonsPHOSList->GetSize(); ev++) {
      TClonesArray* mixPHOS = static_cast<TClonesArray*>(fCaloPhotonsPHOSList->At(ev));
      for (Int_t i2 = 0; i2 < mixPHOS->GetEntriesFast(); i2++) {
        AliCaloPhoton* ph2 = (AliCaloPhoton*)mixPHOS->At(i2);

        p12 = *ph1 + *ph2;
        m12 = p12.M();
        energy = ph2->Energy();

        if (!(ph1->IsCPVOK() && ph2->IsCPVOK() && ph1->IsDispOK() && ph2->IsDispOK()))
          continue; // PhIDCut

        fHistMixMggTOFCutEffBase->Fill(m12, energy);
        if (ph2->IsTOFOK()){
          fHistMixMggTOFCutEffProbe->Fill(m12, energy);
        }

      } // ph2 loop for mixed events
    }   // loop over fNMixEvents
  }     // ph1 loop for mixed events
} // end of EstimateTOFCutEfficiency()
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::ProcessMCParticles()
{

  if (!fIsMC)
    return;

  Long_t multMC = fMCStack->GetEntriesFast();
  for (Long_t iMC = 0; iMC < multMC; iMC++) {

    AliAODMCParticle* particle = (AliAODMCParticle*)fMCStack->At(iMC);

    Double_t pt = particle->Pt();
    Double_t phi = particle->Phi();

    // PHOS phi region
    const Double_t phiMax = 320 * TMath::Pi() / 180;
    const Double_t phiMin = 250 * TMath::Pi() / 180;

    while (phi < 0)
      phi += TMath::TwoPi();
    while (phi > TMath::TwoPi())
      phi -= TMath::TwoPi();

    Int_t pdg = particle->GetPdgCode();
    Int_t mpdg;
    Int_t iPart = 0;

    Double_t MCWeight = 1.;

    TF1 *f1Pi0Weight   = (TF1*)GetAdditionalPi0PtWeightFunction(fCentrality);
    TF1 *f1EtaWeight   = (TF1*)GetAdditionalEtaPtWeightFunction(fCentrality);
    TF1 *f1GammaWeight = (TF1*)GetAdditionalGammaPtWeightFunction(fCentrality);
    TF1 *f1K0SWeight   = (TF1*)GetAdditionalK0SPtWeightFunction(fCentrality);

    Double_t PrimPt = -1.;

    if (pdg == 111){
      MCWeight = f1Pi0Weight->Eval(pt);
      if (MCMotherIsFound(iMC, PrimPt, 321)) MCWeight = f1Pi0Weight->Eval(PrimPt)*f1K0SWeight->Eval(PrimPt); // K+-
      else if (MCMotherIsFound(iMC, PrimPt, 221)) MCWeight = f1Pi0Weight->Eval(PrimPt)*f1EtaWeight->Eval(PrimPt); // eta
      else if (MCMotherIsFound(iMC, PrimPt, 310)) MCWeight = f1Pi0Weight->Eval(PrimPt)*f1K0SWeight->Eval(PrimPt); //K0s
      else if (MCMotherIsFound(iMC, PrimPt, 130)) MCWeight = f1Pi0Weight->Eval(PrimPt)*f1K0SWeight->Eval(PrimPt); //K0L
    }
    else if (TMath::Abs(pdg) == 211)
      MCWeight = f1Pi0Weight->Eval(pt);
    else if (pdg == 22){
      if (MCMotherIsFound(iMC, PrimPt, 111)) MCWeight = f1Pi0Weight->Eval(PrimPt); // pi0
      else if (MCMotherIsFound(iMC, PrimPt, 221)) MCWeight = f1Pi0Weight->Eval(PrimPt)*f1EtaWeight->Eval(PrimPt); // eta
      else if (MCMotherIsFound(iMC, PrimPt, 321)) MCWeight = f1Pi0Weight->Eval(PrimPt)*f1K0SWeight->Eval(PrimPt); // K+-
      else if (MCMotherIsFound(iMC, PrimPt, 310)) MCWeight = f1Pi0Weight->Eval(PrimPt)*f1K0SWeight->Eval(PrimPt); //K0s
      else if (MCMotherIsFound(iMC, PrimPt, 130)) MCWeight = f1Pi0Weight->Eval(PrimPt)*f1K0SWeight->Eval(PrimPt); //K0L
    }
    else if (pdg == 221) MCWeight = f1Pi0Weight->Eval(pt)*f1EtaWeight->Eval(pt); // eta
    else if (pdg == 321) MCWeight = f1Pi0Weight->Eval(pt)*f1K0SWeight->Eval(pt); // K+-
    else if (pdg == 310) MCWeight = f1Pi0Weight->Eval(pt)*f1K0SWeight->Eval(pt); // K0s
    else if (pdg == 130) MCWeight = f1Pi0Weight->Eval(pt)*f1K0SWeight->Eval(pt); // K0L

    if (pdg == 111){
      Int_t primlb = FindMCPrimaryMotherLabel(iMC);
      // Int_t primlb = particle->GetMother();
      AliAODMCParticle *PrimMother = nullptr;
      if (primlb > -1){
        AliAODMCParticle *PrimMother = (AliAODMCParticle*)fMCStack->At(primlb);
        fHistTruePi0MotherIDvsTruePt->Fill(pt, GetPi0SourceIndex(PrimMother->GetPdgCode()), MCWeight);
      }
    }

    if (fClearGammaWeight && pdg == 22)
      MCWeight = f1GammaWeight->Eval(pt);

    if (!IsPrimaryMCParticle(particle))
      continue; // particles in R > 1.0 from prim vertex are concidered
    if (TMath::Abs(particle->Y()) > 0.5)
      continue; // only partciles in |y| < 0.5 window are concidered

    Int_t ParticleIndex = 15;

    if (TMath::Abs(pdg) == 310)
      ParticleIndex = 1; // K0s
    else if (TMath::Abs(pdg) == 3122)
      ParticleIndex = 2; // Lambda
    else if (TMath::Abs(pdg) == 130)
      ParticleIndex = 3; // K0L
    else if (TMath::Abs(pdg) == 2212)
      ParticleIndex = 4; // proton
    else if (TMath::Abs(pdg) == 2112)
      ParticleIndex = 5; // neutron
    else if (TMath::Abs(pdg) == 211)
      ParticleIndex = 6; // pion
    else if (TMath::Abs(pdg) == 321)
      ParticleIndex = 7; // kaon
    else if (TMath::Abs(pdg) == 113 || TMath::Abs(pdg) == 213)
      ParticleIndex = 8; // rho 0,+,-
    else if (TMath::Abs(pdg) == 3222 || TMath::Abs(pdg) == 3212 || TMath::Abs(pdg) == 3112)
      ParticleIndex = 9; // Sigma
    else if (TMath::Abs(pdg) == 2224 || TMath::Abs(pdg) == 2214 || TMath::Abs(pdg) == 2114 || TMath::Abs(pdg) == 1114)
      ParticleIndex = 10; // Delta
    else if (TMath::Abs(pdg) == 313 || TMath::Abs(pdg) == 323)
      ParticleIndex = 11; // K*
    else
      ParticleIndex = 15;

    fHistMCPartIDvsPt->Fill(pt, ParticleIndex+0.5);
    fHistMCPartIDvsPtW->Fill(pt, ParticleIndex+0.5, MCWeight);

    if (pdg == 111) { // pi0
      fHistPrimPi0Pt->Fill(pt);
      fHistPrimPi0PtW->Fill(pt, MCWeight);
      if (InPHOSAcc(particle)){
        fHistPrimPi0InAccPt->Fill(pt);
        fHistPrimPi0InAccPtW->Fill(pt, MCWeight);
      }

      if (particle->GetDaughterLabel(0) > -1 && particle->GetDaughterLabel(1) > -1){

        AliAODMCParticle* ph1 = (AliAODMCParticle*)fMCStack->At(particle->GetDaughterLabel(0));
        AliAODMCParticle* ph2 = (AliAODMCParticle*)fMCStack->At(particle->GetDaughterLabel(1));

        if (ph1->GetPdgCode() == 22 && ph2->GetPdgCode() == 22) {
          if (InPHOSAcc(ph1) && InPHOSAcc(ph2)){
            fHistPrimPi0BothPhInAccPt->Fill(pt);
            fHistPrimPi0BothPhInAccPtW->Fill(pt, MCWeight);
          }
          if (InPHOSAcc(ph1) || InPHOSAcc(ph2)) {
            Double_t ptPh = (InPHOSAcc(ph1)) ? ph1->Pt() : ph2->Pt();
            fHistPrimPi0OnePhInAccPt->Fill(pt);
            fHistPrimPi0PtvsPhInAccPt->Fill(pt, ptPh / pt);
          }
        }
      }
    }

    else if (pdg == 221) {
      fHistPrimEtaPt->Fill(pt);
      fHistPrimEtaPtW->Fill(pt, MCWeight);
      if (InPHOSAcc(particle)){
        fHistPrimEtaInAccPt->Fill(pt);
      }

      if (particle->GetDaughterLabel(0) > -1 && particle->GetDaughterLabel(1) > -1){

        AliAODMCParticle* ph1 = (AliAODMCParticle*)fMCStack->At(particle->GetDaughterLabel(0));
        AliAODMCParticle* ph2 = (AliAODMCParticle*)fMCStack->At(particle->GetDaughterLabel(1));

        if (ph1->GetPdgCode() == 22 && ph2->GetPdgCode() == 22) {
          if (InPHOSAcc(ph1) && InPHOSAcc(ph2))
            fHistPrimEtaBothPhInAccPt->Fill(pt);
          if (InPHOSAcc(ph1) || InPHOSAcc(ph2)) {
            Double_t ptPh = (InPHOSAcc(ph1)) ? ph1->Pt() : ph2->Pt();
            fHistPrimEtaOnePhInAccPt->Fill(pt);
            fHistPrimEtaPtvsPhInAccPt->Fill(pt, ptPh / pt);
          }
        }
      }
    }

    else if (pdg == 22) {
      fHistPrimGammaPt->Fill(pt);
      fHistPrimGammaPtW->Fill(pt, MCWeight);
      if (InPHOSAcc(particle)){
        fHistPrimGammaInAccPt->Fill(pt);
        fHistPrimGammaInAccPtW->Fill(pt, MCWeight);
      }
    }
  }
}
//_____________________________________________________________________________
Bool_t AliAnalysisPHOSNeutralMesonsAndPhotons::PassSTDCut(AliAODCluster *cluster)
{
  if(cluster->GetM20() > 2.0) return kFALSE;
  if(cluster->E() > 1.0 && cluster->GetM02() < 0.1) return kFALSE;
  if(cluster->E() > 2.0 && cluster->GetM20() < 0.1) return kFALSE;
  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisPHOSNeutralMesonsAndPhotons::InPHOSAcc(AliAODMCParticle* particle)
{

  // PHOS phi region
  const Double_t phiMax = 320 * TMath::Pi() / 180;
  const Double_t phiMin = 250 * TMath::Pi() / 180;
  const Double_t yMax = 0.125; // PHOS rapidity window

  Double_t phi = particle->Phi();

  while (phi < 0)
    phi += TMath::TwoPi();
  while (phi > TMath::TwoPi())
    phi -= TMath::TwoPi();

  if (phi <= phiMax && phi >= phiMin && TMath::Abs(particle->Y()) <= yMax)
    return kTRUE;
  else
    return kFALSE;
}
//_____________________________________________________________________________
Int_t AliAnalysisPHOSNeutralMesonsAndPhotons::FindClusterPrimary(AliCaloPhoton *ph,  Bool_t&sure)
{
  //Finds primary and estimates if it unique one?
  //First check can it be photon/electron
  AliVCluster *clu = (AliVCluster*)ph->GetCluster();
  const Double_t emFraction=0.9; //part of energy of cluster to be assigned to EM particle
  Int_t n = clu->GetNLabels();

  for(Int_t i=0;  i<n;  i++){
    Int_t label = clu->GetLabelAt(i);
    AliAODMCParticle *p =  (AliAODMCParticle*)fMCStack->At(label);
    Int_t pdg = p->PdgCode() ;
    if(pdg==22  ||  pdg==11 || pdg == -11){
      if(p->E()>emFraction*clu->E()){
        sure=kTRUE ;
        return label;
      }
    }
  }

  Double_t *Ekin = new Double_t[n];

  for(Int_t i=0;  i<n;  i++){
    Int_t label = clu->GetLabelAt(i);
    AliAODMCParticle* p = (AliAODMCParticle*)fMCStack->At(label);
    Ekin[i]=p->P() ;  // estimate of kinetic energy
    if(p->PdgCode()==-2212  ||  p->PdgCode()==-2112){
      Ekin[i]+=1.8  ;  //due to annihilation
    }
  }

  Int_t iMax=0;
  Double_t eMax=0.,eSubMax=0. ;
  for(Int_t i=0;  i<n;  i++){
    if(Ekin[i]>eMax){
      eSubMax=eMax;
      eMax=Ekin[i];
      iMax=i;
    }
  }
  if(eSubMax>0.8*eMax)//not obvious primary
    sure=kFALSE;
  else
    sure=kTRUE;

  delete[]  Ekin;

  return clu->GetLabelAt(iMax);

}
//_____________________________________________________________________________
Bool_t AliAnalysisPHOSNeutralMesonsAndPhotons::IsFrom(Int_t label, Double_t &TruePt, const Int_t target_pdg)
{

  Int_t motherid = -1;
  Int_t pdg=0;
  Double_t pT = 0;

    AliAODMCParticle *p1 = (AliAODMCParticle*)fMCStack->At(label);
    AliAODMCParticle *mp = 0x0;
    motherid = p1->GetMother();

    while(motherid > -1){
      mp = (AliAODMCParticle*)fMCStack->At(motherid);
      pT = mp->Pt();
      pdg = mp->PdgCode(); 

      //if(TMath::Abs(pdg) == target_pdg && Rho(mp) < 1.0){//pi0 from primary vertex
      if(TMath::Abs(pdg) == target_pdg && IsPrimaryMCParticle(mp)){ //pi0 from primary vertex
        TruePt = pT;
        return kTRUE;
      }
      motherid = mp->GetMother();
    }

  return kFALSE;

} // end of IsFrom()
//_____________________________________________________________________________
AliAODMCParticle* AliAnalysisPHOSNeutralMesonsAndPhotons::FindMCPrimaryMother(AliAODMCParticle* particle)
{

  Int_t primlb = particle->GetMother();
  if (primlb == -1)
    return particle;

  AliAODMCParticle* prim = (AliAODMCParticle*)fMCStack->At(primlb);
  AliAODMCParticle* mprim = nullptr;

  while (!IsPrimaryMCParticle(prim) && (primlb > -1)) {
    if (prim->GetMother() < 0)
      break;
    primlb = prim->GetMother();
    mprim = (AliAODMCParticle*)fMCStack->At(primlb);
    if (TMath::Abs(mprim->GetPdgCode()) <= 6 || TMath::Abs(mprim->GetPdgCode()) == 21) {
      break;
    }
    prim = mprim;
  }

  if (IsPrimaryMCParticle(prim))
    return prim;
  else
    return nullptr;
}
//_____________________________________________________________________________
AliAODMCParticle* AliAnalysisPHOSNeutralMesonsAndPhotons::FindMCPrimaryMother(Int_t inputlb)
{

  Int_t primlb = inputlb;
  if (primlb < 0)
    return nullptr;

  AliAODMCParticle* prim = (AliAODMCParticle*)fMCStack->At(primlb);
  AliAODMCParticle* mprim = nullptr;

  while (!IsPrimaryMCParticle(prim) && (primlb > -1)) {
    if (prim->GetMother() < 0)
      break;
    primlb = prim->GetMother();
    mprim = (AliAODMCParticle*)fMCStack->At(primlb);
    if (TMath::Abs(mprim->GetPdgCode()) <= 6 || TMath::Abs(mprim->GetPdgCode()) == 21) {
      break;
    }
    prim = mprim;
  }
  if (IsPrimaryMCParticle(prim))
    return prim;
  else
    return nullptr;
}
//_____________________________________________________________________________
Int_t AliAnalysisPHOSNeutralMesonsAndPhotons::FindMCPrimaryMotherLabel(Int_t inputlb)
{

  if (inputlb < 0)
    return -1;

  AliAODMCParticle* prim = (AliAODMCParticle*)fMCStack->At(inputlb);
  Int_t mprimlb = inputlb;

  while (!IsPrimaryMCParticle(prim) && (mprimlb > -1)) {
    if (TMath::Abs(prim->GetPdgCode()) <= 6 || TMath::Abs(prim->GetPdgCode()) == 21) {
      break;
    }
    if (prim->GetMother() < 0) break;
    mprimlb = prim->GetMother();
    prim = (AliAODMCParticle*)fMCStack->At(mprimlb);
  }

  return mprimlb;
}
//_____________________________________________________________________________
Bool_t AliAnalysisPHOSNeutralMesonsAndPhotons::MCMotherIsFound(Int_t inputlb, Double_t &primpt, Int_t primpdg){

  Int_t primlb = inputlb;

  while (primlb > -1) {
    AliAODMCParticle* prim = (AliAODMCParticle*)fMCStack->At(primlb);
    // if (TMath::Abs(prim->GetPdgCode()) <= 6 || TMath::Abs(prim->GetPdgCode()) == 21) {  
    //   return kFALSE;
    // }
    if (TMath::Abs(prim->GetPdgCode()) == primpdg && IsPrimaryMCParticle(prim)){
      primpt = prim->Pt();
      return kTRUE;
    }
    primlb = prim->GetMother();
  }

  return kFALSE;

}
//_____________________________________________________________________________
Bool_t AliAnalysisPHOSNeutralMesonsAndPhotons::MCMotherIsFound(Int_t inputlb, Double_t &primpt, Int_t &primlb, Int_t primpdg){

  Int_t motherlb = inputlb;

  while (motherlb > -1) {
    AliAODMCParticle* prim = (AliAODMCParticle*)fMCStack->At(motherlb);
    // if (TMath::Abs(prim->GetPdgCode()) <= 6 || TMath::Abs(prim->GetPdgCode()) == 21) {
    //   return kFALSE;
    // }
    if (TMath::Abs(prim->GetPdgCode()) == primpdg && IsPrimaryMCParticle(prim)){
      primpt = prim->Pt();
      primlb = motherlb;
      return kTRUE;
    }
    motherlb = prim->GetMother();
  }

  return kFALSE;

}
//_____________________________________________________________________________
Int_t AliAnalysisPHOSNeutralMesonsAndPhotons::FindCommonParent(Int_t iPart, Int_t jPart)
{
  //check if there is a common parent for particles i and j
  // -1: no common parent or wrong iPart/jPart

  Int_t ntrack = fMCStack->GetEntriesFast();
  if(iPart==-1 || iPart>=ntrack || jPart==-1 || jPart>=ntrack) return -1;

  Int_t iprim1 = iPart;

  while(iprim1>-1){
    AliAODMCParticle *p1 = dynamic_cast<AliAODMCParticle*>(fMCStack->At(iprim1));
    Int_t pdg1 = p1->GetPdgCode();
    if(TMath::Abs(pdg1) <= 6 || TMath::Abs(pdg1) == 21){//reject quarks and gluons as parents.
      iprim1 = p1->GetMother();
      continue;
    }
    if(iprim1 == 0 || iprim1 == 1){//colliding protons
      iprim1 = p1->GetMother();
      continue;
    }

    Int_t iprim2 = jPart;
    while(iprim2>-1){
      AliAODMCParticle *p2 = dynamic_cast<AliAODMCParticle*>(fMCStack->At(iprim2));
      Int_t pdg2 = p2->GetPdgCode();

      if(TMath::Abs(pdg2) <= 6 || TMath::Abs(pdg2) == 21){//reject quarks and gluons as parents.
        iprim2 = p2->GetMother();
        continue;
      }
      if(iprim2 == 0 || iprim2 == 1){//colliding protons
        iprim2 = p2->GetMother();
        continue;
      }

      if(iprim1==iprim2) return iprim1;
      iprim2 = p2->GetMother();

    }

    iprim1 = p1->GetMother();
  }

  return -1;
}
//_____________________________________________________________________________
Bool_t AliAnalysisPHOSNeutralMesonsAndPhotons::IsPrimaryMCParticle(AliAODMCParticle* particle)
{

  Bool_t isPrimary = kFALSE;

  const Double_t rcut = 1.;

  Double_t r2 = (particle->Xv() - fVertex[0]) * (particle->Xv() - fVertex[0]) +
                (particle->Yv() - fVertex[1]) * (particle->Yv() - fVertex[1]);

  if (TMath::Sqrt(r2) <= rcut)
    isPrimary = kTRUE;

  return isPrimary;
}
//_____________________________________________________________________________
Bool_t AliAnalysisPHOSNeutralMesonsAndPhotons::IsPrimaryMCParticle3D(AliAODMCParticle* particle)
{

  Bool_t isPrimary = kFALSE;

  const Double_t rcut = 1;
  // const AliVVertex *primVtxMC = fMCEvent->GetPrimaryVertex();
  const AliVVertex* primVtxMC = fEvent->GetPrimaryVertex();
  Double_t mcVtxX = primVtxMC->GetX();
  Double_t mcVtxY = primVtxMC->GetY();
  Double_t mcVtxZ = primVtxMC->GetZ();

  Double_t r2 = (particle->Xv() - mcVtxX) * (particle->Xv() - mcVtxX) +
                (particle->Yv() - mcVtxY) * (particle->Yv() - mcVtxY) +
                (particle->Zv() - mcVtxZ) * (particle->Zv() - mcVtxZ);

  if (TMath::Sqrt(r2) <= rcut)
    isPrimary = kTRUE;

  return isPrimary;
}
//_____________________________________________________________________________
Int_t AliAnalysisPHOSNeutralMesonsAndPhotons::GetPi0SourceIndex(Int_t MotherPdg){

  if (TMath::Abs(MotherPdg) == 310) return 1; // K0s
  else if (TMath::Abs(MotherPdg) == 3122) return 2; // Lambda
  else if (TMath::Abs(MotherPdg) == 130) return 3; // K0L
  else if (TMath::Abs(MotherPdg) == 2212) return 4; // proton
  else if (TMath::Abs(MotherPdg) == 2112) return 5; // neutron
  else if (TMath::Abs(MotherPdg) == 211) return 6; // pion
  else if (TMath::Abs(MotherPdg) == 321) return 7; // kaon
  else if (TMath::Abs(MotherPdg) == 113 || TMath::Abs(MotherPdg) == 213 ) return 8; // rho 0,+,-
  else if (TMath::Abs(MotherPdg) == 3222 || TMath::Abs(MotherPdg) == 3212 || 
            TMath::Abs(MotherPdg) == 3112  ) return 9; // Sigma
  else if (TMath::Abs(MotherPdg) == 2224 || TMath::Abs(MotherPdg) == 2214 || 
            TMath::Abs(MotherPdg) == 2114 || TMath::Abs(MotherPdg) == 1114  ) return 10; // Delta
  else if (TMath::Abs(MotherPdg) == 313 || TMath::Abs(MotherPdg) == 323   ) return 11; // K*
  else if (TMath::Abs(MotherPdg) == 111) return 12; //pi0
  else return 15;

}
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::SelectCentrality()
{

  Float_t fCentralityV0M = -1.;
  Float_t fCentralityCL0 = -1.;
  Float_t fCentralityCL1 = -1.;
  Float_t fCentralityV0A = -1.;
  Float_t fCentralityV0C = -1.;
  Float_t fCentralityZNA = -1.;
  Float_t fCentralityZNC = -1.;

  if (fCentralityEstimator.Contains("V0") || fCentralityEstimator.Contains("ZN") || fCentralityEstimator.Contains("CL")) {

    // Get Centrality
    fMultSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
    if (!fMultSelection) {
      // If you get this warning (and fCentralityV0M 300) please check that the AliMultSelectionTask actually ran (before your task)
      AliWarning("AliMultSelection object not found!");
      return;
    } else {
      fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");
      fCentralityCL0 = fMultSelection->GetMultiplicityPercentile("CL0");
      fCentralityCL1 = fMultSelection->GetMultiplicityPercentile("CL1");
      fCentralityV0A = fMultSelection->GetMultiplicityPercentile("V0A");
      fCentralityV0C = fMultSelection->GetMultiplicityPercentile("V0C");
      fCentralityZNA = fMultSelection->GetMultiplicityPercentile("ZNA");
      fCentralityZNC = fMultSelection->GetMultiplicityPercentile("ZNC");
      fCentrality = fMultSelection->GetMultiplicityPercentile(fCentralityEstimator);

      if (fCentrality < 0) {
        AliInfo("Negative centrality!");
        return;
      }

      fHistCentMain->Fill(fCentrality);
      fHistCentV0MvsCL0->Fill(fCentralityV0M, fCentralityCL0);
      fHistCentV0MvsCL1->Fill(fCentralityV0M, fCentralityCL1);
      fHistCentV0MvsV0C->Fill(fCentralityV0A, fCentralityV0C);
      fHistCentCL0vsCL1->Fill(fCentralityCL0, fCentralityCL1);
      fHistCentZNAvsZNC->Fill(fCentralityZNA, fCentralityZNC);

    }

  } else {
    AliInfo(Form("%s is not supported.", fCentralityEstimator.Data()));
    return;
  }
}
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::InitPHOSGeometry()
{
  // Rotation matrixes are set with Tender

  if (fPHOSGeo)
    return;

  fPHOSGeo = AliPHOSGeometry::GetInstance();

  if (!fPHOSGeo) {             // Geometry not yet constructed with Tender
    if (fRunNumber < 209122) { // Run1
      AliError("Can not get geometry from TENDER, creating PHOS geometry for Run1\n");
      fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP", "");
    } else {
      AliError("Can not get geometry from TENDER, creating PHOS geometry for Run2\n");
      fPHOSGeo = AliPHOSGeometry::GetInstance("Run2", "");
    }
    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root", "PHOSRotationMatrixes");
    TObjArray* matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber, "PHOSRotationMatrixes");
    for (Int_t mod = 0; mod < 5; mod++) {
      if (!matrixes->At(mod))
        continue;
      fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)), mod);
    }
  }

  // Read BadMap for MC simulations
  AliOADBContainer badmapContainer(Form("phosBadMap"));
  badmapContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSBadMaps.root", "phosBadMap");
  TObjArray* maps = (TObjArray*)badmapContainer.GetObject(fRunNumber, "phosBadMap");
  if (!maps) {
    AliError("Can not read Bad map\n");
  } else {
    AliInfo(Form("Setting PHOS bad map with name %s \n", maps->GetName()));
    for (Int_t mod = 0; mod < 5; mod++) {
      if (fPHOSBadMap[mod])
        delete fPHOSBadMap[mod];
      TH2I* h = (TH2I*)maps->At(mod);
      if (h)
        fPHOSBadMap[mod] = new TH2I(*h);
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::FillHistogram(const char* key, Double_t x) const
{
  // FillHistogram
  TH1I* tmpI = dynamic_cast<TH1I*>(fOutputContainer->FindObject(key));
  if (tmpI) {
    tmpI->Fill(x);
    return;
  }
  TH1F* tmpF = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key));
  if (tmpF) {
    tmpF->Fill(x);
    return;
  }
  TH1D* tmpD = dynamic_cast<TH1D*>(fOutputContainer->FindObject(key));
  if (tmpD) {
    tmpD->Fill(x);
    return;
  }
  AliInfo(Form("can not find histogram <%s> ", key));
}
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::FillHistogram(const char* key, Double_t x, Double_t y) const
{
  // FillHistogram
  TObject* tmp = fOutputContainer->FindObject(key);
  if (!tmp) {
    AliInfo(Form("can not find histogram <%s> ", key));
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH1F")) {
    ((TH1F*)tmp)->Fill(x, y);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH1D")) {
    ((TH1D*)tmp)->Fill(x, y);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH2F")) {
    ((TH2F*)tmp)->Fill(x, y);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH2D")) {
    ((TH2D*)tmp)->Fill(x, y);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TProfile")) {
    ((TProfile*)tmp)->Fill(x, y);
    return;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s", key, tmp->IsA()->GetName()));
}
//_____________________________________________________________________________
void AliAnalysisPHOSNeutralMesonsAndPhotons::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z) const
{
  // Fills 1D histograms with key
  TObject* tmp = fOutputContainer->FindObject(key);
  if (!tmp) {
    AliInfo(Form("can not find histogram <%s> ", key));
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH2F")) {
    ((TH2F*)tmp)->Fill(x, y, z);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH2D")) {
    ((TH2D*)tmp)->Fill(x, y, z);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH3F")) {
    ((TH3F*)tmp)->Fill(x, y, z);
    return;
  }
  if (tmp->IsA() == TClass::GetClass("TH3D")) {
    ((TH3D*)tmp)->Fill(x, y, z);
    return;
  }
}