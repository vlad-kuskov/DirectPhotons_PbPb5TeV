#include "AliAnalysisPHOSNeutralMesonsAndPhotons.h"

AliAnalysisPHOSNeutralMesonsAndPhotons* AddTaskPHOSNeutralMesonsAndPhotons(
  const char* name = "PHOSNeutralMesonsAndPhotons",
  const TString CollisionSystem = "PbPb",
  const Bool_t isMC = kFALSE,
  const UInt_t trigger = AliVEvent::kINT7,
  const AliCaloTriggerMimicHelper::phosTriggerType PHOSTrigType = AliCaloTriggerMimicHelper::kPHOSAny,
  const Double_t DispCut = 2.5,
  const Double_t CPVCut = 2.5,
  const Double_t TOFCut = 30.,
  const Double_t Emin = 0.2,
  const Bool_t useCorrE = kTRUE,
  const Bool_t usePHOSClusterCuts = kFALSE,
  const Bool_t useMinBiasLowEClustCut = kTRUE,
  const Bool_t doNonlinCorr = kTRUE,
  const TString period = "LHC18q",
  const Int_t NMix = 10,
  const TString OptionSet = "ClustersQA_CellsQA_TracksQA",
  const Bool_t doCentralityStudy = kTRUE,
  const TString CentralityEstimator = "V0M",
  const Float_t CentralityMin = 0.,
  const Float_t CentralityMax = 10.,
  const char* subname = "")
{

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskPHOSNeutralMesonsAndPhotons", "No analysis manager to connect to");
    return NULL;
  }

  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskPHOSNeutralMesonsAndPhotons", "This task requires an input event handler");
    return NULL;
  }

  TString taskname = TString(name);
  taskname += "_" + CollisionSystem;

  if (isMC) {
    taskname += "_MC";
  }
  TString TrigName = "";
  Bool_t MBTrigFlag = kFALSE;
  Bool_t PHOSTrigFlag = kFALSE;
  if (trigger == (UInt_t)AliVEvent::kAny)
    TrigName = "kAny";
  else if (trigger == (UInt_t)AliVEvent::kINT7)
    TrigName = "kINT7";
  else if (trigger == (UInt_t)AliVEvent::kCentral)
    TrigName = "kCentral";
  else if (trigger == (UInt_t)AliVEvent::kSemiCentral)
    TrigName = "kSemiCentral";
  else if (trigger == (UInt_t)AliVEvent::kMB)
    TrigName = "kMB";
  else if (trigger == (UInt_t)AliVEvent::kPHI7) {
    switch (PHOSTrigType) {
      case AliCaloTriggerMimicHelper::kPHOSAny:
        TrigName = "kPHI7Any";
        break;
      case AliCaloTriggerMimicHelper::kPHOSL0:
        TrigName = "kPHI7L0";
        break;
      case AliCaloTriggerMimicHelper::kPHOSL1low:
        TrigName = "kPHI7L1L";
        break;
      case AliCaloTriggerMimicHelper::kPHOSL1med:
        TrigName = "kPHI7L1M";
        break;
      case AliCaloTriggerMimicHelper::kPHOSL1high:
        TrigName = "kPHI7L1H";
        break;
      default:
        break;
    }
  } else
    TrigName = "UnknownTrigger";

  taskname += "_" + TrigName;

  UInt_t MBTriggerMask = AliVEvent::kINT7 | AliVEvent::kCentral | AliVEvent::kSemiCentral | AliVEvent::kMB;
  if (trigger & MBTriggerMask)
    MBTrigFlag = kTRUE;
  else if (trigger & AliVEvent::kPHI7)
    PHOSTrigFlag = kTRUE;

  Int_t systemID = -1;
  if (CollisionSystem == "pp")
    systemID = 0;
  else if (CollisionSystem == "PbPb")
    systemID = 1;
  else if (CollisionSystem == "pPb" || CollisionSystem == "Pbp")
    systemID = 2;

  Bool_t doClusterQA = kFALSE;
  Bool_t doCellQA = kFALSE;

  if (OptionSet.Contains("ClustersQA")) {
    doClusterQA = kTRUE;
  }
  if (OptionSet.Contains("CellsQA")) {
    doCellQA = kTRUE;
  }

  taskname += OptionSet;

  if (doCentralityStudy) {
    taskname += "_" + CentralityEstimator;
    taskname += Form("_%d_%d", (Int_t)CentralityMin, (Int_t)CentralityMax);
  }

  AliAnalysisPHOSNeutralMesonsAndPhotons* task = new AliAnalysisPHOSNeutralMesonsAndPhotons(taskname.Data());

  task->SelectCollisionCandidates(trigger);
  task->SetMBTrigger(MBTrigFlag);
  task->SetPHOSTrigger(PHOSTrigFlag, PHOSTrigType);
  task->SetMC(isMC);
  task->SetClusterCuts(useCorrE, Emin, DispCut, CPVCut, TOFCut, useMinBiasLowEClustCut, usePHOSClusterCuts);
  task->SetQAStudy(doClusterQA, doCellQA);
  task->SetCentralityStudy(doCentralityStudy, CentralityEstimator, CentralityMin, CentralityMax);
  task->SetNEventsForMixing(NMix);

  if (doNonlinCorr) {
    if (isMC){
      TF1* f1NonlinFunc = new TF1("f1NonlinFunc", "[2]*(1.+[0]/(1.+TMath::Power(x/[1],2)))", 0, 200);

      if (period.Contains("LHC18q") || period.Contains("LHC18r")){
        f1NonlinFunc->FixParameter(0, -0.06); // for core E at ZS 20 MeV with only MIP cut
        f1NonlinFunc->FixParameter(1, 0.7);   // for core E at ZS 20 MeV with only MIP cut
        f1NonlinFunc->FixParameter(2, 1.013); // for core E at ZS 20 MeV with only MIP cut
      }
      else{
        f1NonlinFunc->FixParameter(0, 0.);
        f1NonlinFunc->FixParameter(1, 1.);
        f1NonlinFunc->FixParameter(2, 1.);
      }

      task->SetUserNonlinearityFunction(f1NonlinFunc);
    }
    else{
      if (period.Contains("LHC18q") || period.Contains("LHC18r")){
        Double_t PHOSModWeights[4] = {1.025, 1., 1., 1.02}; //Pi0 rec. mass in M1 and M4 is low
        task->SetPHOSModEnergyCorr(PHOSModWeights);
      }
    }
  }

  if (isMC){

    if (OptionSet.Contains("Pi0Weight")){

      // const Int_t Ncen_Pi0 = 11;
      // const Double_t centrality_Pi0[Ncen_Pi0] = {0,5,10,20,30,40,50,60,70,80,100};
      // TArrayD *centarray_Pi0 = new TArrayD(Ncen_Pi0, centrality_Pi0);

      // TObjArray *farray_Pi0 = new TObjArray(Ncen_Pi0-1);
      // TF1 *f1weightPi0[Ncen_Pi0-1];
      // Double_t Pi0Amp = 2.;
      // const Double_t p0_Pi0[Ncen_Pi0-1] = {Pi0Amp * 8.52796e-02, Pi0Amp * 9.57970e-02, Pi0Amp * 1.09042e-01, Pi0Amp * 1.28762e-01, Pi0Amp * 1.51087e-01,
      //                                      Pi0Amp * 1.82705e-01, Pi0Amp * 2.16360e-01, Pi0Amp * 2.37666e-01, Pi0Amp * 2.52706e-01, Pi0Amp * 3.34001e-01};

      // const Double_t p1_Pi0[Ncen_Pi0-1] = {Pi0Amp * 8.34243e-01, Pi0Amp * 8.11715e-01, Pi0Amp * 7.73274e-01, Pi0Amp * 7.28962e-01, Pi0Amp * 6.77506e-01,
      //                                      Pi0Amp * 6.06502e-01, Pi0Amp * 5.31093e-01, Pi0Amp * 4.52193e-01, Pi0Amp * 3.86976e-01, Pi0Amp * 3.22488e-01};

      // const Double_t p2_Pi0[Ncen_Pi0-1] = {9.27577e-01, 9.53380e-01, 9.52280e-01, 9.78872e-01, 9.82192e-01, 1.01124e+00, 1.08236e+00,
      //                                      1.14572e+00, 1.12243e+00, 2.16920e+00};
      // const Double_t p3_Pi0[Ncen_Pi0-1] = {2.13453e-01, 2.09818e-01, 2.03573e-01, 2.00238e-01, 1.94211e-01, 1.87993e-01, 1.94509e-01,
      //                                      1.95069e-01 ,1.75698e-01 ,3.62140e-01};

      // for(Int_t icen = 0; icen < Ncen_Pi0-1; icen++){
      //   f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen), "[0] + [1]*pow(x,[2])*exp(-[3]*x*x)", 0, 100);
      //   f1weightPi0[icen]->SetParameters(p0_Pi0[icen], p1_Pi0[icen], p2_Pi0[icen], p3_Pi0[icen]);
      //   farray_Pi0->Add(f1weightPi0[icen]);
      // }

      //for pi0 pT weighting
      const Int_t Ncen_Pi0 = 11;
      const Double_t centrality_Pi0[Ncen_Pi0] = {0,5,10,20,30,40,50,60,70,80,100};
      TArrayD *centarray_Pi0 = new TArrayD(Ncen_Pi0,centrality_Pi0);

      TObjArray *farray_Pi0 = new TObjArray(Ncen_Pi0-1);
      TF1 *f1weightPi0[Ncen_Pi0-1];

      const Double_t p0[Ncen_Pi0-1] = {2.39991e+02, 1.78111e+02, 1.21109e+02, 6.91786e+01, 4.07880e+01, 2.05577e+01, 1.08079e+01, 4.98463e+00, 2.23119e+00, 1.16590e+00};//Ae
      const Double_t p1[Ncen_Pi0-1] = {3.84238e-01, 3.90424e-01, 3.96450e-01, 4.05978e-01, 4.11701e-01, 4.22538e-01, 4.23961e-01, 4.31410e-01, 4.28510e-01, 3.85667e-01};//Te
      const Double_t p2[Ncen_Pi0-1] = {1.16561e+03, 9.25084e+02, 7.15782e+02, 5.01254e+02, 3.55730e+02, 2.36175e+02, 1.49593e+02, 8.91887e+01, 4.56253e+01, 1.97146e+01};//A
      const Double_t p3[Ncen_Pi0-1] = {3.06202e-01, 3.12875e-01, 3.13599e-01, 3.14247e-01, 3.05660e-01, 3.01051e-01, 2.89748e-01, 2.75489e-01, 2.64481e-01, 2.50078e-01};//T
      const Double_t p4[Ncen_Pi0-1] = {2.73068e+00, 2.72595e+00, 2.70812e+00, 2.68939e+00, 2.64558e+00, 2.61370e+00, 2.56506e+00, 2.53088e+00, 2.51428e+00, 2.43636e+00};//n

      //printf("reading...alien:///alice/cern.ch/user/d/dsekihat/InputPtSpectra/InputPtSpectra_Embedding_PbPb_5.02TeV.root\n");
      //TFile *rootfile_pi0 = TFile::Open("alien:///alice/cern.ch/user/d/dsekihat/InputPtSpectra/InputPtSpectra_Embedding_PbPb_5.02TeV.root","READ");

      for(Int_t icen=0;icen<Ncen_Pi0-1;icen++){
        f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen),"[0] * TMath::Exp(-(TMath::Sqrt(x*x + 0.139*0.139) - 0.139) / [1]) + [2] * TMath::Power(1 + (x*x)/([3]*[3]*[4]) , -[4])",0,100);//TCM fit to PbPb
        f1weightPi0[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        f1weightPi0[icen]->SetNpx(1000);
        farray_Pi0->Add(f1weightPi0[icen]);
      }

      task->SetAdditionalPi0PtWeightFunction(centarray_Pi0, farray_Pi0);

    }

    if (OptionSet.Contains("EtaWeight")){

      // const Int_t Ncen_Eta = 2;
      // const Double_t centrality_Eta[Ncen_Eta] = {0,9999};
      // TArrayD *centarray_Eta = new TArrayD(Ncen_Eta, centrality_Eta);

      // TObjArray *farray_Eta = new TObjArray(Ncen_Eta-1);
      // TF1 *f1weightEta[Ncen_Eta-1];
      // const Double_t p0_eta[Ncen_Eta-1] = { 1.21};
      // const Double_t p1_eta[Ncen_Eta-1] = {-0.499};
      // const Double_t p2_eta[Ncen_Eta-1] = { 1.44};

      // for(Int_t icen = 0; icen < Ncen_Eta-1; icen++){
      //   f1weightEta[icen] = new TF1(Form("f1weightEta_%d", icen), "[0]/(exp(-(x-[1])/[2]) + 1)", 0, 100);
      //   f1weightEta[icen]->SetParameters(p0_eta[icen], p1_eta[icen], p2_eta[icen]);
      //   farray_Eta->Add(f1weightEta[icen]);
      // }

      const Int_t Ncen_Eta = 11;
      const Double_t centrality_Eta[Ncen_Eta] = {0,5,10,20,30,40,50,60,70,80,100};
      TArrayD *centarray_Eta = new TArrayD(Ncen_Eta,centrality_Eta);

      TObjArray *farray_Eta = new TObjArray(Ncen_Eta-1);
      TF1 *f1weightEta[Ncen_Eta-1];

      const Double_t p0[Ncen_Eta-1] = {2.39991e+02, 1.78111e+02, 1.21109e+02, 6.91786e+01, 4.07880e+01, 2.05577e+01, 1.08079e+01, 4.98463e+00, 2.23119e+00, 1.16590e+00};//Ae
      const Double_t p1[Ncen_Eta-1] = {3.84238e-01, 3.90424e-01, 3.96450e-01, 4.05978e-01, 4.11701e-01, 4.22538e-01, 4.23961e-01, 4.31410e-01, 4.28510e-01, 3.85667e-01};//Te
      const Double_t p2[Ncen_Eta-1] = {1.16561e+03, 9.25084e+02, 7.15782e+02, 5.01254e+02, 3.55730e+02, 2.36175e+02, 1.49593e+02, 8.91887e+01, 4.56253e+01, 1.97146e+01};//A
      const Double_t p3[Ncen_Eta-1] = {3.06202e-01, 3.12875e-01, 3.13599e-01, 3.14247e-01, 3.05660e-01, 3.01051e-01, 2.89748e-01, 2.75489e-01, 2.64481e-01, 2.50078e-01};//T
      const Double_t p4[Ncen_Eta-1] = {2.73068e+00, 2.72595e+00, 2.70812e+00, 2.68939e+00, 2.64558e+00, 2.61370e+00, 2.56506e+00, 2.53088e+00, 2.51428e+00, 2.43636e+00};//n

      for(Int_t icen=0;icen<Ncen_Eta-1;icen++){
        f1weightEta[icen] = new TF1(Form("f1weightEta_%d",icen),"0.48 * ([0] * TMath::Exp(-(TMath::Sqrt(x*x + 0.547*0.547) - 0.547) / [1]) + [2] * TMath::Power(1 + (x*x)/([3]*[3]*[4]) , -[4]))",0,100);//TCM fit to PbPb
        f1weightEta[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        f1weightEta[icen]->SetNpx(1000);
        farray_Eta->Add(f1weightEta[icen]);
      }

      task->SetAdditionalEtaPtWeightFunction(centarray_Eta, farray_Eta);

    }

    if (OptionSet.Contains("K0SWeight")){
      const Int_t Ncen_K0S = 7;
      const Double_t centrality_K0S[Ncen_K0S] = {0,5,10,20,40,60,100};
      TArrayD *centarray_K0S = new TArrayD(Ncen_K0S,centrality_K0S);

      TObjArray *farray_K0S = new TObjArray(Ncen_K0S-1);
      TF1 *f1weightK0S[Ncen_K0S-1];
      const Double_t p0_K0S[Ncen_K0S-1] = {  1.81,   1.83,   1.81,   1.70,   1.91,   1.85};
      const Double_t p1_K0S[Ncen_K0S-1] = {  1.38,   1.31,   1.27,   1.30,   1.61,   1.06};
      const Double_t p2_K0S[Ncen_K0S-1] = {-0.373, -0.398, -0.424, -0.565, -0.640, -0.714};
      const Double_t p3_K0S[Ncen_K0S-1] = {  3.93,   3.70,   4.14,   4.88,  0.542,   1.30};
      const Double_t p4_K0S[Ncen_K0S-1] = {-0.235, -0.279, -0.435,  -2.56, -0.356, -0.543};

      for(Int_t icen = 0; icen < Ncen_K0S-1; icen++){
        f1weightK0S[icen] = new TF1(Form("f1weightK0S_%d", icen), "[0] * (2/(1+exp(-[1]*x)) - 1) - ( 0 + [2]/(exp( -(x-[3]) / [4] ) + 1) )",0,100);
        f1weightK0S[icen]->SetParameters(p0_K0S[icen], p1_K0S[icen], p2_K0S[icen], p3_K0S[icen], p4_K0S[icen]);
        farray_K0S->Add(f1weightK0S[icen]);
      }

      task->SetAdditionalK0SPtWeightFunction(centarray_K0S, farray_K0S);
    }

    if (OptionSet.Contains("GammaWeight")){

      // const Int_t Ncen_Pi0 = 4;
      // const Double_t centrality_Pi0[Ncen_Pi0] = {0,10,30,50};
      // TArrayD *centarray_Pi0 = new TArrayD(Ncen_Pi0, centrality_Pi0);

      // TObjArray *farray_Pi0 = new TObjArray(Ncen_Pi0-1);
      // TF1 *f1weightPi0[Ncen_Pi0-1];

      // const Double_t p0_Pi0[Ncen_Pi0-1] = {0.15, 0.15, 0.499};
      // const Double_t p1_Pi0[Ncen_Pi0-1] = {1.345, 1.345, 2.29047};
      // const Double_t p2_Pi0[Ncen_Pi0-1] = {0.998, 0.998, 0.967643};
      // const Double_t p3_Pi0[Ncen_Pi0-1] = {0.374, 0.374, 0.369492};

      // for(Int_t icen = 0; icen < Ncen_Pi0-1; icen++){
      //   f1weightPi0[icen] = new TF1(Form("f1weightPi0_%d",icen), "[0] + [1]*pow(x,[2])*exp(-[3]*x*x)", 0, 100);
      //   f1weightPi0[icen]->SetParameters(p0_Pi0[icen], p1_Pi0[icen], p2_Pi0[icen], p3_Pi0[icen]);
      //   farray_Pi0->Add(f1weightPi0[icen]);
      // }

      //for gamma pT weighting
      const Int_t Ncen_Gamma = 11;
      const Double_t centrality_Gamma[Ncen_Gamma] = {0,5,10,20,30,40,50,60,70,80,100};
      TArrayD *centarray_Gamma = new TArrayD(Ncen_Gamma,centrality_Gamma);

      TObjArray *farray_Gamma = new TObjArray(Ncen_Gamma-1);
      TF1 *f1weightGamma[Ncen_Gamma-1];

      const Double_t p0[Ncen_Gamma-1] = {4.12097e+02 , 2.72720e+02, 1.85512e+02, 6.40594e+01, 6.40594e+01, 6.71162e+00, 6.71162e+00, 5.94487e+02, 5.94487e+02, 5.94487e+02};//Ae
      const Double_t p1[Ncen_Gamma-1] = {3.19655e-01 , 3.28095e-01, 3.31559e-01, 3.42014e-01, 3.42014e-01, 4.27787e-02, 4.27787e-02, 1.16145e-01, 1.16145e-01, 1.16145e-01};//Te
      const Double_t p2[Ncen_Gamma-1] = {5.26532e+03 , 3.86198e+03, 2.22781e+03, 8.59703e+02, 8.59703e+02, 1.74914e+02, 1.74914e+02, 1.32598e+01, 1.32598e+01, 1.32598e+01};//A
      const Double_t p3[Ncen_Gamma-1] = {2.25277e-01 , 2.38582e-01, 2.56551e-01, 2.93980e-01, 2.93980e-01, 3.75784e-01, 3.75784e-01, 4.84576e-01, 4.84576e-01, 4.84576e-01};//T
      const Double_t p4[Ncen_Gamma-1] = {2.76231e+00 , 2.79069e+00, 2.80401e+00, 2.86491e+00, 2.86491e+00, 3.02951e+00, 3.02951e+00, 3.09311e+00, 3.09311e+00, 3.09311e+00};//n

      for(Int_t icen=0;icen<Ncen_Gamma-1;icen++){
        f1weightGamma[icen] = new TF1(Form("f1weightGamma_%d",icen),"[0] * TMath::Exp(-(TMath::Sqrt(x*x + 0*0) - 0) / [1]) + [2] * TMath::Power(1 + (x*x)/([3]*[3]*[4]) , -[4])",0,100);//TCM fit to PbPb
        f1weightGamma[icen]->SetParameters(p0[icen],p1[icen],p2[icen],p3[icen],p4[icen]);
        f1weightGamma[icen]->SetNpx(1000);
        farray_Gamma->Add(f1weightGamma[icen]);
      }

      task->SetAdditionalGammaPtWeightFunction(centarray_Gamma, farray_Gamma);//do not change pi0 spectra in MC
    }

  }

  mgr->AddTask(task);
  mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());

  TString outputFile = AliAnalysisManager::GetCommonFileName();
  TString outputList = taskname;
  if (TString(subname).Contains("_"))
    outputList += Form("%s", subname);
  else
    outputList += Form("_%s", subname);

  AliAnalysisDataContainer* coutput = mgr->CreateContainer(outputList.Data(),
                                                           TList::Class(),
                                                           AliAnalysisManager::kOutputContainer,
                                                           outputFile.Data());
  mgr->ConnectOutput(task, 1, coutput);

  return task;
}