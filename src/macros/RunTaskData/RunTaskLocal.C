#if defined(__CLING__)
R__ADD_INCLUDE_PATH($ALICE_ROOT)
R__ADD_INCLUDE_PATH($ALICE_PHYSICS)
R__ADD_INCLUDE_PATH(./)
// #include <TAlienCollection.h>
// #include <TAlienResult.h>
#include "AliAODInputHandler.h"
#include "AliVEventHandler.h"
#include <OADB/macros/AddTaskPhysicsSelection.C>
#include <OADB/macros/AddTaskCentrality.C>
#include <OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C>
#include <ANALYSIS/macros/AddTaskPIDResponse.C>
#include <PWGGA/PHOSTasks/PHOS_PbPb/AddAODPHOSTender.C>
#include <PWGGA/GammaConv/macros/AddTask_PHOSTender_PCMconfig.C>
#include <AddTaskPHOSNeutralMesonsAndPhotons.C>
#include "AliAnalysisPHOSNeutralMesonsAndPhotons.h"
#include "PWGGA/PHOSTasks/PHOS_PbPb/AddTaskPHOSEtaPhi.C"
#include <PWGGA/PHOSTasks/PHOS_Run2/macros/AddTaskPHOSPi0EtaToGammaGamma_PbPb_5TeV.C>
#include <PWGGA/PHOSTasks/PHOS_Run2/macros/AddTaskPHOSObjectCreator.C>

R__LOAD_LIBRARY(AliAnalysisPHOSNeutralMesonsAndPhotons_cxx)

#endif

void RunTaskLocal(const char* dataset = "LHC18q_296621.xml", Int_t imin = 1, Int_t imax = 1){

  gErrorIgnoreLevel = 2001;  
  
  TChain* chain = new TChain("aodTree");

  // Connect to alien
  TGrid::Connect("alien://");

  Int_t FileCount = 0;

  TXMLEngine a ;
  XMLDocPointer_t doc = a.ParseFile(dataset);
  XMLNodePointer_t rn = a.DocGetRootElement(doc);
  XMLNodePointer_t nn = a.GetChild(rn);
  XMLNodePointer_t nn2 = a.GetChild(nn);
  while(nn2){
    XMLNodePointer_t nn3 = a.GetChild(nn2);
    const char *rawFile = a.GetAttr(nn3, "turl");
    printf("Processing %s\n", rawFile);
    if(!rawFile) break;
    FileCount++;
    if (FileCount < imin) continue;
    if (FileCount > imax) break;
    chain->Add(rawFile);
    nn2 = a.GetNext(nn2);
  }

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("DirectPhotons_PbPb5TeV");
  
  // ESD input handler
  AliAODInputHandler* esdH = new AliAODInputHandler();
  //  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  AddTaskPhysicsSelection(false, true);
  AddTaskMultSelection();
  AddTaskPIDResponse();


  // AddTask_PHOSTender_PCMconfig("PHOSTenderTask", //MC
  //                          "PHOSTender",
  //                          "Run2Default",
  //                          1,
  //                          kTRUE,
  //                          0,
  //                          "",
  //                          "Default",
  //                          kTRUE);

  AddTask_PHOSTender_PCMconfig("PHOSTenderTask", "PHOSTender", "Run2Default", 1, //Data
                               kFALSE, 0, "", "Run2Tune", kFALSE);

  AddTaskPHOSNeutralMesonsAndPhotons("PHOSDirectPhotons", //name
                                     "PbPb", //collision system
                                     kFALSE, //MC flag
                                     AliVEvent::kCentral, //trigger
                                     AliCaloTriggerMimicHelper::kPHOSAny, //specific PHOS trigger
                                     2.5, //Disp cut
                                     2.5, //CPV cut
                                     150., //TOF cut
                                     0.1, //Clust Emin cut
                                     kTRUE, //Clust core E flag
                                     kTRUE, //Use PHOSClusterCuts
                                     kTRUE, //Min bias low E clust cut flag
                                     kFALSE, //Nonlinearity flag
                                     "LHC18q", //LHC period (for nonlin)
                                     20, //N events for mixing
                                     "ClustersQA_CellsQA", //QA options
                                     kTRUE, //do or do not centrality study
                                     "V0M", //centrality estimator
                                     0., //min. centrality
                                     10., //max. centrality
                                     "_CoreE_MBClustCut_PHOSClusterCuts_NonlinCorr"); //subname

  AddTaskPHOSNeutralMesonsAndPhotons("PHOSDirectPhotons", //name
                                     "PbPb", //collision system
                                     kFALSE, //MC flag
                                     AliVEvent::kSemiCentral, //trigger
                                     AliCaloTriggerMimicHelper::kPHOSAny, //specific PHOS trigger
                                     2.5, //Disp cut
                                     2.5, //CPV cut
                                     150., //TOF cut
                                     0.1, //Clust Emin cut
                                     kTRUE, //Clust core E flag
                                     kTRUE, //Use PHOSClusterCuts
                                     kTRUE, //Min bias low E clust cut flag
                                     kFALSE, //Nonlinearity flag
                                     "LHC18q", //LHC period (for nonlin)
                                     20, //N events for mixing
                                     "ClustersQA_CellsQA", //QA options
                                     kTRUE, //do or do not centrality study
                                     "V0M", //centrality estimator
                                     30., //min. centrality
                                     50., //max. centrality
                                     "_CoreE_MBClustCut_PHOSClusterCuts_NonlinCorr"); //subname

  AddTaskPHOSObjectCreator("PHOSObjectCreatorl_kCentral", AliVEvent::kCentral, kTRUE, kFALSE, 150., kFALSE, kFALSE, "LHC18q", kFALSE);

  AddTaskPHOSPi0EtaToGammaGamma_PbPb_5TeV("Pi0EtaToGammaGamma", //name
                                               AliVEvent::kCentral, // trigger
                                               "PbPb", // system
                                               kFALSE, // MC flag
                                               -1, //L1H,L1M,L1L PHOS trigger
                                               -1, //L0 PHOS trigger
                                               0, // min. centrality
                                               10, // max. centrality
                                               10, // N events for mixing
                                               kFALSE, // flow task flag
                                               -1, // harmonics
                                               -1, // flow method
                                               -1, // Qn detectors
                                               kTRUE, // use core E
                                               kTRUE, // use core Disp
                                               2.5, // CPV cut
                                               2.5, // Disp cut
                                               kTRUE, // use PHOS tender
                                               kFALSE, // TOFEff correction
                                               kFALSE, // TrigEff correction
                                               kFALSE, // nolin study
                                               150., //TOF cut
                                               -1, // minimum distance to bad channel
                                               0.2, // min clust E
                                               kFALSE, // isJJMC
                                               "MBMC", // MC type
                                               kFALSE, // ForceActiveTRU
                                               kFALSE,  // ApplyTOFTrigger
                                               AliPHOSEventCuts::kMultiVertexer, // PileupFinder
                                               "kCentral"); // subname

  AddTaskPHOSObjectCreator("PHOSObjectCreatorl_kSemiCentral", AliVEvent::kSemiCentral, kTRUE, kFALSE, 150., kFALSE, kFALSE, "LHC18q", kFALSE);

  AddTaskPHOSPi0EtaToGammaGamma_PbPb_5TeV("Pi0EtaToGammaGamma", //name
                                               AliVEvent::kSemiCentral, // trigger
                                               "PbPb", // system
                                               kFALSE, // MC flag
                                               -1, //L1H,L1M,L1L PHOS trigger
                                               -1, //L0 PHOS trigger
                                               30, // min. centrality
                                               50, // max. centrality
                                               10, // N events for mixing
                                               kFALSE, // flow task flag
                                               -1, // harmonics
                                               -1, // flow method
                                               -1, // Qn detectors
                                               kTRUE, // use core E
                                               kTRUE, // use core Disp
                                               2.5, // CPV cut
                                               2.5, // Disp cut
                                               kTRUE, // use PHOS tender
                                               kFALSE, // TOFEff correction
                                               kFALSE, // TrigEff correction
                                               kFALSE, // nolin study
                                               150., //TOF cut
                                               -1, // minimum distance to bad channel
                                               0.2, // min clust E
                                               kFALSE, // isJJMC
                                               "MBMC", // MC type
                                               kFALSE, // ForceActiveTRU
                                               kFALSE,  // ApplyTOFTrigger
                                               AliPHOSEventCuts::kMultiVertexer, // PileupFinder
                                               "kSemiCentral"); // subname

  // AddTaskPHOSEtaPhi("PHOSggHBT_kCentral",AliVEvent::kCentral);
 
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }


}