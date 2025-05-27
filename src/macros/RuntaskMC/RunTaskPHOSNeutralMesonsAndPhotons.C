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

R__LOAD_LIBRARY(AliAnalysisPHOSNeutralMesonsAndPhotons_cxx)
#endif

void RunTaskPHOSNeutralMesonsAndPhotons(const char* dataset = "collection.xml"){

  gErrorIgnoreLevel = 2001;  
  
  TChain* chain = new TChain("aodTree");

  // // Connect to alien
  TGrid::Connect("alien://");

  // // cout << "Pi0Analysis: processing collection " << dataset << endl;

  TXMLEngine a ;
  XMLDocPointer_t doc = a.ParseFile(dataset);
  XMLNodePointer_t rn= a.DocGetRootElement(doc);
  XMLNodePointer_t nn= a.GetChild(rn);
  XMLNodePointer_t nn2= a.GetChild(nn);
  while(nn2){
    XMLNodePointer_t nn3= a.GetChild(nn2);
    const char * rawFile = a.GetAttr(nn3, "turl");
    printf("Processing %s\n", rawFile);
    if(!rawFile) break;
    chain->Add(rawFile);
    nn2= a.GetNext(nn2);
  }

  // Make the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("DirectPhotons_PbPb5TeV");
  
  // ESD input handler
  AliAODInputHandler* esdH = new AliAODInputHandler();
//  esdH->SetReadFriends(kFALSE);
  mgr->SetInputEventHandler(esdH);
  
  // Debug level
//  mgr->SetDebugLevel(2);

  // AddTaskPhysicsSelection(false, true); //Data
  AddTaskPhysicsSelection(true); //MC
  AddTaskMultSelection();
  // AddTaskPIDResponse();


  AddTask_PHOSTender_PCMconfig("PHOSTenderTask", //MC
                               "PHOSTender",
                               "Run2NoNcellCut", //decalib Run2Default
                               1,
                               kTRUE,
                               0,
                               "",
                               "Run2TuneMCNoNcell", //nonlin
                               kTRUE);


  // AddTask_PHOSTender_PCMconfig("PHOSTenderTask", "PHOSTender", "Run2", 1, //Data
  //                              kFALSE, 0, "", "Run2Tune", kFALSE);

  AddTaskPHOSNeutralMesonsAndPhotons("PHOSDirectPhotons", //name
                                     "PbPb", //collision system
                                     kTRUE, //MC flag
                                     AliVEvent::kINT7, //trigger
                                     AliCaloTriggerMimicHelper::kPHOSAny, //specific PHOS trigger
                                     2.5, //Disp cut
                                     2.5, //CPV cut
                                     150., //TOF cut
                                     0.1, //Clust Emin cut
                                     kTRUE, //Clust core E flag
                                     kFALSE, //Use PHOSClusterCuts
                                     kTRUE, //Min bias low E clust cut flag
                                     kFALSE, //Nonlinearity flag
                                     "LHC18q", //LHC period (for nonlin)
                                     20, //N events for mixing
                                     "ClustersQA_CellsQA_Pi0Weight_EtaWeight", //QA options
                                     kTRUE, //do or do not centrality study
                                     "V0M", //centrality estimator
                                     0., //min. centrality
                                     10., //max. centrality
                                     "_CoreE_MBClustCut"); //subname

  AddTaskPHOSNeutralMesonsAndPhotons("PHOSDirectPhotons", //name
                                     "PbPb", //collision system
                                     kTRUE, //MC flag
                                     AliVEvent::kINT7, //trigger
                                     AliCaloTriggerMimicHelper::kPHOSAny, //specific PHOS trigger
                                     2.5, //Disp cut
                                     2.5, //CPV cut
                                     150., //TOF cut
                                     0.1, //Clust Emin cut
                                     kTRUE, //Clust core E flag
                                     kFALSE, //Use PHOSClusterCuts
                                     kTRUE, //Min bias low E clust cut flag
                                     kFALSE, //Nonlinearity flag
                                     "LHC18q", //LHC period (for nonlin)
                                     20, //N events for mixing
                                     "ClustersQA_CellsQA_GammaWeight", //QA options
                                     kTRUE, //do or do not centrality study
                                     "V0M", //centrality estimator
                                     0., //min. centrality
                                     10., //max. centrality
                                     "_CoreE_MBClustCut_GammaWeight"); //subname

  AddTaskPHOSNeutralMesonsAndPhotons("PHOSDirectPhotons", //name
                                     "PbPb", //collision system
                                     kTRUE, //MC flag
                                     AliVEvent::kINT7, //trigger
                                     AliCaloTriggerMimicHelper::kPHOSAny, //specific PHOS trigger
                                     2.5, //Disp cut
                                     2.5, //CPV cut
                                     150., //TOF cut
                                     0.1, //Clust Emin cut
                                     kTRUE, //Clust core E flag
                                     kFALSE, //Use PHOSClusterCuts
                                     kTRUE, //Min bias low E clust cut flag
                                     kFALSE, //Nonlinearity flag
                                     "LHC18q", //LHC period (for nonlin)
                                     20, //N events for mixing
                                     "ClustersQA_CellsQA_Pi0Weight_EtaWeight", //QA options
                                     kTRUE, //do or do not centrality study
                                     "V0M", //centrality estimator
                                     30., //min. centrality
                                     50., //max. centrality
                                     "_CoreE_MBClustCut"); //subname

  AddTaskPHOSNeutralMesonsAndPhotons("PHOSDirectPhotons", //name
                                     "PbPb", //collision system
                                     kTRUE, //MC flag
                                     AliVEvent::kINT7, //trigger
                                     AliCaloTriggerMimicHelper::kPHOSAny, //specific PHOS trigger
                                     2.5, //Disp cut
                                     2.5, //CPV cut
                                     150., //TOF cut
                                     0.1, //Clust Emin cut
                                     kTRUE, //Clust core E flag
                                     kFALSE, //Use PHOSClusterCuts
                                     kTRUE, //Min bias low E clust cut flag
                                     kFALSE, //Nonlinearity flag
                                     "LHC18q", //LHC period (for nonlin)
                                     20, //N events for mixing
                                     "ClustersQA_CellsQA_GammaWeight", //QA options
                                     kTRUE, //do or do not centrality study
                                     "V0M", //centrality estimator
                                     30., //min. centrality
                                     50., //max. centrality
                                     "_CoreE_MBClustCut_GammaWeight"); //subname

  // AddTaskPHOSEtaPhi("PHOSggHBT_kCentral",AliVEvent::kCentral);
 
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local", chain);
  }


}