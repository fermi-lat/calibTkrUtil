#include <map>
#include <algorithm>
#include <utility>
#include <iterator>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <string>
//#include "ElecToGeo.h"
#include "RootTestAnalyzer.h"
#include "TROOT.h"
//#include "ToString.h"
//#include "enums/TriggerBits.h"
//#include "timestamps.h"

using std::string;
using std::vector;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;

RootAnalyzer::RootAnalyzer() : m_outputFile(0),
                               m_mcChain(0), m_mcBranch(0), 
                               //m_mcEvent(0),  m_tree(0), m_branch(0),
                               m_reconChain(0), m_reconBranch(0), 
                               m_reconEvent(0), m_digiChain(0),m_digiBranch(0),
                               m_digiEvent(0)
{
  // make sure unsigned int hold 32 bit data 
  assert(sizeof(unsigned int) == 4);

    // initialize ROOT if not already done
  if(gROOT == 0) {
    static TROOT gRoot("RootAnalyzer", "RootAnalyzer");
  }

  m_mcChain = new TChain("Mc");
  m_digiChain = new TChain("Digi");
  m_reconChain = new TChain("Recon");

  m_tkrHits = new TkrHits( true );

  m_tkrNoiseOcc = new TkrNoiseOcc();

}

RootAnalyzer::~RootAnalyzer()
{
  // since root will do the garbage collection automatically for objects such
  // as TTree and TH1, we don't want to deallocate them once more
  delete m_tkrHits;
  delete m_tkrNoiseOcc;
}

void RootAnalyzer::produceOutputFile()
{
  //  TDirectory* saveDir = gDirectory;

  if(m_outputFile) {
    m_outputFile->cd();
    m_outputFile->cd("TkrHits");
    m_tkrHits->saveAllHist();
    m_tkrNoiseOcc->writeAnaToHis(m_tkrNoiseOcc_dir);
    //m_outputFile->cd();
    //m_outputFile->Write(0, TObject::kOverwrite);
    m_outputFile->Close();
  }

  /*
  if(m_histFile) {
    m_histFile->cd();
    m_histFile->Write(0, TObject::kOverwrite);
    m_histFile->Close();
  }
  */

  //  saveDir->cd();
}


bool RootAnalyzer::isRootFile(const string& f)
{
  TFile *rootF = TFile::Open( f.c_str() );
  if(rootF->IsZombie()) {
    return false;
  }
  else {
    return true;
  }
}

void RootAnalyzer::makeTChain(const string& line, TChain* chain)
{
  std::istringstream stream(line);
  string f;
  while(stream >> f) {
    if(isRootFile(f)) { 
      chain->Add(f.c_str());
    }
    else {
      cout << f << " does not exist or is not a root file! It will not be used to make the SVAC ntuple and histogram files!" << endl;
    } 
  }
}


void RootAnalyzer::setOutputRootFile( const char* rootF ){
  cout << "Output TKR root file: " << rootF << endl;
  m_outputFile = new TFile(rootF, "RECREATE");
  m_outputFile->mkdir("TkrHits");
  m_tkrNoiseOcc_dir = m_outputFile->mkdir("TkrNoiseOcc");
}


void RootAnalyzer::setInputRootFiles( const string& mcF, const string& digiF, const string& reconF ){
  cout << "Input mc file(s): " << mcF << endl;
  if( mcF != "none" ) makeTChain(mcF, m_mcChain);
  cout << "Input digi file(s): " << digiF << endl;
  makeTChain(digiF, m_digiChain);
  cout << "Input recon file(s): " << reconF << endl;
  makeTChain(reconF, m_reconChain);
}


void RootAnalyzer::parseLine(const string& line, string& str)
{
  std::istringstream stream(line);
  stream >> str;
}

void RootAnalyzer::parseOptionFile(const char* f)
{
  ifstream optF(f);
  string line, mcF, digiF, reconF;

  while( getline(optF, line) ) {
    if(!isEmptyOrCommentStr(line)) break;
  }
  mcF = line;
  
  while( getline(optF, line) ) {
    if(!isEmptyOrCommentStr(line)) break;
  }
  digiF = line;
  
  while( getline(optF, line) ) {
    if(!isEmptyOrCommentStr(line)) break;
  }
  reconF = line;

  setInputRootFiles( mcF, digiF, reconF );

  while( getline(optF, line) ) {
    if(!isEmptyOrCommentStr(line)) break;
  }
  string rootF;
  parseLine(line, rootF);
  setOutputRootFile( rootF.c_str() );

  /*
  m_tree = new TTree("Output", "Root Analyzer");

  // Set max file size to 500 GB:
  Long64_t maxTreeSize = 5000000000000;
  m_tree->SetMaxTreeSize(maxTreeSize);

  // create branches for each ntuple variable
  //  createBranches();

  while( getline(optF, line) ) {
    if(!isEmptyOrCommentStr(line)) break;
  }
  string histF;
  parseLine(line, histF);
  cout << "Output hist file: " << histF << endl;
  m_histFile = new TFile(histF.c_str(), "RECREATE");

  // create histograms for strip hits

  for(int iTower = 0; iTower != g_nTower; ++iTower) {
    for(int iLayer = 0; iLayer != g_nTkrLayer; ++iLayer) {
      for(int iView = 0; iView != g_nView; ++iView) {

        char name1[] = "hit00000";
        sprintf(name1, "hit%02d%02d%1d", iTower, iLayer, iView);

        m_stripHits[iTower][iLayer][iView] = 
          new TH1F(name1, name1, g_nStripsPerLayer, 0, g_nStripsPerLayer);

        char name2[] = "map00000";
        sprintf(name2, "map%02d%02d%1d", iTower, iLayer, iView);

        m_stripMap[iTower][iLayer][iView] = 
          new TH2F(name2, name2, g_nFEC, 0, g_nFEC, g_nStripsPerLayer/g_nFEC,
                   0, g_nStripsPerLayer/g_nFEC);

      }
    }
  }
  */
}

bool RootAnalyzer::isEmptyOrCommentStr(const string& s)
{
  int len = s.length();
  if(len == 0) return true;
  bool empty = true;
  for(int i = 0; i != len; ++i) {
    if(s[0] != ' ') empty = false;
  }
  if(empty) return true;

  if(len >= 2 && s[0] == '/' && s[1] == '/') {
    return true;
  }
  else {
    return false;
  }
}

void RootAnalyzer::analyzeData( int numEvents = 10000 )
{
  Long64_t nMc = m_mcChain->GetEntries();
  cout << "No. of Mc events to be processed: " << nMc << endl;

  Long64_t nDigi = m_digiChain->GetEntries();
  cout << "No. of Digi events to be processed: " << nDigi << endl;

  Long64_t nRecon = m_reconChain->GetEntries();
  cout << "No. of Recon events to be processed: " << nRecon << endl;

  Long64_t nEvent = std::max(std::max(nMc, nDigi), nRecon);
  if( (nEvent != nMc && nMc != 0) || (nEvent != nDigi && nDigi != 0) ||
      (nEvent != nRecon && nRecon != 0) ) {
    cout << "No. of events in mc, digi or recon files are not identical with each other, stop processing!" << endl;
    exit(1);
  }

//   if(nMc != 0) {
//     m_mcEvent = 0;  
//     // what is stored in the root tree is actually a pointer rather than
//     // mc event
//     m_mcChain->SetBranchAddress("McEvent", &m_mcEvent);
//   }

  if(nRecon != 0) {
    m_reconEvent = 0;
    m_reconChain->SetBranchAddress("ReconEvent", &m_reconEvent);
    //m_reconChain->SetBranchStatus("m_acd",0);
  }

  if(nDigi != 0) {
    m_digiEvent = 0;
    m_digiChain->SetBranchAddress("DigiEvent", &m_digiEvent);
    m_digiChain->SetBranchStatus("*",0);
    m_digiChain->SetBranchStatus("m_eventId", 1);
    m_digiChain->SetBranchStatus("m_runId", 1);
    m_digiChain->SetBranchStatus("m_timeStamp", 1);
    m_digiChain->SetBranchStatus("m_metaEvent", 1);
    m_digiChain->SetBranchStatus("m_tkrDigiCol", 1);
    m_digiChain->SetBranchStatus("m_gem", 1);
  }

  if( numEvents>0 && numEvents<nEvent  ) nEvent = numEvents;
  cout << "# of events to be processed: " << nEvent << std::endl;

  m_tkrHits->setNevents(nEvent);
  m_tkrHits->setOutputFile(m_outputFile);
  //Load current event pointers in TkrCalibManager
  m_tkrHits->setEventPtrs(m_digiEvent, m_reconEvent);
  
  //TkrNoiseOcc::initAnalysis(int nEvent, int evt_interval, int coincidence_cut, int multi_ld, int multi_hd)
  m_tkrNoiseOcc->initAnalysis();
  m_tkrNoiseOcc->setDigiEvtPtr(m_digiEvent);

  time_t startTime, endTime, cTime;
  time( &startTime );
  Long64_t dEvent = 100;

  for(Long64_t  iEvent = 0; iEvent != nEvent; ++iEvent) {

    //    m_ntuple.reset();  
    //    if(m_mcEvent) m_mcEvent->Clear();
    if(m_digiEvent) m_digiEvent->Clear();
    if(m_reconEvent) m_reconEvent->Clear();

    if(nMc != 0) {
      m_mcChain->GetEvent(iEvent);
//       analyzeMcTree();
    }

     if(nDigi != 0) {
       m_digiChain->GetEvent(iEvent);
       //       analyzeDigiTree();
     }

     if(nRecon != 0) {
       m_reconChain->GetEvent(iEvent);
       //       analyzeReconTree();
     }

//     analyzeTot();
    
    //Tracker Calibration Analysis
    m_tkrHits->analyzeEvent();
    m_tkrNoiseOcc->anaDigiEvt();

    //fillOutputTree();
//     if(m_mcEvent) m_mcEvent->Clear();
//     if(m_digiEvent) m_digiEvent->Clear();
//     if(m_reconEvent) m_reconEvent->Clear();

    //
    // show status
    //
    if( iEvent%dEvent == 0 ){
      if( iEvent/dEvent >= 10 )
        if( nEvent/dEvent > 100 ) dEvent *= 10;
        else if( nEvent/dEvent > 50 ) dEvent *= 5;
        else if( nEvent/dEvent > 20 ) dEvent *= 2;
      time( &cTime );
      float dt = cTime - startTime;
      if( dt < 10 ) continue;
      std::cout << iEvent << " events processed in "
                << dt << " s, "
                << int(iEvent/dt) << " events/s, estimated time: "
                << int( (nEvent-iEvent)*dt/iEvent ) << " s"
                << std::endl;
    }


  }  
  time( &endTime );
  if( endTime==startTime )
    std::cout << "total # of events: " << nEvent 
                << " in " << (endTime-startTime) << " s, "
                << std::endl;
  else
    std::cout << "total # of events: " << nEvent 
              << " in " << (endTime-startTime) << " s, "
              << nEvent/(endTime-startTime) << " events/s"
              << std::endl;
  
}

//****************************************************************
//
// following are remant from original code.
//
//****************************************************************


// void RootAnalyzer::analyzeMcTree()
// {
//     m_ntuple.m_runId = (int) m_mcEvent->getRunId();

//     m_ntuple.m_seqNo = (int) m_mcEvent->getSequence();

//     // get info about primary particle
//     McParticle* particle = m_mcEvent->getMcParticle(0);
//     if(particle) {
//       m_ntuple.m_parId = particle->getParticleId();
//       m_ntuple.m_mcEnergy = particle->getInitialFourMomentum().E();
//       const TVector3& pos = particle->getInitialPosition();
//       m_ntuple.m_startPos[0] = (float) pos.X();
//       m_ntuple.m_startPos[1] = (float) pos.Y();
//       m_ntuple.m_startPos[2] = (float) pos.Z();
//       const TLorentzVector& lor = particle->getInitialFourMomentum();
//       float norm = sqrt(lor(0)*lor(0) + lor(1)*lor(1) + lor(2)*lor(2));
//       m_ntuple.m_startDir[0] = lor(0)/norm;
//       m_ntuple.m_startDir[1] = lor(1)/norm;
//       m_ntuple.m_startDir[2] = lor(2)/norm;
//     }

//     // get gamma conversion (pair production) point
//     int nPar = m_mcEvent->getMcParticleCount();

//     // first particle produced by pair production, used in filling m_pairEne
//     int iPair = 0;

//     double x0=-9999.0, y0=-9999.0, z0=-9999.0, x1=-9999.0, y1=-999.0, z1=-9999.0;

//     for(int iPar = 0; iPar != nPar; ++iPar) {
//       McParticle* par = m_mcEvent->getMcParticle(iPar);
//       if(par && par->getProcess() == "conv" && iPair == 0) {
//           const TVector3& pos = par->getInitialPosition();
//           m_ntuple.m_convPos[0] = (float) pos.X();
//           m_ntuple.m_convPos[1] = (float) pos.Y();
//           m_ntuple.m_convPos[2] = (float) pos.Z();
//           m_ntuple.m_pairEne[iPair++] = par->getInitialFourMomentum().Energy();
//           x0 = par->getInitialFourMomentum().Px();
//           y0 = par->getInitialFourMomentum().Py();
//           z0 = par->getInitialFourMomentum().Pz();
//           continue;
//       }
//       if(par && par->getProcess() == "conv" && iPair == 1) { 
//         //second particle produced by pair production
//           m_ntuple.m_pairEne[iPair++] = par->getInitialFourMomentum().Energy();
//           x1 = par->getInitialFourMomentum().Px();
//           y1 = par->getInitialFourMomentum().Py();
//           z1 = par->getInitialFourMomentum().Pz();
//       }
//     }

//     if(iPair == 2) { // there is a conversion
//       double a = sqrt(x0*x0 + y0*y0 + z0*z0);
//       double b = sqrt(x1*x1 + y1*y1 + z1*z1);
//       m_ntuple.m_convAngle = acos((x0*x1 + y0*y1 + z0*z1) / (a * b) );
//     }

//     // fill in energy deposited in each TKR layer
//     int nPosHit = m_mcEvent->getMcPositionHitCount();
//     for(int iPosHit = 0; iPosHit != nPosHit; ++iPosHit) {
//       McPositionHit* posHit = m_mcEvent->getMcPositionHit(iPosHit);
//       if(posHit) {
//         VolumeIdentifier id = posHit->getVolumeId();
//         int iTower, iLayer, iView;
//         if( extractTowerLayerView(id, iTower, iLayer, iView) ) {
//             m_ntuple.m_depositEne[iTower][iLayer][iView] += 
//               posHit->getDepositedEnergy();
//         }
//       }
//     } 

//     // fill in energy deposited in CAL
//     int nIntHit = m_mcEvent->getIntegratingHitCount();
//     m_ntuple.m_mcCalEnergy = 0.;
//     for(int iIntHit = 0; iIntHit != nIntHit; ++iIntHit) {
//       McIntegratingHit* intHit = m_mcEvent->getMcIntegratingHit(iIntHit);
//       if(intHit) {
//         VolumeIdentifier id = intHit->getVolumeId();
//         // the 4th id must be 0 if the volume is a CAL crystal
//         if(id[3] == 0) {
//           m_ntuple.m_mcCalEnergy += intHit->getTotalEnergy();
//         }
//       }
//     }
// }
 

// void RootAnalyzer::analyzeReconTree()
// {
//   TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();
 
//   if (tkrRecon == 0) return;

//   TObjArray* siClusterCol = tkrRecon->getClusterCol();
//   if(siClusterCol) {
//     int nClusters = siClusterCol->GetLast()+1;
//     // fill in no of clusters in each layer
//     for(int i = 0; i != nClusters; ++i) {
//       TkrCluster* cluster = dynamic_cast<TkrCluster*>(siClusterCol->At(i));
//       if(cluster) {

//         TowerId tId(cluster->getTkrId().getTowerX(), cluster->getTkrId().getTowerY());
//         int iTower = tId.id();
//         int iLayer = cluster->getLayer();
//         int iView = cluster->getTkrId().getView();

//         assert(iLayer >= 0 && iLayer <= g_nTkrLayer - 1);

//         ++m_ntuple.m_nTkrClusters[iTower][iLayer][iView];
//         ++m_ntuple.m_totalClusters[iTower];

//       }
//     }
//   }

//   m_ntuple.m_nTkrTracks = tkrRecon->getTrackCol()->GetLast()+1;
//   m_ntuple.m_nTkrVertices = tkrRecon->getVertexCol()->GetLast()+1;

//   // fill in info stored in the very first vertex
//   TObjArray* vertices = tkrRecon->getVertexCol();
//   if(m_ntuple.m_nTkrVertices >= 1) {
//     TkrVertex* tkrVertex = dynamic_cast<TkrVertex*>(vertices->At(0));
//     if(tkrVertex) {
//       const TVector3& pos = tkrVertex->getPosition();
//       const TVector3& dir = tkrVertex->getDirection();
//       m_ntuple.m_pos[0] = pos.X();
//       m_ntuple.m_pos[1] = pos.Y();
//       m_ntuple.m_pos[2] = pos.Z();
//       m_ntuple.m_dir[0] = dir.X();
//       m_ntuple.m_dir[1] = dir.Y();
//       m_ntuple.m_dir[2] = dir.Z();
//       m_ntuple.m_fitTotalEnergy = tkrVertex->getEnergy();
//       m_ntuple.m_vtxTrks = tkrVertex->getNumTracks();
//     }
//   }

//   // fill in info stored in the very first and second track
//   TObjArray* tracks = tkrRecon->getTrackCol();

//   for(int i = 0; i != 2; ++i) {
//     TkrTrack* tkrTrack = dynamic_cast<TkrTrack*>(tracks->At(i));
//     if(tkrTrack) {
//       m_ntuple.m_nFit[i] += (int) tkrTrack->getNumFitHits();
//       m_ntuple.m_chi2[i] = tkrTrack->getChiSquareFilter();
//       m_ntuple.m_chi2Smooth[i] = tkrTrack->getChiSquareSmooth();
//       m_ntuple.m_rms[i] = tkrTrack->getScatter();
//       m_ntuple.m_msAngle[i] = tkrTrack->getKalThetaMS();
//       m_ntuple.m_tkrEnergy[i] = tkrTrack->getKalEnergy();


//       // End-of-track parameters:
//       TkrTrackHit* hit = (TkrTrackHit*) tkrTrack->Last();
//       TVector3 endPos = hit->getPoint(TkrTrackHit::SMOOTHED);
//       TVector3 endSlopeTmp = hit->getDirection(TkrTrackHit::SMOOTHED);
//       TVector3 endSlope = endSlopeTmp.Unit();

//       if (i == 0) {
//         m_ntuple.m_tkr1EndPos[0] = endPos.x();
//         m_ntuple.m_tkr1EndPos[1] = endPos.y();
//         m_ntuple.m_tkr1EndPos[2] = endPos.z();
//         // Need -1 here .....
//         m_ntuple.m_tkr1EndDir[0] = -1.0*endSlope.x();
//         m_ntuple.m_tkr1EndDir[1] = -1.0*endSlope.y();
//         m_ntuple.m_tkr1EndDir[2] = -1.0*endSlope.z();
//       } else {
//         m_ntuple.m_tkr2EndPos[0] = endPos.x();
//         m_ntuple.m_tkr2EndPos[1] = endPos.y();
//         m_ntuple.m_tkr2EndPos[2] = endPos.z();
//         // Need -1 here .....
//         m_ntuple.m_tkr2EndDir[0] = -1.0*endSlope.x();
//         m_ntuple.m_tkr2EndDir[1] = -1.0*endSlope.y();
//         m_ntuple.m_tkr2EndDir[2] = -1.0*endSlope.z();
//       }
//     }
//   }

//   // calculate energy measured in calorimeter
//   CalRecon* calRecon = m_reconEvent->getCalRecon();
//   if(calRecon) {
//     TObjArray* calClusterCol = calRecon->getCalClusterCol();
//     if(calClusterCol) {
//       // currently there is just one cluster in CAL
//       CalCluster* calCluster = dynamic_cast<CalCluster*>(calClusterCol->At(0));
//       if (calCluster) {
//         m_ntuple.m_calPos[0] = calCluster->getParams().getCentroid().x();
//         m_ntuple.m_calPos[1] = calCluster->getParams().getCentroid().y();
//         m_ntuple.m_calPos[2] = calCluster->getParams().getCentroid().z();
//         m_ntuple.m_calEnergy = calCluster->getParams().getEnergy();
//       }
//     }

//     TObjArray* calXtalRecCol = calRecon->getCalXtalRecCol();

//     if(calXtalRecCol) {
//       int nCalRec = calXtalRecCol->GetLast() + 1;

//       for(int i = 0; i != nCalRec; ++i) {
//         CalXtalRecData* calData = 
//           dynamic_cast<CalXtalRecData*>(calXtalRecCol->At(i));
//         if(calData) {

//           CalXtalId id = calData->getPackedId();
//           int iTower = id.getTower();
//           int iLayer = id.getLayer();
//           int iCol = id.getColumn();

//           float ene = calData->getEnergy();
          
//           float eneNeg = calData->getEnergy(0,CalXtalId::NEG);
//           float enePos = calData->getEnergy(0,CalXtalId::POS);

//           if(ene >= 0) ++(m_ntuple.m_nCrystalHit[iTower]);

//           // CAL layer end energies:
//           m_ntuple.m_xtalEne[iTower][iLayer][iCol][0] = enePos;
//           m_ntuple.m_xtalEne[iTower][iLayer][iCol][1] = eneNeg;

//           if(ene > m_ntuple.m_maxCalEnergy) m_ntuple.m_maxCalEnergy = ene;

//           TVector3 pos =  calData->getPosition();
//           m_ntuple.m_xtalPos[iTower][iLayer][iCol][0] = pos.x();
//           m_ntuple.m_xtalPos[iTower][iLayer][iCol][1] = pos.y();
//           m_ntuple.m_xtalPos[iTower][iLayer][iCol][2] = pos.z();

//         }
//       }
//     }

//     TObjArray* calMipTrackCol = calRecon->getCalMipTrackCol();
//     if (calMipTrackCol) {
//       int nCalMip = calMipTrackCol->GetLast() + 1;
//       m_ntuple.m_calMipNum = nCalMip;

//       for (int i = 0; i< std::min(2,nCalMip); i++) {
//         CalMipTrack* calMip = dynamic_cast<CalMipTrack*>(calMipTrackCol->At(i));
//         if (calMip) {
//           TVector3 pos =  calMip->getPoint();
//           TVector3 dir =  calMip->getDir();

//           if (i == 0) {
//             m_ntuple.m_calMip1Pos[0] = pos.x();
//             m_ntuple.m_calMip1Pos[1] = pos.y();
//             m_ntuple.m_calMip1Pos[2] = pos.z();

//             m_ntuple.m_calMip1Dir[0] = dir.x();
//             m_ntuple.m_calMip1Dir[1] = dir.y();
//             m_ntuple.m_calMip1Dir[2] = dir.z();

//             m_ntuple.m_calMip1Chi2    = calMip->getChi2();
//             m_ntuple.m_calMip1D2edge  = calMip->getD2Edge();
//             m_ntuple.m_calMip1ArcLen  = calMip->getArcLen();
//             m_ntuple.m_calMip1Ecor    = calMip->getEcor();
//             m_ntuple.m_calMip1EcorRms = calMip->getEcorRms();
//             m_ntuple.m_calMip1Erm     = calMip->getErm();
//           }
//           if (i == 1) {
//             m_ntuple.m_calMip2Pos[0] = pos.x();
//             m_ntuple.m_calMip2Pos[1] = pos.y();
//             m_ntuple.m_calMip2Pos[2] = pos.z();

//             m_ntuple.m_calMip2Dir[0] = dir.x();
//             m_ntuple.m_calMip2Dir[1] = dir.y();
//             m_ntuple.m_calMip2Dir[2] = dir.z();

//             m_ntuple.m_calMip2Chi2    = calMip->getChi2();
//             m_ntuple.m_calMip2D2edge  = calMip->getD2Edge();
//             m_ntuple.m_calMip2ArcLen  = calMip->getArcLen();
//             m_ntuple.m_calMip2Ecor    = calMip->getEcor();
//             m_ntuple.m_calMip2EcorRms = calMip->getEcorRms();
//             m_ntuple.m_calMip2Erm     = calMip->getErm();
//           }
//         }
//       }
//     }
//   }  // calRecon


//   //
//   // ACD recon:
//   //
//   AcdRecon* acdRecon = m_reconEvent->getAcdRecon();

//   //
//   if (acdRecon) {
//     UInt_t nAcdHit = acdRecon->nAcdHit();
//     for ( UInt_t iAcdHit(0); iAcdHit < nAcdHit; iAcdHit++ ) {
//       const AcdHit* acdHit = acdRecon->getAcdHit(iAcdHit);

//       const AcdId& acdId = acdHit->getId();
//       int acdID = acdId.getId(); 
//       m_ntuple.m_acdMips[acdID][0] = acdHit->getMips(AcdHit::A);
//       m_ntuple.m_acdMips[acdID][1] = acdHit->getMips(AcdHit::B);
//     }
//   }

//   //
//   if (acdRecon) {
//     m_ntuple.m_acdEnergy     = acdRecon->getEnergy();
//     m_ntuple.m_acdTileCount  = acdRecon->getTileCount();
//   }
// }
  
// void RootAnalyzer::analyzeDigiTree()
// {
//   m_ntuple.m_runId = m_digiEvent->getRunId();
//   m_ntuple.m_eventId = m_digiEvent->getEventId();

//   unsigned int word = m_digiEvent->getL1T().getTriggerWord();
//   unsigned bitMask = 0;
//   int ibit = enums::number_of_trigger_bits;
//   while(ibit--) { bitMask |= 1<<ibit; }
//   m_ntuple.m_trigger = word & bitMask;

//   m_ntuple.m_timeStamp = m_digiEvent->getTimeStamp();

//   m_ntuple.m_ebfSecond = m_digiEvent->getEbfTimeSec();
//   m_ntuple.m_ebfNanoSecond =  m_digiEvent->getEbfTimeNanoSec();

//   m_ntuple.m_upperTime   = m_digiEvent->getEbfUpperPpcTimeBase();
//   m_ntuple.m_lowerTime   = m_digiEvent->getEbfLowerPpcTimeBase();
//   m_ntuple.m_timeSeconds = m_digiEvent->getEbfPpcTimeSeconds();

//   m_ntuple.m_summaryWord = m_digiEvent->getEventSummaryData().summary();
//   m_ntuple.m_eventSize   = m_digiEvent->getEventSummaryData().eventSizeInBytes();


//   // GEM information:
//   m_ntuple.m_gemConditionsWord = m_digiEvent->getGem().getConditionSummary();

//   m_ntuple.m_gemLiveTime             = m_digiEvent->getGem().getLiveTime();
//   m_ntuple.m_gemTriggerTime          = m_digiEvent->getGem().getTriggerTime();
//   m_ntuple.m_gemDeltaEventTime       = m_digiEvent->getGem().getDeltaEventTime();
//   m_ntuple.m_gemOnePpsSeconds        = m_digiEvent->getGem().getOnePpsTime().getSeconds();
//   m_ntuple.m_gemOnePpsTime           = m_digiEvent->getGem().getOnePpsTime().getTimebase();
//   m_ntuple.m_gemPrescaled            = m_digiEvent->getGem().getPrescaled();
//   m_ntuple.m_gemDiscarded            = m_digiEvent->getGem().getDiscarded();
//   m_ntuple.m_gemCondArrivalTimeWord  = m_digiEvent->getGem().getCondArrTime().condArr();
//   m_ntuple.m_gemCondArrivalTimeExt   = m_digiEvent->getGem().getCondArrTime().external();
//   m_ntuple.m_gemCondArrivalTimeCno   = m_digiEvent->getGem().getCondArrTime().cno();
//   m_ntuple.m_gemCondArrivalTimeCalLe = m_digiEvent->getGem().getCondArrTime().calLE();
//   m_ntuple.m_gemCondArrivalTimeCalHe = m_digiEvent->getGem().getCondArrTime().calHE();
//   m_ntuple.m_gemCondArrivalTimeTkr   = m_digiEvent->getGem().getCondArrTime().tkr();
//   m_ntuple.m_gemCondArrivalTimeRoi   = m_digiEvent->getGem().getCondArrTime().roi();
//   m_ntuple.m_gemDeltaWindowOpenTime  = m_digiEvent->getGem().getDeltaWindowOpenTime();
//   m_ntuple.m_gemDeadZone             = m_digiEvent->getGem().getMissed();

//   m_ntuple.m_gemAcdTilesXzp = m_digiEvent->getGem().getTileList().getXzp();
//   m_ntuple.m_gemAcdTilesXzm = m_digiEvent->getGem().getTileList().getXzm();
//   m_ntuple.m_gemAcdTilesYzp = m_digiEvent->getGem().getTileList().getYzp();
//   m_ntuple.m_gemAcdTilesYzm = m_digiEvent->getGem().getTileList().getYzm();
//   m_ntuple.m_gemAcdTilesXy  = m_digiEvent->getGem().getTileList().getXy();
//   m_ntuple.m_gemAcdTilesRbn = m_digiEvent->getGem().getTileList().getRbn();
//   m_ntuple.m_gemAcdTilesNa  = m_digiEvent->getGem().getTileList().getNa();

//   unsigned tmpGemTkr   = m_digiEvent->getGem().getTkrVector();
//   unsigned tmpGemRoi   = m_digiEvent->getGem().getRoiVector();
//   unsigned tmpGemCalLe = m_digiEvent->getGem().getCalLeVector();
//   unsigned tmpGemCalHe = m_digiEvent->getGem().getCalHeVector();
//   unsigned tmpGemCno   = m_digiEvent->getGem().getCnoVector();
  
//   m_ntuple.m_triggerTicks = evtTicks(m_ntuple.m_gemTriggerTime, 
//                                                                         m_ntuple.m_gemOnePpsSeconds, m_ntuple.m_gemOnePpsTime,
//                                                                         m_ntuple.m_ebfSecond, m_ntuple.m_ebfNanoSecond);

//   for (int iTower = 0; iTower<g_nTower; iTower++) {
//     m_ntuple.m_gemTkrVector[iTower]   = ((tmpGemTkr >> iTower) & 1) ;      
//     m_ntuple.m_gemRoiVector[iTower]   = ((tmpGemRoi >> iTower) & 1) ;      
//     m_ntuple.m_gemCalLeVector[iTower] = ((tmpGemCalLe >> iTower) &1 ) ;      
//     m_ntuple.m_gemCalHeVector[iTower] = ((tmpGemCalHe >> iTower) &1 ) ;      
//   }
  
//   int iCno = 0;
//   m_ntuple.m_gemCnoVector = ((tmpGemCno >> iCno) & 1) ;      
  
 
//   // Luis's three-in-a-row trigger bits:
//   for (int iTower = 0; iTower<g_nTower; iTower++) {
//     m_ntuple.m_digiTriRowBits[iTower] = m_digiEvent->getL1T().getDigiTriRowBits(iTower);
//     m_ntuple.m_trgReqTriRowBits[iTower] = m_digiEvent->getL1T().getTrgReqTriRowBits(iTower);
//   }
  
//   // Event sizes:
//   for (int iTower = 0; iTower<g_nTower; iTower++) {
//     m_ntuple.m_temLength[iTower] = m_digiEvent->getEventSummaryData().temLength(iTower);
//   }

//   m_ntuple.m_gemLength  = m_digiEvent->getEventSummaryData().gemLength();
//   m_ntuple.m_oswLength  = m_digiEvent->getEventSummaryData().oswLength();
//   m_ntuple.m_aemLength  = m_digiEvent->getEventSummaryData().aemLength();
  
//   for (int iTower = 0; iTower<g_nTower; iTower++) {
//     m_ntuple.m_errLength[iTower]  = m_digiEvent->getEventSummaryData().errLength(iTower);
//     m_ntuple.m_diagLength[iTower] = m_digiEvent->getEventSummaryData().diagLength(iTower);
//   }

//   // Event quality:
//   m_ntuple.m_eventFlags = m_digiEvent->getEventSummaryData().eventFlags();

//   parseDiagnosticData(); 
 

//   //
//   // ACD digi:
//   //
//   const TObjArray* acdDigiCol = m_digiEvent->getAcdDigiCol();
//   assert(acdDigiCol != 0);

//   // Number of digis:
//   int nAcdDigi = acdDigiCol->GetLast()+1;
//   m_ntuple.m_acdNumDigis = nAcdDigi;


//   // Loop over all digis:
//   for(int iDigi = 0; iDigi != nAcdDigi; ++iDigi) {

//     const AcdDigi* acdDigi = dynamic_cast<const AcdDigi*>(acdDigiCol->At(iDigi));
//     assert(acdDigi != 0);

//     // Tile ID:
//     int AcdID = acdDigi->getId().getId();


//     // Attached tile and ID out of bounds?
//     if (acdDigi->getId().getNa()==0 && (AcdID>603 || AcdID<0)) {
//       std::cout << "ACD tile ID for attached tile is >603! - " << AcdID << std::endl;
//     }


//     // Attached tile?
//     if (acdDigi->getId().getNa()==0 && AcdID>-1 && AcdID<604) {

//       m_ntuple.m_acdPha[AcdID][0] = acdDigi->getPulseHeight(AcdDigi::A);
//       m_ntuple.m_acdPha[AcdID][1] = acdDigi->getPulseHeight(AcdDigi::B);

//       m_ntuple.m_acdHitMap[AcdID][0] = acdDigi->getHitMapBit(AcdDigi::A);
//       m_ntuple.m_acdHitMap[AcdID][1] = acdDigi->getHitMapBit(AcdDigi::B);

//       m_ntuple.m_acdRange[AcdID][0] = acdDigi->getRange(AcdDigi::A);
//       m_ntuple.m_acdRange[AcdID][1] = acdDigi->getRange(AcdDigi::B);

//       m_ntuple.m_acdOddParityError[AcdID][0] = acdDigi->getOddParityError(AcdDigi::A);
//       m_ntuple.m_acdOddParityError[AcdID][1] = acdDigi->getOddParityError(AcdDigi::B);

//       m_ntuple.m_acdHeaderParityError[AcdID][0] = acdDigi->getHeaderParityError(AcdDigi::A);
//       m_ntuple.m_acdHeaderParityError[AcdID][1] = acdDigi->getHeaderParityError(AcdDigi::B);

//       m_ntuple.m_acdLowDisc[AcdID][0] = acdDigi->getAcceptMapBit(AcdDigi::A);
//       m_ntuple.m_acdLowDisc[AcdID][1] = acdDigi->getAcceptMapBit(AcdDigi::B);

//       m_ntuple.m_acdTileNumber[AcdID] = acdDigi->getTileNumber();
//       m_ntuple.m_acdMCEnergy[AcdID]   = acdDigi->getEnergy();
//     }
//   }





//   // fill in no of Tkr digis and TOTs
//   m_ntuple.m_nTkrDigis = m_digiEvent->getTkrDigiCol()->GetLast()+1;

//   for(int i = 0; i != m_ntuple.m_nTkrDigis; ++i) {
//     const TkrDigi* tkrDigi = m_digiEvent->getTkrDigi(i);

//     assert(tkrDigi != 0);

//    // fill hit patterns for each plane
//     fillStripHits(tkrDigi);

//     int iTower = tkrDigi->getTower().id();
//     int iLayer = tkrDigi->getBilayer();
//     GlastAxis::axis iView = tkrDigi->getView();
  
//     int nStrips = tkrDigi->getNumHits();

//     m_ntuple.m_totalStripHits[iTower] += nStrips;

//     if(iView == GlastAxis::X) {
//       m_ntuple.m_nStrips[iTower][iLayer][0] = nStrips;
//       m_ntuple.m_tot[iTower][iLayer][0][0] = tkrDigi->getToT(0);
//       m_ntuple.m_tot[iTower][iLayer][0][1] = tkrDigi->getToT(1);
//     }
//     else if(iView == GlastAxis::Y) {
//       m_ntuple.m_nStrips[iTower][iLayer][1] = nStrips;
//       m_ntuple.m_tot[iTower][iLayer][1][0] = tkrDigi->getToT(0);
//       m_ntuple.m_tot[iTower][iLayer][1][1] = tkrDigi->getToT(1);
//     }

//       // fill in corrected tot
//     /*
//       if(m_mcFile == 0) {
//         correctTotDataLinear(tkrDigi);
//         correctTotDataQuad(tkrDigi);
//       }
//     */
//   }

//   /*
//   // ACD digi
//   const TObjArray* acdDigiCol = m_digiEvent->getAcdDigiCol();
//   if (!acdDigiCol) return;

//   int nAcdDigi = acdDigiCol->GetLast() + 1;
//   cout << nAcdDigi << endl;
//   */
// }

/*
void RootAnalyzer::fillOutputTree()
{
  if(m_outputFile) {
    TDirectory* saveDir = gDirectory;
    m_outputFile->cd();
    m_tree->Fill();
    saveDir->cd();
  }
}
*/


// bool RootAnalyzer::extractTowerLayerView(const VolumeIdentifier& id, 
//                                          int& iTower, int& iLayer, 
//                                          int& iView) const
// {
//   // id code is formatted as the following (if only one tower is simulated):
//   // Tower(0)/TowerY(0)/TowerX(0)/TKR(1)/TrayNo(0-18)/Measure(0,1)/View(0,1)/
//   // Ladder/Waffer

//   iTower = id[0];

//   // not Tracker
//   if(id[3] != 1) return 0;

//   int iTray = id[4];
//   iView = id[5];  // Measure X: 0, Measure Y: 1
//   int iBotTop = id[6];  // Bott: 0, Top: 1

//   iLayer = iTray -1 + iBotTop;

//   // check boundary condition
//   if(iTower<0 || iTower>15 || iLayer<0 || iLayer>17 || iView<0 || 
//      iView>1) return 0;

//   return 1;
// }

/*
void RootAnalyzer::fillStripHits(const TkrDigi* tkrDigi)
{

  int tower = tkrDigi->getTower().id();
  int layer = tkrDigi->getBilayer(); 
  int view = tkrDigi->getView();

  int nStrips = tkrDigi->getNumHits();

  TDirectory* saveDir = gDirectory;
  m_histFile->cd();

  for(int i = 0; i != nStrips; ++i) {

    int stripId = tkrDigi->getStrip(i);
    m_stripHits[tower][layer][view]->Fill(stripId);

    int nStripPerFEC = g_nStripsPerLayer / g_nFEC;

    int fec = stripId / nStripPerFEC;
    int res = stripId - fec * nStripPerFEC;

    m_stripMap[tower][layer][view]->Fill(fec, res);

  }

  saveDir->cd();
}
*/


// void RootAnalyzer::analyzeTot()
// {
//   for(int iTower = 0; iTower != g_nTower; ++iTower) {
//     for(int iLayer = (g_nTkrLayer-1); iLayer >= 0; --iLayer) {

//       if(m_ntuple.m_nStrips[iTower][iLayer][0] > 0 &&
//          m_ntuple.m_nStrips[iTower][iLayer][1] > 0) {

//         float totX = std::max(m_ntuple.m_tot[iTower][iLayer][0][0], 
//                               m_ntuple.m_tot[iTower][iLayer][0][1]);
//         float totY = std::max(m_ntuple.m_tot[iTower][iLayer][1][0], 
//                               m_ntuple.m_tot[iTower][iLayer][1][1]);

//         m_ntuple.m_topTot[iTower] = std::max(totX, totY);
//         break;
//       }
//     }

//   }

//   if(m_ntuple.m_nTkrTracks > 0) {

//     TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();
//     TObjArray* vertices = tkrRecon->getVertexCol();
//     if(m_ntuple.m_nTkrVertices >= 1) {
//       TkrVertex* tkrVertex = dynamic_cast<TkrVertex*>(vertices->At(0));
//       if(tkrVertex) {

//         TkrTrack* trk = (TkrTrack*) tkrVertex->getTrack(0);
//         TkrTrackHit* hit = (TkrTrackHit*) trk->First();
//         commonRootData::TkrId id = hit->getClusterPtr()->getTkrId();
//         TowerId tId(id.getTowerX(), id.getTowerX());
//         int convLayer = hit->getClusterPtr()->getLayer();
//         int convTower = tId.id();

//         assert(convLayer >= 0 && convLayer <g_nTkrLayer);
//         float totX = std::max(m_ntuple.m_tot[convTower][convLayer][0][0], 
//                               m_ntuple.m_tot[convTower][convLayer][0][1]);
//         float totY = std::max(m_ntuple.m_tot[convTower][convLayer][1][0], 
//                               m_ntuple.m_tot[convTower][convLayer][1][1]);

//         m_ntuple.m_convTot = std::max(totX, totY);
//       }
//     }
//   }
 
//}

// void RootAnalyzer::parseDiagnosticData()
// {
//   if(m_digiEvent->getTkrDiagnosticCol()) {
//     int nTkrDiag = m_digiEvent->getTkrDiagnosticCol()->GetLast()+1;

//     for(int i = 0; i != nTkrDiag; ++i) {

//       const TkrDiagnosticData* pDiag = m_digiEvent->getTkrDiagnostic(i);
//       int iTower = pDiag->tower();
//       if(m_ntuple.m_diagLength[iTower]) {
//         m_ntuple.m_tpTkr[iTower][pDiag->gtcc()] = pDiag->getDataWord();
//       }
//     }

//     ElecToGeo::getInstance()->decodeTkrTp(m_ntuple.m_tpTkr,m_ntuple.m_tkrReq);
//   }

//   if(m_digiEvent->getCalDiagnosticCol()) {

//     int nCalDiag = m_digiEvent->getCalDiagnosticCol()->GetLast()+1;

//     for(int i = 0; i != nCalDiag; ++i) {

//       const CalDiagnosticData* pDiag = m_digiEvent->getCalDiagnostic(i);

//       int iTower = pDiag->tower();
//       if(m_ntuple.m_diagLength[iTower]) {
//         m_ntuple.m_tpCal[pDiag->tower()][pDiag->layer()] = pDiag->getDataWord();
//       }
//     }

//     ElecToGeo::getInstance()->decodeCalTp(m_ntuple.m_tpCal, m_ntuple.m_calReq,
//                                           m_ntuple.m_calLogAccepts);

//   }
// }

// void RootAnalyzer::createBranches()
// {

//   // Event information:
//   m_tree->Branch("RunID", &(m_ntuple.m_runId), "RunID/i");
//   m_tree->Branch("EventID", &(m_ntuple.m_eventId), "EventID/i");
//   m_tree->Branch("EventSize", &(m_ntuple.m_eventSize), "EventSize/i");
//   m_tree->Branch("EventFlags", &(m_ntuple.m_eventFlags), "EventFlags/i");
//   m_tree->Branch("EvtTime", &(m_ntuple.m_timeStamp), "EvtTime/D");
//   m_tree->Branch("EvtSecond", &(m_ntuple.m_ebfSecond), "EvtSecond/i");
//   m_tree->Branch("EvtNanoSecond", &(m_ntuple.m_ebfNanoSecond), "EvtNanoSecond/i");
//   m_tree->Branch("EvtUpperTime", &(m_ntuple.m_upperTime), "EvtUpperTime/i");
//   m_tree->Branch("EvtLowerTime", &(m_ntuple.m_lowerTime), "EvtLowerTime/i");
//   m_tree->Branch("EvtTimeSeconds", &(m_ntuple.m_timeSeconds),"EvtTimeSeconds/D");
//   m_tree->Branch("EvtTicks", &(m_ntuple.m_triggerTicks),"EvtTicks/D");
//   m_tree->Branch("EvtSummary", &(m_ntuple.m_summaryWord), "EvtSummary/i");

//   // MC information:
//   m_tree->Branch("McSeqNo", &(m_ntuple.m_seqNo), "McSeqNo/i");
//   m_tree->Branch("McId", &(m_ntuple.m_parId), "McId/I");
//   m_tree->Branch("McTotalEnergy", &(m_ntuple.m_mcEnergy), "McTotalEnergy/F");
//   m_tree->Branch("McX0", &(m_ntuple.m_startPos[0]), "McX0/F");
//   m_tree->Branch("McY0", &(m_ntuple.m_startPos[1]), "McY0/F");
//   m_tree->Branch("McZ0", &(m_ntuple.m_startPos[2]), "McZ0/F");
//   m_tree->Branch("McXDir", &(m_ntuple.m_startDir[0]), "McXDir/F");
//   m_tree->Branch("McYDir", &(m_ntuple.m_startDir[1]), "McYDir/F");
//   m_tree->Branch("McZDir", &(m_ntuple.m_startDir[2]), "McZDir/F");
//   m_tree->Branch("McConvPointX", &(m_ntuple.m_convPos[0]), "McConvPointX/F");
//   m_tree->Branch("McConvPointY", &(m_ntuple.m_convPos[1]), "McConvPointY/F");
//   m_tree->Branch("McConvPointZ", &(m_ntuple.m_convPos[2]), "McConvPointZ/F");
//   m_tree->Branch("McCalEneSum", &(m_ntuple.m_mcCalEnergy), "McCalEneSum/F");
//   m_tree->Branch("McTkr1Ene", &(m_ntuple.m_pairEne[0]), "McTkr1Ene/F");
//   m_tree->Branch("McTkr2Ene", &(m_ntuple.m_pairEne[1]), "McTkr2Ene/F");
//   m_tree->Branch("McConvAngle", &(m_ntuple.m_convAngle), "McConvAngle/F");

//   // TKR information:
//   m_tree->Branch("TkrNumDigis", &(m_ntuple.m_nTkrDigis), "TkrNumDigis/I");
//   m_tree->Branch("TkrNumStrips", &(m_ntuple.m_nStrips), "TkrNumStrips[4][18][2]/I");
//   m_tree->Branch("tot", &(m_ntuple.m_tot), "tot[4][18][2][2]/I");
//   m_tree->Branch("TkrDepositEne", &(m_ntuple.m_depositEne), "TkrDepositEne[4][18][2]/F");
//   m_tree->Branch("TkrNumClusters", &(m_ntuple.m_nTkrClusters), "TkrNumClusters[4][18][2]/I");
//   m_tree->Branch("TkrNumTracks", &(m_ntuple.m_nTkrTracks), "TkrNumTracks/I");
//   m_tree->Branch("TkrNumVertices", &(m_ntuple.m_nTkrVertices), "TkrNumVertices/I");
//   m_tree->Branch("TkrTotalHits", &(m_ntuple.m_totalStripHits), "TkrTotalHits[4]/i");
//   m_tree->Branch("TkrTotalClusters", &(m_ntuple.m_totalClusters), "TkrTotalClusters[4]/i");
//   m_tree->Branch("Tkr1NumHits", &(m_ntuple.m_nFit[0]), "Tkr1NumHits/I");
//   m_tree->Branch("Tkr2NumHits", &(m_ntuple.m_nFit[1]), "Tkr2NumHits/I");
//   m_tree->Branch("Tkr1Chisq", &(m_ntuple.m_chi2[0]), "Tkr1Chisq/F");
//   m_tree->Branch("Tkr2Chisq", &(m_ntuple.m_chi2[1]), "Tkr2Chisq/F");
//   m_tree->Branch("Tkr1ChisqS", &(m_ntuple.m_chi2Smooth[0]), "Tkr1ChisqS/F");
//   m_tree->Branch("Tkr2ChisqS", &(m_ntuple.m_chi2Smooth[1]), "Tkr2ChisqS/F");
//   m_tree->Branch("Tkr1Rms", &(m_ntuple.m_rms[0]), "Tkr1Rms/F");
//   m_tree->Branch("Tkr2Rms", &(m_ntuple.m_rms[1]), "Tkr2Rms/F");
//   m_tree->Branch("Tkr1KalThetaMs", &(m_ntuple.m_msAngle[0]), "Tkr1KalThetaMs/F");
//   m_tree->Branch("Tkr2KalThetaMs", &(m_ntuple.m_msAngle[1]), "Tkr1Ka2ThetaMs/F");
//   m_tree->Branch("Tkr1KalEne", &(m_ntuple.m_tkrEnergy[0]), "Tkr1KalEne/F");
//   m_tree->Branch("Tkr2KalEne", &(m_ntuple.m_tkrEnergy[1]), "Tkr2KalEne/F");
//   m_tree->Branch("Tkr1EndPos", &(m_ntuple.m_tkr1EndPos), "Tkr1EndPos[3]/F");
//   m_tree->Branch("Tkr2EndPos", &(m_ntuple.m_tkr2EndPos), "Tkr2EndPos[3]/F");
//   m_tree->Branch("Tkr1EndDir", &(m_ntuple.m_tkr1EndDir), "Tkr1EndDir[3]/F");
//   m_tree->Branch("Tkr2EndDir", &(m_ntuple.m_tkr2EndDir), "Tkr2EndDir[3]/F");
//   m_tree->Branch("TkrTopTot", &(m_ntuple.m_topTot), "TkrTopTot[4]/F");
//   m_tree->Branch("Tkr1ConvTot", &(m_ntuple.m_convTot), "Tkr1ConvTot/F");


//   // Vertex information:
//   m_tree->Branch("VtxX0", &(m_ntuple.m_pos[0]), "VtxX0/F");
//   m_tree->Branch("VtxY0", &(m_ntuple.m_pos[1]), "VtxY0/F");
//   m_tree->Branch("VtxZ0", &(m_ntuple.m_pos[2]), "VtxZ0/F");
//   m_tree->Branch("VtxXDir", &(m_ntuple.m_dir[0]), "VtxXDir/F");
//   m_tree->Branch("VtxYDir", &(m_ntuple.m_dir[1]), "VtxYDir/F");
//   m_tree->Branch("VtxZDir", &(m_ntuple.m_dir[2]), "VtxZDir/F");
//   m_tree->Branch("Vtx1Energy", &(m_ntuple.m_fitTotalEnergy), "Vtx1Energy/F");
//   m_tree->Branch("Vtx1NumTkrs", &(m_ntuple.m_vtxTrks), "Vtx1NumTkrs/I");


//   // CAL information:
//   m_tree->Branch("CalEneSum", &(m_ntuple.m_calEnergy), "CalEneSum/F");
//   m_tree->Branch("CalXEcentr", &(m_ntuple.m_calPos[0]), "CalXEcentr/F");
//   m_tree->Branch("CalYEcentr", &(m_ntuple.m_calPos[1]), "CalYEcentr/F");
//   m_tree->Branch("CalZEcentr", &(m_ntuple.m_calPos[2]), "CalZEcentr/F");
//   m_tree->Branch("CalXtalEne", &(m_ntuple.m_xtalEne), "CalXtalEne[4][8][12][2]/F");
//   m_tree->Branch("CalMaxEne", &(m_ntuple.m_maxCalEnergy), "CalMaxEne/F");
//   m_tree->Branch("CalNumHit", &(m_ntuple.m_nCrystalHit), "CalNumHit[4]/I");
//   m_tree->Branch("CalXtalPos", &(m_ntuple.m_xtalPos), "CalXtalPos[4][8][12][3]/F");

//   m_tree->Branch("CalMipNum", &(m_ntuple.m_calMipNum),"CalMipNum/I");

//   m_tree->Branch("CalMip1Pos", &(m_ntuple.m_calMip1Pos),"CalMip1Pos[3]/F");
//   m_tree->Branch("CalMip1Dir", &(m_ntuple.m_calMip1Dir),"CalMip1Dir[3]/F");
//   m_tree->Branch("CalMip1Chi2", &(m_ntuple.m_calMip1Chi2),"CalMip1Chi2/F");
//   m_tree->Branch("CalMip1D2edge", &(m_ntuple.m_calMip1D2edge),"CalMip1D2edge");
//   m_tree->Branch("CalMip1ArcLen", &(m_ntuple.m_calMip1ArcLen),"CalMip1ArcLen/F");
//   m_tree->Branch("CalMip1Ecor", &(m_ntuple.m_calMip1Ecor),"CalMip1Ecor/F");
//   m_tree->Branch("CalMip1EcorRms", &(m_ntuple.m_calMip1EcorRms),"CalMip1EcorRms/F");
//   m_tree->Branch("CalMip1Erm", &(m_ntuple.m_calMip1Erm),"CalMip1Erm/F");

//   m_tree->Branch("CalMip2Pos", &(m_ntuple.m_calMip2Pos),"CalMip2Pos[3]/F");
//   m_tree->Branch("CalMip2Dir", &(m_ntuple.m_calMip2Dir),"CalMip2Dir[3]/F");
//   m_tree->Branch("CalMip2Chi2", &(m_ntuple.m_calMip2Chi2),"CalMip2Chi2/F"); 
//   m_tree->Branch("CalMip2D2edge", &(m_ntuple.m_calMip2D2edge),"CalMip2D2edge");
//   m_tree->Branch("CalMip2ArcLen", &(m_ntuple.m_calMip2ArcLen),"CalMip2ArcLen/F");
//   m_tree->Branch("CalMip2Ecor", &(m_ntuple.m_calMip2Ecor),"CalMip2Ecor/F");
//   m_tree->Branch("CalMip2EcorRms", &(m_ntuple.m_calMip2EcorRms),"CalMip2EcorRms/F");
//   m_tree->Branch("CalMip2Erm", &(m_ntuple.m_calMip2Erm),"CalMip2Erm/F");

//   // GLT information:
//   m_tree->Branch("GltWord", &(m_ntuple.m_trigger), "GltWord/i");


//   // GEM information:
//   m_tree->Branch("GemConditionsWord", &(m_ntuple.m_gemConditionsWord), "GemConditionsWord/I");
//   m_tree->Branch("GemTkrVector", &(m_ntuple.m_gemTkrVector), "GemTkrVector[4]/I");
//   m_tree->Branch("GemRoiVector", &(m_ntuple.m_gemRoiVector), "GemRoiVector[4]/I");
//   m_tree->Branch("GemCalLeVector", &(m_ntuple.m_gemCalLeVector), "GemCalLeVector[4]/I");
//   m_tree->Branch("GemCalHeVector", &(m_ntuple.m_gemCalHeVector), "GemCalHeVector[4]/I");
//   m_tree->Branch("GemCnoVector", &(m_ntuple.m_gemCnoVector), "GemCnoVector/I");
//   m_tree->Branch("GemLiveTime", &(m_ntuple.m_gemLiveTime), "GemLiveTime/i");
//   m_tree->Branch("GemTriggerTime", &(m_ntuple.m_gemTriggerTime), "GemTriggerTime/i");
//   m_tree->Branch("GemDeltaEventTime", &(m_ntuple.m_gemDeltaEventTime), "GemDeltaEventTime/i");
//   m_tree->Branch("GemOnePpsSeconds", &(m_ntuple.m_gemOnePpsSeconds), "GemOnePpsSeconds/i");
//   m_tree->Branch("GemOnePpsTime", &(m_ntuple.m_gemOnePpsTime), "GemOnePpsTime/i");
//   m_tree->Branch("GemPrescaled", &(m_ntuple.m_gemPrescaled), "GemPrescaled/i");
//   m_tree->Branch("GemDiscarded", &(m_ntuple.m_gemDiscarded), "GemDiscarded/i");
//   m_tree->Branch("GemCondArrivalTimeWord",&(m_ntuple.m_gemCondArrivalTimeWord), "GemCondArrivalTimeWord/i");
//   m_tree->Branch("GemCondArrivalTimeExt",&(m_ntuple.m_gemCondArrivalTimeExt), "GemCondArrivalTimeExt/i");
//   m_tree->Branch("GemCondArrivalTimeCno",&(m_ntuple.m_gemCondArrivalTimeCno), "GemCondArrivalTimeCno/i");
//   m_tree->Branch("GemCondArrivalTimeCalLe",&(m_ntuple.m_gemCondArrivalTimeCalLe), "GemCondArrivalTimeCalLe/i");
//   m_tree->Branch("GemCondArrivalTimeCalHe",&(m_ntuple.m_gemCondArrivalTimeCalHe), "GemCondArrivalTimeCalHe/i");
//   m_tree->Branch("GemCondArrivalTimeTkr",&(m_ntuple.m_gemCondArrivalTimeTkr), "GemCondArrivalTimeTkr/i");
//   m_tree->Branch("GemCondArrivalTimeRoi",&(m_ntuple.m_gemCondArrivalTimeRoi), "GemCondArrivalTimeRoi/i");
//   m_tree->Branch("GemDeltaWindowOpenTime",&(m_ntuple.m_gemDeltaWindowOpenTime), "GemDeltaWindowOpenTime/i");
//   m_tree->Branch("GemDeadZone",&(m_ntuple.m_gemDeadZone), "GemDeadZone/i");
//   m_tree->Branch("GemAcdTilesXzp", &(m_ntuple.m_gemAcdTilesXzp), "GemAcdTilesXzp/i");
//   m_tree->Branch("GemAcdTilesXzm", &(m_ntuple.m_gemAcdTilesXzm), "GemAcdTilesXzm/i");
//   m_tree->Branch("GemAcdTilesYzp", &(m_ntuple.m_gemAcdTilesYzp), "GemAcdTilesYzp/i");
//   m_tree->Branch("GemAcdTilesYzm", &(m_ntuple.m_gemAcdTilesYzm), "GemAcdTilesYzm/i");
//   m_tree->Branch("GemAcdTilesXy", &(m_ntuple.m_gemAcdTilesXy), "GemAcdTilesXy/i");
//   m_tree->Branch("GemAcdTilesRbn", &(m_ntuple.m_gemAcdTilesRbn), "GemAcdTilesRbn/i");
//   m_tree->Branch("GemAcdTilesNa", &(m_ntuple.m_gemAcdTilesNa), "GemAcdTilesNa/i");


//   // TME diagnostic information:
//   m_tree->Branch("DigiTriRowBits",&(m_ntuple.m_digiTriRowBits),"DigiTriRowBits[4]/i");
//   m_tree->Branch("TrgReqTriRowBits",&(m_ntuple.m_trgReqTriRowBits),"TrgReqTriRowBits[4]/i");
//   m_tree->Branch("TkrReq", &(m_ntuple.m_tkrReq), "TkrReq[4][18][2][2]/i");
//   m_tree->Branch("TkrTp", &(m_ntuple.m_tpTkr), "TkrTp[4][8]/i");
//   m_tree->Branch("CalReq", &(m_ntuple.m_calReq), "CalReq[4][8][2]/i");
//   m_tree->Branch("CalTp", &(m_ntuple.m_tpCal), "CalTp[4][8]/i");


//   // Contribution lenghts:
//   m_tree->Branch("DiagLength", &(m_ntuple.m_diagLength), "DiagLength[4]/i");
//   m_tree->Branch("TemLength", &(m_ntuple.m_temLength), "TemLength[4]/i");
//   m_tree->Branch("GemLength", &(m_ntuple.m_gemLength), "GemLength/i");
//   m_tree->Branch("OswLength", &(m_ntuple.m_oswLength), "OswLength/i");
//   m_tree->Branch("AemLength", &(m_ntuple.m_aemLength), "AemLength/i");
//   m_tree->Branch("ErrLength", &(m_ntuple.m_errLength), "ErrLength[4]/i");


//   // ACD digi:
//   m_tree->Branch("AcdNumDigis", &(m_ntuple.m_acdNumDigis), "AcdNumDigis/I");
//   m_tree->Branch("AcdPha", &(m_ntuple.m_acdPha), "AcdPha[8][2]/I");
//   m_tree->Branch("AcdHitMap", &(m_ntuple.m_acdHitMap), "AcdHitMap[8][2]/I");
//   m_tree->Branch("AcdRange", &(m_ntuple.m_acdRange), "AcdRange[8][2]/I");
//   m_tree->Branch("AcdOddParityError", &(m_ntuple.m_acdOddParityError), "AcdOddParityError[8][2]/I");
//   m_tree->Branch("AcdHeaderParityError", &(m_ntuple.m_acdHeaderParityError), "AcdHeaderParityError[8][2]/I");
//   m_tree->Branch("AcdLowDisc", &(m_ntuple.m_acdLowDisc), "AcdLowDisc[8][2]/I");
//   m_tree->Branch("AcdTileNumber", &(m_ntuple.m_acdTileNumber), "AcdTileNumber[8]/I");
//   m_tree->Branch("AcdMCEnergy", &(m_ntuple.m_acdMCEnergy), "AcdMCEnergy[8]/F");

//   // ACD recon:
//   m_tree->Branch("AcdEnergy", &(m_ntuple.m_acdEnergy),"AcdEnergy/F");
//   m_tree->Branch("AcdTileCount", &(m_ntuple.m_acdTileCount),"AcdTileCount/I");

//   // ACD MIPs:
//   m_tree->Branch("AcdMips", &(m_ntuple.m_acdMips),"AcdMips[8][2]/F");
// }
