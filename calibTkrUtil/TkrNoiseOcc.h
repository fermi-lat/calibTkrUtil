#include "TNtuple.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCollection.h"  // Declares TIter
#include "TObjArray.h"
#include "digiRootData/DigiEvent.h"
//#include "reconRootData/ReconEvent.h"
//#include "mcRootData/McEvent.h"
#include <iostream>

#include "GeoConstants.h"

/*
#define NUMTOWER 16
#define NUMLAYER 18
#define NUMVIEW  2
#define NUMSTRIP 1536
*/

#define TOT_MAX 300
#define MUL_MAX 150
#define TOT_INI_VAL 295

class TkrNoiseOcc {
 public:
  TkrNoiseOcc();
  ~TkrNoiseOcc();  

  //void initAnalysis(int nEvent, int evt_interval, int coincidence_cut,
  //		    int multi_ld, int multi_hd, int periodic_trig);
  void initAnalysis(int nEvent, int evt_interval);
  void setCoincidenceCut(int coincidence_cut);
  void setMultiRange(int multi_ld, int multi_hd);
  void setPeriodicTrigCut(int periodic_trig);
  void clearAnalysis();

  void setDigiEvtPtr(DigiEvent *digiEvt);
  void anaDigiEvt();

  void openHistFile(char* histFileName);
  void closeHistFile();
  void writeAnaToHis(TDirectory* tkrNoiseOcc_dir);
  void saveAnaToHis(char* histFileName);
  
  void setCritStripRate(float crit_strip_rate);
  void setCritLayerRate(float crit_layer_rate);
  void setTrigCut(int trig_cut);

 private:

  DigiEvent *m_digiEvt;
  
  /// output ROOT file
  TFile *m_histFile;
  
  /// analysis parameter
  int m_coincidence_cut, m_multi_ld, m_multi_hd, m_periodic_trig, m_trig_cut;
  float m_crit_strip_rate, m_crit_layer_rate;
  /// data parameter
  int m_nEvent, m_evt_interval, m_nx, m_event_counter;

  /// data array definition
  float  *vTkrExposure[g_nTower][g_nTkrLayer];
  float  *vTkrStripOcc[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrLayerOcc[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrHitMap[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrNoiseMul[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrNoiseTot0[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrNoiseTot1[g_nTower][g_nTkrLayer][g_nView];

};


