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

#define NUMTOWER 16
#define NUMLAYER 18
#define NUMVIEW  2
#define NUMSTRIP 1536

#define TOT_MAX 300
#define MUL_MAX 150
#define TOT_INI_VAL 295

class TkrNoiseOcc {
 public:
  TkrNoiseOcc();
  ~TkrNoiseOcc();  

  //void initPar(int coincidence_cut, int multi_ld, int multi_hd);
  //void initArrDat( int nEvent, int evt_interval );
  //void initPar(int coincidence_cut, int multi_ld, int multi_hd);
  void initAnalysis(int nEvent, int evt_interval, int coincidence_cut, int multi_ld, int multi_hd);
  void clearAnalysis();

  void setDigiEvtPtr(DigiEvent *digiEvt);
  void anaDigiEvt();

  void openHistFile(char* histFileName);
  void closeHistFile();
  void writeAnaToHis(TDirectory* tkrNoiseOcc_dir);
  void saveAnaToHis(char* histFileName);
  
  void setCritStripRate(float crit_strip_rate);
  void setCritLayerRate(float crit_layer_rate);
 private:

  DigiEvent *m_digiEvt;
  
  /// output ROOT file
  TFile *m_histFile;
  
  /// analysis parameter
  int m_coincidence_cut, m_multi_ld, m_multi_hd;
  float m_crit_strip_rate, m_crit_layer_rate;
  /// data parameter
  int m_nEvent, m_evt_interval, m_nx, m_event_counter;

  /// data array definition
  float  *vTkrExposure[NUMTOWER][NUMLAYER];
  float  *vTkrStripOcc[NUMTOWER][NUMLAYER][NUMVIEW];
  float  *vTkrLayerOcc[NUMTOWER][NUMLAYER][NUMVIEW];
  float  *vTkrHitMap[NUMTOWER][NUMLAYER][NUMVIEW];
  float  *vTkrNoiseMul[NUMTOWER][NUMLAYER][NUMVIEW];
  float  *vTkrNoiseTot0[NUMTOWER][NUMLAYER][NUMVIEW];
  float  *vTkrNoiseTot1[NUMTOWER][NUMLAYER][NUMVIEW];

};


