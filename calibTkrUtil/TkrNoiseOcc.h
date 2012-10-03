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

#define TOT_MAX 300
#define MUL_MAX 150
#define TOT_INI_VAL 295


class paramTimeDep{
 public:
  paramTimeDep( UInt_t id, UInt_t startT, UInt_t duration, UInt_t nx );
  ~paramTimeDep(){;};
  void clear();
  void fill( int tower, int bilayer, int view, Double_t timeStamp,
             float nExp, float nOcc );
 private:
  UInt_t m_id, m_startTime, m_duration, m_nx;
  float m_binSize;
  float *m_Exposure[g_nTower][g_nTkrLayer][g_nView];
  float *m_StripOcc[g_nTower][g_nTkrLayer][g_nView];
  float *m_LayerOcc[g_nTower][g_nTkrLayer][g_nView];

 public:
  inline UInt_t id(){ return m_id; };
  inline UInt_t startTime(){ return m_startTime; };
  inline float* getExposure( UInt_t tower, UInt_t bilayer, UInt_t view){
    return m_Exposure[tower][bilayer][view]; };
  inline float* getStripOcc( UInt_t tower, UInt_t bilayer, UInt_t view){
    return m_StripOcc[tower][bilayer][view]; };
  inline float* getLayerOcc( UInt_t tower, UInt_t bilayer, UInt_t view){
    return m_LayerOcc[tower][bilayer][view]; };
};



class TkrNoiseOcc {
 public:
  TkrNoiseOcc();
  ~TkrNoiseOcc();  

  //void initAnalysis(int nEvent, int evt_interval, int coincidence_cut,
  //                    int multi_ld, int multi_hd, int periodic_trig);
  void initAnalysis(UInt_t duration=120, UInt_t nx=60);
  void setCoincidenceCut(int coincidence_cut);
  void setMultiRange(int multi_ld, int multi_hd);
  void setPeriodicTrigCut(int periodic_trig);
  void setTrigCut(int trig_cut);

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
  int    m_coincidence_cut, m_multi_ld, m_multi_hd, m_periodic_trig, m_trig_cut;
  float  m_crit_strip_rate, m_crit_layer_rate;
  /// data parameter
  //int    m_nEvent, m_evt_interval, m_event_counter;
  UInt_t    m_nx, m_duration;
  Double_t m_evtTime;

  /// data array definition
  float  *vTkrHitMap[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrWHitMap[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrNoiseMul[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrNoiseTot0[g_nTower][g_nTkrLayer][g_nView];
  float  *vTkrNoiseTot1[g_nTower][g_nTkrLayer][g_nView];

  UInt_t m_iP;
  std::vector<paramTimeDep> vParamTimeDep;

};
