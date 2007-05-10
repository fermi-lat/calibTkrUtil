#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TLine.h"
#include "TText.h"
#include "TLegend.h"
#include <iostream>
#include <string.h>
#include <time.h>

class TkrNoiseRep {

 public:
  //TkrNoiseRep(const char* reportDir, const char* prefix, const char* version);
  TkrNoiseRep();
  ~TkrNoiseRep();

  void openSvacFile(const char *filename);
  //void openNewHistFile(const char *filename);
  void makeAncRootFile();
  void setReportDirName(const char *reportDirName);
  void setPrefix(const char *prefix);

  void findNoisyLayer();

  void drawAncRootHist(const char *hname);

  void drawStripOccAve();
  void drawStripOccMax();
  void drawLayerOccAve();
  void drawLayerOccMax();
  void drawLargeMultiRatio();

  void drawTowerAveStripOccHist();
  void drawStripOcc_TowerAve(int tower);
  void drawStripOccGr_TowerAve(int tower);
  void drawStripOccGr_LatAve();

  void drawStripOcc_perLayer(int tower, int layer);
  void drawLayerOcc_perLayer(int tower, int layer);
  void drawStripOccGr_perLayer(int tower, int layer);
  void drawLayerOccGr_perLayer(int tower, int layer);

  void drawStripHist(int tower, int layer);
  void drawMultiHist(int tower, int layer);
  void drawTotHist(int tower, int layer);

  void epstogif(const char *figFileDir, const char *epsFileName, const char *gifFileName, int xsize, int ysize);

  void getTwrOccList();
  void getTimeStr();

  void writeSummaryXML();
  void generateSummary();
  void generateSummaryPage();
  void generateNoisyLayerReport();
  void generateLayerReport(int tower, int layer);


 private:

  TFile *m_svacFile;
  TFile *m_ancRootFile;
  const char  *m_reportDirName;
  const char  *m_prefix;

  int    m_test_status;
  double m_critical_strip_occ;
  double m_critical_layer_occ;
  int    m_critical_multi;
  double m_critical_multi_ratio;
  double m_latAveOccMean, m_latAveOccMax;
  double m_aveTwrOccList[16], m_maxTwrOccList[16];

  int    m_failStripOcc[16][18][2];
  int    m_failLayerOcc[16][18][2];
  int    m_failMultiRatio[16][18][2];

  double m_start_time;
  double m_end_time;
  float  m_minExpos;
  char   m_startTimeStr[80];
  char   m_durationStr[80];
  char   m_timeAxTitle[80];
};
