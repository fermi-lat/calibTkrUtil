#ifndef RootAnalyzer_Class
#define RootAnalyzer_Class

#include <string>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
//#include "mcRootData/McEvent.h"
#include "reconRootData/ReconEvent.h"
#include "digiRootData/DigiEvent.h"
//#include "idents/CalXtalId.h"
//#include "NtupleDef.h"
#include "calibTkrUtil/TkrHits.h"
#include "calibTkrUtil/TkrNoiseOcc.h"

class RootAnalyzer {
 public:

  RootAnalyzer();
  ~RootAnalyzer();

  /// parse option file to the RootAnalyzer
  void parseOptionFile(const char* f);
  void setInputRootFiles( const string&, const string&, const string& );
  void setOutputRootFile( const char* );

  void analyzeData( int );

  /// produce output ntuple and histogram files based on analysis results
  void produceOutputFile();

  //enum {g_nTower=16, g_nTkrLayer=18, g_nView=2, g_nCalLayer=8, g_nCol=12, 
  //        g_nTot=2, g_nStripsPerLayer=1536, g_nFace=2, g_nTP=8, g_nFEC=24, 
  //        g_nCno=12};

 private:

  //  void analyzeMcTree();

  //  void analyzeReconTree();

  //  void analyzeDigiTree();

  /// determine whether input string is empty or a comment  
  /// a comment  is defined as a string with "//" in front
  bool isEmptyOrCommentStr(const std::string& s);

  /// determine whether a file is a root file
  bool isRootFile(const std::string& f);

  /// make a root TChain based on file names in a line, the file names should
  /// be separated by one or more empty spaces
  void  makeTChain(const std::string& line, TChain* chain);

  /// input a line of string containing sub strings separated by empty spaces,
  /// output first such sub string
  void parseLine(const std::string& line, std::string& str);

  // parse diagnostic data (trigger primitive)
  //  void parseDiagnosticData();

  // extract tower number (from 0 to 16) layer number (from 0 to 17) and view 
  // number (X:0, Y:1) from id
  // return false if the volume identified by id is not a tracker si wafer
  // or tower or layer or view number fall outside the boundary.
  // Convertion is hardwired so it may not work if identifier in xml files
  // changed
/*   bool extractTowerLayerView(const VolumeIdentifier& id, int& iTower,  */
/*                              int& iLayer, int& iView) const; */

  // fill stripHits and map histogrammes
  //void fillStripHits(const TkrDigi* tkrDigi);

  //void fillOutputTree();

  // calculate a couple of tot values based on event topology
  //  void analyzeTot();

  // create branches for each ntuple variable
  //  void createBranches();

  TFile* m_outputFile;
  /*  TTree* m_tree;
      TBranch* m_branch; */
  //  NtupleDef m_ntuple;
  TChain* m_mcChain;
  TBranch* m_mcBranch;
  //  McEvent* m_mcEvent;

  TChain* m_reconChain;
  TBranch* m_reconBranch;
  ReconEvent* m_reconEvent;

  TChain* m_digiChain;
  TBranch* m_digiBranch;
  DigiEvent* m_digiEvent;

  /*  TFile* m_histFile;
  TH1F* m_stripHits[g_nTower][g_nTkrLayer][g_nView];
  TH2F* m_stripMap[g_nTower][g_nTkrLayer][g_nView];
  */

  TkrHits* m_tkrHits;

  TkrNoiseOcc* m_tkrNoiseOcc;
  TDirectory* m_tkrNoiseOcc_dir;
};

#endif
