#include <string>
#include <fstream>
#include "TkrNoiseRep.h"

int main(int argc, char** argv)
{
  
  std::string svacRootFile, reportDir, prefix;
  std::string optionFileName;

  if (argc==4) {
    svacRootFile = argv[1];
    reportDir    = argv[2];
    prefix       = argv[3];

  } 
  else {
    if (argc==1) {
      optionFileName = "../src/ReportOption.txt";
    } else {
      optionFileName = argv[1];
    }
    std::ifstream optionFile(optionFileName.c_str());

    optionFile >> svacRootFile >> reportDir >> prefix;

  } 

  std::cout << "svacRootFile = " << svacRootFile.c_str() << std::endl;
  std::cout << "reportDir    = " << reportDir.c_str() << std::endl;
  std::cout << "prefix       = " << prefix.c_str() << std::endl;
  

  TkrNoiseRep *tkrNoiseRep = new TkrNoiseRep();

  tkrNoiseRep->openSvacFile(svacRootFile.c_str());
  tkrNoiseRep->setReportDirName(reportDir.c_str());
  tkrNoiseRep->setPrefix(prefix.c_str());

  tkrNoiseRep->generateSummary();
  tkrNoiseRep->writeSummaryXML();
  tkrNoiseRep->generateNoisyLayerReport();


}
