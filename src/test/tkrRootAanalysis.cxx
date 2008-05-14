#include <iostream>
#include <fstream>
#include <string>
#include "RootTestAnalyzer.h"

// If there is no argument, use AnalyzerOption.dat as input, otherwise use 
// first argument as the input

int main(int argc, char** argv)
{
  std::string optionFileName;
  int numEvents = 10000;
  if(argc == 2 || argc==3) numEvents = atoi( argv[1] );

  RootAnalyzer analyzer;
  if(argc == 3)
    analyzer.parseOptionFile( argv[2] );
  else if(argc == 4) {
    numEvents = -1; // analyze all events
    analyzer.setInputRootFiles( "none", argv[1], argv[2] );
    analyzer.setOutputRootFile( argv[3] );
  }
  else {
    optionFileName = "../src/test/AnalyzerOption.dat";
    analyzer.parseOptionFile(optionFileName.c_str());
  }

  analyzer.analyzeData( numEvents );

  analyzer.produceOutputFile();

}

