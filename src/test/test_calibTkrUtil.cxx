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
  if(argc > 1) numEvents = atoi( argv[1] );
  if(argc > 2) {
    optionFileName = argv[2];
  }
  else {
    optionFileName = "../src/test/AnalyzerOption.dat";
  }

  RootAnalyzer analyzer;

  analyzer.parseOptionFile(optionFileName.c_str());

  analyzer.analyzeData( numEvents );

  analyzer.produceOutputFile();

}

