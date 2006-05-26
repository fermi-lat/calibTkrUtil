#include <iostream>
#include <fstream>
#include <string>
#include "RootTestAnalyzer.h"

// If there is no argument, use AnalyzerOption.dat as input, otherwise use 
// first argument as the input

int main(int argc, char** argv)
{
  std::string optionFileName;
  if(argc == 1) {
    optionFileName = "../src/test/AnalyzerOption.dat";
  }
  else {
    optionFileName = argv[1];
  }

  RootAnalyzer analyzer;

  analyzer.parseOptionFile(optionFileName.c_str());

  analyzer.analyzeData();

  analyzer.produceOutputFile();

}

