
#ifndef tkrPyRoot_h
#define tkrPyRoot_h

#include <string>
#include <vector>

#include "TMath.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"


double langaufun(double *, double *);
double langau2fun(double *, double *);
double add( double, double );
TF1* defLangau( char*, Double_t, Double_t);
TF1* defLangau2( char*, Double_t, Double_t);

#endif
