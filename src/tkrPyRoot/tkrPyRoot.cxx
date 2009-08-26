
//-----------------------------------------------------------------------
//
//      Convoluted Landau and Gaussian Fitting Funcion
//         (using ROOT's Landau and Gauss functions)
//
//  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
//  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
//  Markus Friedl (Markus.Friedl@cern.ch)
//  Modified by Hiro Tajima (tajima@stanford.edu) for speed up.
//
//-----------------------------------------------------------------------

#ifndef SCons
#  include "tkrPyRoot.h"
#else
#  include "tkrPyRoot/tkrPyRoot.h"
#endif

TF1* defLangau( char* name, Double_t xmin = 0, Double_t xmax = 0 ){
  return new TF1( name, langaufun, xmin, xmax, 4 );
}

TF1* defLangau2( char* name, Double_t xmin = 0, Double_t xmax = 0 ){
  return new TF1( name, langau2fun, xmin, xmax, 6 );
}

double add( double x, double y ){ return x*2 + y; }

bool test(const std::string &x, int n, std::vector<double>& z){

  std::string y(x);
  y += "+test";

  for( int i=0;i<n;i++){
    z.push_back( i*2.0 + 1.0 );
  }
  return true;
}


double langaufun(double *x, double *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 40.0;   // number of convolution steps
  Double_t sc =  2.0;   // convolution extends to +-sc Gaussian sigmas
  Double_t np2 = 20.0;  // number of convolution steps
  Double_t sc2 =  4.0;  // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland, fgauss;
  Double_t sum = 0.0;
  Double_t xm, dx;
  Double_t step;
  Double_t sigma = par[3];
  Double_t width = par[0];
  if( sigma<=0.0 || width<=0.0 ) return 0.0;

  // MP shift correction
  mpc = par[1] - mpshift * width;
  if( mpc < 0.0 ) return 0.0;
  
  // Convolution integral parameters
  xm = x[0];
  Double_t dxmax = sc * sigma;
  step = dxmax * 2 / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(dx=0.5*step; dx<=dxmax; dx+=step) {
    fgauss = TMath::Gaus(0.0,dx,sigma);
    
    xx = xm + dx;
    fland = TMath::Landau( xx, mpc, width );
    sum += fland * fgauss;
    
    xx = xm - dx;
    fland = TMath::Landau( xx, mpc, width );
    sum += fland * fgauss;
  }
  sum *= step/width;
  
  Double_t dxmax2 = sc2 * sigma;
  step = (dxmax2-dxmax) * 2 / np2;
  Double_t sum2 = 0.0;
  
  // Convolution integral of Landau and Gaussian by sum
  for(dx=dxmax+0.5*step; dx<=dxmax2; dx+=step) {
    fgauss = TMath::Gaus(0.0,dx,sigma);
    
    xx = xm + dx;
    fland = TMath::Landau( xx, mpc, width );
    sum2 += fland * fgauss;
    
    xx = xm - dx;
    fland = TMath::Landau( xx, mpc, width );
    sum2 += fland * fgauss;
  }
  sum2 *= step/width;
  
  return (par[2] * (sum+sum2) * invsq2pi / sigma);
}

//
// Two Gaussian convolved Lanbdau function
// it also employs variable integration step.
//
Double_t langau2fun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //par[4]=Width2/Width of convoluted Gaussian function
   //par[5]=fraction of narrow Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation),
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  // variable interation step
  // integration step increases as points move to outside.
  Double_t np = 20.0;   // number of convolution steps in a region
  Double_t sc =  1.0;   // convolution extends to +-sc Gaussian sigmas
  Double_t sc2 =  4.0;  // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland, fgauss;
  Double_t sum = 0.0;
  Double_t xm, dx;
  Double_t step;
  Double_t sigma = par[3];
  Double_t sigma2 = par[4]*sigma;
  Double_t frac = par[5];
  Double_t width = par[0];
  if( sigma<=0.0 || width<=0.0 ) return 0.0;

  // MP shift correction
  mpc = par[1] - mpshift * width;
  if( mpc < 0.0 ) return 0.0;
  
  // Convolution integral parameters
  xm = x[0];
  Double_t dxmax = sc * sigma;
  step = dxmax * 2 / np;
  Double_t dxmax2 = sc2 * sigma2;
  
  // Convolution integral of Landau and Gaussian by sum
  for(dx=0.5*step; dx<=dxmax2; dx+=step) {
    fgauss = step * ( frac * TMath::Gaus(0.0,dx,sigma) / sigma
                      + (1-frac) * TMath::Gaus(0.0,dx,sigma2) / sigma2 );
    
    xx = xm + dx;
    fland = TMath::Landau( xx, mpc, width );
    sum += fland * fgauss;
    
    xx = xm - dx;
    fland = TMath::Landau( xx, mpc, width );
    sum += fland * fgauss;

    if( dx >= dxmax ){ // getting out of the curent region, change step
      dxmax *= 2;
      step *= 2;
    }
  }
  sum /= width;
  
  return ( par[2] * sum * invsq2pi );
}
