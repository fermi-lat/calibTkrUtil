
#include <cmath>
#include <ctime>
#include <cassert>

#include "facilities/Util.h"
#include "commonRootData/idents/TowerId.h"
#include "commonRootData/idents/CalXtalId.h"

#include "calibTkrUtil/TkrHits.h"
#include "src/tkrPyRoot/tkrPyRoot.h"

using std::string;
using std::cout;
using std::endl;

//#define PRINT_DEBUG 1



bool Cluster::addStrip( int strip ){

  if( strip == lastStrip+1 ){ // assume strip number is sorted.
    lastStrip = strip;
    return true;
  }
  else return false;

}

//
// layerId class implementation
//
layerId::layerId( int lyr, int vw, int twr ){ 
  setLayer( lyr, vw, twr ); }
layerId::layerId( int tr, std::string wh, int twr ){ 
  setTray( tr, wh, twr ); }
layerId::layerId( int unp ){ 
  setUniPlane( unp ); }
layerId::layerId( std::string tspt, int twr ){ 
  setSPT( tspt, twr ); }

void layerId::setLayer( int lyr, int vw, int twr ){
  layer = lyr; view = vw; tower = twr;
  layerToTray();
  trayToUniPlane();
  layerToSPT();
}

void layerId::setUniPlane( int unp, int twr ){
  uniPlane = unp; tower = twr;
  uniPlaneToTray();
  trayToLayer();
  layerToSPT();
}

void layerId::setTray( int tr, std::string wh, int twr ){
  tray = tr; which=wh; tower=twr;
  trayToUniPlane();
  trayToLayer();
  layerToSPT();
}

void layerId::setSPT( std::string tspt, int twr ){
  spt = tspt; tower = twr;
  sptToLayer();
  layerToTray();
  trayToUniPlane();
}

void layerId::trayToUniPlane(){
  uniPlane = tray * 2;
  if( which == "bot" ) uniPlane--;
}

void layerId::trayToLayer(){
  view = (tray+1) % 2;
  layer = tray;
  if( which == "bot" ) layer--;
}

void layerId::layerToTray(){
  if(view==0){
    tray = 2 * ( layer/2 ) + 1;
    if(layer%2==0) which = "bot";
    else which = "top";
  }
  else{
    tray = 2 * ( (layer+1)/2 );
    if(layer%2==0) which = "top";
    else which = "bot";
  }
}

void layerId::uniPlaneToTray(){
  tray = (uniPlane+1) / 2;
  if( uniPlane%2 == 0 ) which ="top";
  else which = "bot";
}

void layerId::layerToSPT(){
  char cspt[] = "+x0";
  if( view==0 )
    if( layer%2==0 ) sprintf( cspt, "+y%d", layer/2 );
    else sprintf( cspt, "-y%d", layer/2 );
  else
    if( layer%2==0 ) sprintf( cspt, "+x%d", layer/2 );
    else sprintf( cspt, "-x%d", layer/2 );
  spt = cspt;
}

void layerId::sptToLayer(){
  // carefull this is different from normal definition.
  if( spt.substr(1,1) == 'y' ) view = 0; 
  else view = 1;
  if( spt.substr(0,1) == "+" ) layer = 2*atoi(spt.substr(2,1).c_str());
  else layer = 2*atoi(spt.substr(2,1).c_str()) + 1;
}

std::string layerId::getLayerName(){
  char lname[] = "X17";
  std::string layerName;
  if( view==0 ) sprintf(lname, "X%d", layer );
  else sprintf(lname, "Y%d", layer );
  layerName = lname;
  return layerName;
}

//
// towerVar class implementation
//
//
// towerVar class implementation
//
towerVar::towerVar( int twr, bool badStrips ){
  towerId = twr;
  hwserial = "None";
  runid = "-1";
  bsVar.clear();
  tcVar.clear();

  for( int unp=0; unp!=g_nUniPlane; unp++){
    badStripVar bsv;
    totCalibVar tcv;
    bsv.hLayer = 0;
    bsv.tLayer = 0;
    for( int strip=0; strip!=g_nStrip; strip++){
      rHits[unp][strip] = 0;
      dHits[unp][strip] = 0;
      if( badStrips ){
        bsv.eHits[strip] = 0;
        bsv.tHits[strip] = 0;
        bsv.lHits[strip] = 0;
        for(int iWafer = 0; iWafer != g_nWafer; ++iWafer)
          for( int tDiv = 0; tDiv != g_nTime; tDiv++)
            bsv.nHits[strip][iWafer][tDiv] = 0;
      }
      tcv.totQuad[strip] = -1.0;
      tcv.totGain[strip] = -1.0;
      tcv.totThreshold[strip] = -1.0;
    }
    for( int idiv=0; idiv!=g_nDiv; idiv++){
      tcv.chargeScale[idiv] = 0.0;
      for( int ibin=0; ibin!=nTotHistBin; ibin++)
        tcv.chargeDist[idiv][ibin] = 0;
    }
    if( badStrips ) bsVar.push_back( bsv );
    tcVar.push_back( tcv );
  }

  char name[] = "numHitGTRCT00";
  sprintf(name,"numHitGTRCT%d", towerId);
  numHitGTRC = new histogram(name, name, 65, 0, 65);
  sprintf(name,"resXT%d", towerId);
  resX = new histogram( name, name, 100, -3.0, 3.0 );
  sprintf(name,"resYT%d", towerId);
  resY = new histogram( name, name, 100, -3.0, 3.0 );
#ifndef ROOT_PROFILE
  resProf = new profile( g_nUniPlane, 0, g_nUniPlane );
#else
  sprintf(name,"resProfT%d", towerId);
  resProf = new TProfile( name, name, g_nUniPlane, 0, g_nUniPlane );
#endif
#ifdef PRINT_DEBUG
  std::cout << "towerVar constructer " << twr << std::endl;
#endif
}


void towerVar::saveHists( bool saveTimeOcc ){ 

  TH1F* hist, *rhist, *dhist, *ehist, *thist, *lhist;
  char name[] = "roccT00X17w3t4";
  char dname[] = "T00";

  sprintf(dname,"T%d", towerId);
  gDirectory->mkdir( dname, dname );
  gDirectory->cd( dname );

  for( int unp=0; unp!=g_nUniPlane; unp++){
    layerId lid( unp );
    std::string lname = lid.getLayerName();
    sprintf(name,"roccT%d%s", towerId, lname.c_str());
    rhist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    sprintf(name,"doccT%d%s", towerId, lname.c_str());
    dhist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    for( int strip=0; strip!=g_nStrip; strip++){
      rhist->Fill( strip+0.1, rHits[unp][strip] );
      dhist->Fill( strip+0.1, dHits[unp][strip] );
    }
    rhist->Write(0, TObject::kOverwrite);
    dhist->Write(0, TObject::kOverwrite);
    delete rhist;
    delete dhist;
  }

  sprintf(name,"lhitT%d", towerId);
  dhist = new TH1F(name, name, g_nUniPlane, 0, g_nUniPlane);
  sprintf(name,"ltrkT%d", towerId);
  rhist = new TH1F(name, name, g_nUniPlane, 0, g_nUniPlane);
  for( UInt_t unp=0; unp!=bsVar.size(); unp++){
    dhist->Fill( unp+0.1, bsVar[unp].hLayer );
    rhist->Fill( unp+0.1, bsVar[unp].tLayer );
    layerId lid( unp );
    std::string lname = lid.getLayerName();
    sprintf(name,"eoccT%d%s", towerId, lname.c_str());
    ehist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    sprintf(name,"toccT%d%s", towerId, lname.c_str());
    thist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    sprintf(name,"loccT%d%s", towerId, lname.c_str());
    lhist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
    for( int strip=0; strip!=g_nStrip; strip++){
      ehist->Fill( strip+0.1, bsVar[unp].eHits[strip] );
      thist->Fill( strip+0.1, bsVar[unp].tHits[strip] );
      lhist->Fill( strip+0.1, bsVar[unp].lHits[strip] );
    }
    ehist->Write(0, TObject::kOverwrite);
    thist->Write(0, TObject::kOverwrite);
    lhist->Write(0, TObject::kOverwrite);
    delete ehist;
    delete thist;
    delete lhist;
    if( saveTimeOcc ){
      for(int iWafer = 0; iWafer != g_nWafer; ++iWafer)
        for( int tDiv = 0; tDiv != g_nTime; tDiv++){
          sprintf(name,"occT%d%sw%dt%d", towerId, lname.c_str(), iWafer, tDiv);
          hist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
          for( int strip=0; strip!=g_nStrip; strip++)
            hist->Fill( strip+0.1, bsVar[unp].nHits[strip][iWafer][tDiv] );
          hist->Write(0, TObject::kOverwrite);
          delete hist;
        }
    }
    else{
      for(int iWafer = 0; iWafer != g_nWafer; ++iWafer){
        sprintf(name,"occT%d%sw%d", towerId, lname.c_str(), iWafer);
        hist = new TH1F(name, name, g_nStrip, 0, g_nStrip);
        for( int strip=0; strip!=g_nStrip; strip++)
          hist->Fill( strip+0.1, bsVar[unp].nHits[strip][iWafer][0] );
        hist->Write(0, TObject::kOverwrite); 
        delete hist;
     }
    }
  }
  gDirectory->cd( "../Towers" );
  sprintf(name,"leffT%d", towerId);
  hist = new TH1F(name, name, g_nUniPlane, 0, g_nUniPlane);
  hist->Divide( dhist, rhist );
  for( int unp=0; unp<g_nUniPlane; unp++){
    // calculate binomial error
    float eff = hist->GetBinContent( unp+1 );
    float entry = rhist->GetBinContent( unp+1 );
    float error = 1.0E-5;
    if( entry > 0 ){
      error = eff*(1-eff)/entry;
      if( error > 0.0 ) error = sqrt( error );
      else error = 1.0E-5;
    }
    hist->SetBinError( unp+1, error );
  }
  dhist->Write(0, TObject::kOverwrite);
  rhist->Write(0, TObject::kOverwrite);
  hist->Write(0, TObject::kOverwrite);
  delete dhist;
  delete rhist;
  delete hist;

  numHitGTRC->save();
  resX->save();
  resY->save();
#ifdef ROOT_PROFILE
  resProf->Write(0, TObject::kOverwrite);
#else
  sprintf(name,"resHistT%d", towerId);
  hist = new TH1F( name, name, g_nUniPlane, 0, g_nUniPlane );
  for( int unp=0; unp!=g_nUniPlane; unp++){
    hist->SetBinContent( unp+1, resProf->getMean( unp ) );
    hist->SetBinError( unp+1, resProf->getRMS( unp ) );
  }
  hist->Write(0, TObject::kOverwrite);
#endif
  gDirectory->cd( ".." );
}

profile::profile( int nbin, float min, float max ){

  updated = false;
  numBin = nbin;
  binMin = min;
  binMax = max;
  binSize = (binMax-binMin) / numBin;

  for( int i=0; i!=numBin; i++){
    num.push_back( 0 );
    sum.push_back( 0.0 );
    sumsq.push_back( 0.0 );
    mean.push_back( 0.0 );
    rms.push_back( 0.0 );
  }
}

void profile::Fill( float x, float val ){

  if( x < binMin || x >= binMax ) return;
  int bin = (int) (( x - binMin ) / binSize);
  if( bin >= numBin ) return;

  updated = false;
  num[bin] += 1;
  sum[bin] += val;
  sumsq[bin] += val*val;

  return;

}

void profile::calculate(){

  for( UInt_t bin=0; bin<num.size(); bin++){
    if( num[bin] == 0 ) continue;
    mean[bin] = sum[bin]/num[bin];
    float rmssq = sumsq[bin]/num[bin] - mean[bin]*mean[bin];
    if( rmssq > 0.0 ) rms[bin] = sqrt( rmssq );
  }
  updated = true;

}

float profile::getMean( int bin ){
  if( ! updated ) calculate();
  return mean[bin];
}
float profile::getRMS( int bin ){
  if( ! updated ) calculate();
  return rms[bin];
}

histogram::histogram( const char* cname, const char* ctitle, const int bin, const float min, const float max ){

  name = std::string(cname);
  title = std::string(ctitle);
  nbins = bin;
  xmin = min;
  xmax = max;
  sum = 0.0;
  binSize = (xmax-xmin) / nbins;
  for( int bin=0; bin<nbins+2; bin++) contents.push_back( 0.0 );

}

int histogram::getBin( const float x ){

  int bin=-1;
  if( x < xmin ) bin = 0;
  else if( x >= xmax ) bin = nbins+1;
  else bin = int( ( x - xmin ) / binSize + 1 );

  return bin;
}

void histogram::Fill( const float x, const float val ){
  int bin = getBin( x );
  contents[bin] += val;
  sum += val;
}

void histogram::Add( const TH1* hist ){

  if( hist->GetNbinsX() != nbins ){ 
    std::cout << "histogram::Add, # of bin mismatch" << nbins << " " 
              << hist->GetNbinsX() << std::endl;
    return;
  }

  for( int bin=0; bin<nbins+2; bin++)
    contents[bin] += hist->GetBinContent( bin );
  sum += hist->Integral();

}

void histogram::save(){

  TH1F* hist = new TH1F( name.c_str(), title.c_str(), nbins, xmin, xmax );
  for( int bin=0; bin<nbins+2; bin++) 
    hist->SetBinContent( bin, contents[bin] );  
  hist->Write(0, TObject::kOverwrite);
}

float histogram::Integral( const int bin1, const int bin2 ){

  float integral=0.0;
  int min = bin1;
  if( min < 0 ) min = 0;
  int max = bin2+1;
  if( max > nbins ) max = nbins;
  if( max < min ){ 
    std::cout << "histogram::Integral, invalid integration region: "
              << min << " " << max << std::endl;
    return 0.0;
  }
  for( int bin=min; bin<max; bin++) integral += contents[bin];
  return integral;
}


void towerVar::readHists( TFile* hfile, UInt_t iRoot, UInt_t nRoot ){

  TH1F* hist, *rhist, *dhist, *ehist, *thist, *lhist;
  char name[] = "roccT00X17w3t4";
  
  sprintf(name,"numHitGTRCT%d", towerId);
  hist = (TH1F*)hfile->FindObjectAny( name );
  if( hist ) numHitGTRC->Add( hist );
  sprintf(name,"resXT%d", towerId);
  hist = (TH1F*)hfile->FindObjectAny( name );
  if( hist ) resX->Add( hist );  
  sprintf(name,"resYT%d", towerId);
  hist = (TH1F*)hfile->FindObjectAny( name );
  if( hist ) resY->Add( hist );  
#ifdef ROOT_PROFILE
  sprintf(name,"resProfT%d", towerId);
  TProfile *prof = (TProfile*)hfile->FindObjectAny( name );
  if( prof ) resProf->Add( prof );  
#endif

  for( int unp=0; unp!=g_nUniPlane; unp++){
    layerId lid( unp );
    std::string lname = lid.getLayerName();
    sprintf(name,"roccT%d%s", towerId, lname.c_str());
    rhist = (TH1F*)hfile->FindObjectAny( name );
    sprintf(name,"doccT%d%s", towerId, lname.c_str());
    dhist = (TH1F*)hfile->FindObjectAny( name );
    for( int strip=0; strip!=g_nStrip; strip++){      
      rHits[unp][strip] += (int)rhist->GetBinContent( strip+1 );
      dHits[unp][strip] += (int)dhist->GetBinContent( strip+1 );
    }
  }

  sprintf(name,"lhitT%d", towerId );
  dhist = (TH1F*)hfile->FindObjectAny( name );
  sprintf(name,"ltrkT%d", towerId );
  rhist = (TH1F*)hfile->FindObjectAny( name );
  for( UInt_t unp=0; unp!=bsVar.size(); unp++){
    bsVar[unp].hLayer += (int)dhist->GetBinContent( unp+1 );
    bsVar[unp].tLayer += (int)rhist->GetBinContent( unp+1 );
    layerId lid( unp );
    std::string lname = lid.getLayerName();
    sprintf(name,"eoccT%d%s", towerId, lname.c_str());
    ehist = (TH1F*)hfile->FindObjectAny( name );
    sprintf(name,"toccT%d%s", towerId, lname.c_str());
    thist = (TH1F*)hfile->FindObjectAny( name );
    sprintf(name,"loccT%d%s", towerId, lname.c_str());
    lhist = (TH1F*)hfile->FindObjectAny( name );
    for( int strip=0; strip!=g_nStrip; strip++){
      bsVar[unp].eHits[strip] += (int)ehist->GetBinContent( strip+1 );
      bsVar[unp].tHits[strip] += (int)thist->GetBinContent( strip+1 );
      bsVar[unp].lHits[strip] += (int)lhist->GetBinContent( strip+1 );
    }
    for(int iWafer = 0; iWafer != g_nWafer; ++iWafer){
      sprintf(name,"occT%d%sw%d", towerId, lname.c_str(), iWafer );
      hist = (TH1F*)hfile->FindObjectAny( name );
      if( hist ){ // check if simple version exist
        int tdiv = (iRoot*g_nTime)/nRoot;
        for( int strip=0; strip!=g_nStrip; strip++)
          bsVar[unp].nHits[strip][iWafer][tdiv] 
            += (int)hist->GetBinContent( strip+1 );
        continue; // no need to get time dependent histograms
      }
      for( int tDiv = 0; tDiv != g_nTime; tDiv++){
        sprintf(name,"occT%d%sw%dt%d", towerId, lname.c_str(), iWafer, tDiv );
        hist = (TH1F*)hfile->FindObjectAny( name );
        int tdiv = (tDiv+iRoot*g_nTime)/nRoot;
        for( int strip=0; strip!=g_nStrip; strip++)
          bsVar[unp].nHits[strip][iWafer][tdiv] 
            += (int)hist->GetBinContent( strip+1 );
      }
    }
  }
}

//
// TkrHits implementation 
//
TkrHits::TkrHits( bool initHistsFlag ): 
  m_numErrors(0), m_startTime(-1.0), m_endTime(-1.0), 
  m_firstRunId(0), m_lastRunId(0),
  m_reconEvent(0), m_digiEvent(0), m_rootFile(0),
  m_peakMIP(4.92), m_totAngleCF(5.48E-1), m_RSigma(4.0), m_GFrac(0.78),
  m_nDiv(1), m_correctedTot(true), m_histMode(true), m_badStrips(true),
  m_maxDirZ(-0.85), m_maxTrackRMS(0.3), m_maxDelta(3.0), m_trackRMS(-1.0)
{

  // get version number from CVS string
  std::string tag = "$Name: calibTkrUtil-02-10-01 $";
  int i = tag.find( " " );
  tag.assign( tag, i+1, tag.size() );
  i = tag.find( " " );
  tag.assign( tag, 0, i ) ;
  m_tag = tag;

  std::string version = "$Revision: 1.23 $";
  i = version.find( " " );
  version.assign( version, i+1, version.size() );
  i = version.find( " " );
  version.assign( version, 0, i ) ;
  m_version = version;
  m_revision = atof( m_version.c_str() );
  std::cout << "TkrHits, Tag: " << m_tag 
            << ", version: " << m_version << std::endl;

  for(int tower = 0; tower != g_nTower; ++tower)
    m_towerPtr[tower] = -1;

  if( initHistsFlag ){
    initCommonHists();
    if( m_badStrips ) initOccHists();
    initTotHists();

  }
#ifdef DEBUG_PRINT
  std::cout << "TkrHits Initialization finished." << std::endl;
#endif
}


void TkrHits::initCommonHists(){
  std::cout << "initialize common hisograms." << std::endl;
  m_nTrackDist = new TH1F("nTrack", "nTrack", 10, 0, 10);
  m_maxHitDist = new TH1F("maxHit", "maxHit", g_nUniPlane, 0, g_nUniPlane);
  m_trkRMS = new TH1F("trkRMS", "trkRMS", 100, 0, 2.0);
  m_numHitGTRC = new TH1F("numHitGTRC", "numHitGTRC", 65, 0, 65);
  m_numHitGTRCHE = new TH1F("numHitGTRCHE", "numHitGTRC (E>10GeV)", 65, 0, 65);
  m_largeMulGTRC = new TH1F("largeMulGTRC", "largeMulGTRC", g_nUniPlane, 0, g_nUniPlane);
  m_rawTOT = new TH1F("rawTOT", "rawTOT", 128, 0, 256);
  m_totT0X7L = new TH1F("totT0X7L", "totT0X7L", 128, 0, 256);
  m_totT0X7H = new TH1F("totT0X7H", "totT0X7H", 128, 0, 256);
  m_totT0X7TL = new TH1F("totT0X7TL", "totT0X7TL", 128, 0, 256);
  m_totT0X7TH = new TH1F("totT0X7TH", "totT0X7TH", 128, 0, 256);
  m_trkRMS1TWR = new TH1F("trkRMS1TWR", "trkRMS1TWR", 100, 0, 2.0);
  m_trkRMS2TWR = new TH1F("trkRMS2TWR", "trkRMS2TWR", 100, 0, 2.0);
  m_rmsProf1TWR = new TProfile("rmsProf1TWR","rmsProf1TWR",g_nTower,0,g_nTower);
  m_rmsProf2TWR = new TProfile("rmsProf2TWR","rmsProf2TWR",g_nTower,0,g_nTower);
  m_tresProfX = new TProfile("tresProfX","tresProfX",g_nTower,0,g_nTower);
  m_tresProfY = new TProfile("tresProfY","tresProfY",g_nTower,0,g_nTower);
  m_numClsDist = new TH1F("numCls", "# of cluster per layer", 10, 0, 10 );
  m_dirzDist = new TH1F("dirZ", "dirZ", 100, -1, 1);
  m_armsDist = new TH1F("arms", "arms", 100, -5, 5);
  m_res = new TH1F("resT0Y2", "residual for T0Y2", 100, -1, 1);
  m_resSel = new TH1F("resT0Y2Sel", "residual for T0Y2 (selection)", 100, -1, 1);
  m_lrec = new TH1F("lrec", "lrec", g_nUniPlane, 0, g_nUniPlane);
  m_ldigi = new TH1F("ldigi", "ldigi", g_nUniPlane, 0, g_nUniPlane);
  m_lcls = new TH1F("lcls", "lcls", g_nUniPlane, 0, g_nUniPlane);

  m_sigDist = new TProfile("msigDist", "msigDist", 100, 0, 600.0);
  m_sigRMS = new TProfile("msigRMS", "msigRMS", 100, 0, 0.5);
  m_sigTrad = new TProfile("msigTrad", "msigTrad", 100, 0, 10.0);


  //
  // MIP filter hists
  //
  Float_t bins[21];
  for( int i=0; i<21; i++) bins[i] = 20.0*pow(Float_t(10.),Float_t(0.2*i));
  m_hAcdTileCount = new TH1F("acdTileCount", "AcdTileCount", 10, 0 , 10 );
  m_hAcdTotalEnergy = new TH1F("acdTotalEnergy", "AcdTotalEnergy", 100, 0 , 10 );
  m_hCalEnergyRaw = new TH1F("calEnergyRaw", "calEnergyRaw", 100, 0, 1000 );
  m_hCalTotalEnergyRaw = new TH1F("calTotalEnergyRaw", "calTotalEnergyRaw", \
                                  20, bins );
  m_hNumCalXtal = new TH1F("numCalXtal", "# of Cal Xtal with energy>0 per layer", 10, 0, 10 );
  m_HAcdTileCount = new TH1F("acdTileCountC", "AcdTileCount after track cut", 10, 0 , 10 );
  m_HAcdTotalEnergy = new TH1F("acdTotalEnergyC", "AcdTotalEnergy after track cut", 100, 0 , 10 );
  m_HCalEnergyRaw = new TH1F("calEnergyRawC", "calEnergyRaw after track cut", 100, 0, 1000 );
  m_HNumCalXtal = new TH1F("numCalXtalC", "# of Cal Xtal with energy>0 per layer after track cut", 10, 0, 10 );

  //
  // siz in a row histograms
  //
  m_sixInARow = new TH1F( "sixInARow", "# of 6-in-a-row per Tower", 
                          g_nTower, 0, g_nTower );
  m_sixInARowWithTrig = new TH1F( "sixInARowWithTrig",
                                  "# of 6-in-a-row with trigger", 
                                  g_nTower, 0, g_nTower );
  m_orphanTrig = new TH1F( "orphanTrig", "# of triggers without 6-in-a-row", 
                           g_nTower, 0, g_nTower );
  m_sixInARowMIP = new TH1F( "sixInARowMIP", 
                             "# of 6-in-a-row per Tower (MIP)", 
                             g_nTower, 0, g_nTower );
  m_sixInARowT0X7 = new TH1F( "sixInARowT0X7", 
                             "# of 6-in-a-row per Tower (T0X7)", 
                             g_nTower, 0, g_nTower );
  m_sixInARowT3X7 = new TH1F( "sixInARowT3X7", 
                             "# of 6-in-a-row per Tower (T3X7)", 
                             g_nTower, 0, g_nTower );
  m_sixInARowWithTrigMIP = new TH1F( "sixInARowWithTrigMIP",
                                     "# of 6-in-a-row with trigger (MIP)", 
                                     g_nTower, 0, g_nTower );
  m_sixInARowWithTrigT0X7 = new TH1F( "sixInARowWithTrigT0X7",
                                     "# of 6-in-a-row with trigger (T0X7)", 
                                     g_nTower, 0, g_nTower );
  m_sixInARowWithTrigT3X7 = new TH1F( "sixInARowWithTrigT3X7",
                                     "# of 6-in-a-row with trigger (T3X7)", 
                                     g_nTower, 0, g_nTower );
  m_orphanTrigMIP = new TH1F( "orphanTrigMIP", 
                              "# of triggers without 6-in-a-row (MIP)", 
                              g_nTower, 0, g_nTower );
  m_sixInARowAll = new TH1F( "sixInARowAll", 
                             "# of 6-in-a-row per Tower (MIP, All)", 
                             g_nTower, 0, g_nTower );
  m_sixInARowWithTrigAll = new TH1F( "sixInARowWithTrigAll",
                                     "# of 6-in-a-row with trig (MIP, All)", 
                                     g_nTower, 0, g_nTower );
  m_sixInARowCut = new TH1F( "sixInARowCut", 
                             "# of 6-in-a-row per Tower (MIP, Cut)", 
                             ncut, 0, ncut );
  m_sixInARowWithTrigCut = new TH1F( "sixInARowWithTrigCut",
                                     "# of 6-in-a-row with trig (MIP, Cut)", 
                                     ncut, 0, ncut );
  m_hitCut = new TH1F( "hitCut", "# of hits with cut", ncut, 0, ncut );
  m_trackCut = new TH1F( "trackCut", "# of tracks with cut", ncut, 0, ncut );

  m_trigComb = new TH1F("trigComb", "trigger combination", \
                        g_nTkrLayer, 0, g_nTkrLayer );
  m_totTrig = new TH1F("totTrig", "totTrig", 256, 0, 256);
  m_totNonTrig = new TH1F("totNonTrig", "totNonTrig", 256, 0, 256);
  m_nBadLayersNonTrig = new TH1F("nBadLayersNonTrig", "nBadLayersNonTrig", 10, 0, 10);
  m_nBadLayersTrig = new TH1F("nBadLayersTrig", "nBadLayersTrig", 10, 0, 10);
  m_deltaWindOpenTime = new TH1F( "deltaWindOpenTime", "DeltaWindowOpenTime", 100, 0, 500 );
}


void TkrHits::initOccHists(){
  std::cout << "initialize occupancy hisograms." << std::endl;
  m_locc = new TH1F("locc", "locc", g_nUniPlane, 0, g_nUniPlane);
  m_ltrk = new TH1F("ltrk", "ltrk", g_nUniPlane, 0, g_nUniPlane);
  m_dist = new TH1F("dist", "distance", 50, 0, 200);
  m_brmsDist[0] = new TH1F("brms0", "brms 0-2", 100, -5, 5);
  m_brmsDist[1] = new TH1F("brms1", "brms 3-5", 100, -5, 5);
  m_brmsDist[2] = new TH1F("brms2", "brms 6-8", 100, -5, 5);
  m_brmsDist[3] = new TH1F("brms3", "brms 9-11", 100, -5, 5);
  m_brmsDist[4] = new TH1F("brms4", "brms 12-14", 100, -5, 5);
  m_brmsDist[5] = new TH1F("brms5", "brms 15-17", 100, -5, 5);
  m_occDist = new TH1F("occDist", "occDist", 200, 0, 200);
  m_poissonDist = new TH1F("poissonDist", "poissonDist", 40, -20, 0);
  m_aPos[0] = new TH1F("apos0", "apos0", 100, -250, 250);
  m_aPos[1] = new TH1F("apos1", "apos1", 100, -250, 250);
  m_aPos[2] = new TH1F("apos2", "apos2", 100, -250, 250);
  m_aPos[3] = new TH1F("apos3", "apos3", 100, -250, 250);
  m_aPos[4] = new TH1F("apos4", "apos4", 100, -250, 250);
}


void TkrHits::initTotHists(){
  std::cout << "initialize TOT hisograms." << std::endl;
  m_fracBatTot = new TH1F("fracBadTot", "fraction of bad TOT", 50, 0, 0.2 );
  m_fracErrDist = new TH1F("fracErrDist", "Peak error", 100, 0, 0.1);
  m_chisqDist = new TH1F("chisqDist", "TOT fit chisq/ndf", 60, 0, 3);
  m_chargeScale = new TH1F("chargeScale", "Charge Scale", 50, 0.5, 1.5);
  m_entries = new TH1F("entries", "Entries", 200, 0, 2000);
  m_langauWidth = new TH1F("langauWidth", "Langau Width", 50, 0.1, 0.6);
  m_langauGSigma = new TH1F("langauGSigma", "Langau GSigma", 50, 0.0, 2.0);
  m_dirProfile = new TProfile("dirProfile", "cons(theta) profile", 10, -1, -0.5);

  m_chist[0] = new TH1F("charge0", "TOT charge distribution (-1<cos<-0.95)", nTotHistBin, 0, maxTot);
  m_chist[1] = new TH1F("charge1", "TOT charge distribution (-0.95<cos<-0.9)", nTotHistBin, 0, maxTot);
  m_chist[2] = new TH1F("charge2", "TOT charge distribution (-0.9<cos<-0.85)", nTotHistBin, 0, maxTot);
  m_chist[3] = new TH1F("charge3", "TOT charge distribution (-0.85<cos<-0.8)", nTotHistBin, 0, maxTot);
  m_chist[4] = new TH1F("chargeAllMIP", "TOT charge distribution with MIP", nTotHistBin, 0, maxTot);
  m_chist[5] = new TH1F("chargeAll", "TOT charge distribution", nTotHistBin, 0, maxTot);
  m_chist[6] = new TH1F("chargeCorrected", "TOT charge distribution", nTotHistBin, 0, maxTot);

}


void TkrHits::saveAllHist( bool saveWaferOcc, bool runFitTot )
{
  if( m_rootFile == 0 ) return;
  std::cout << "save histograms" << std::endl;
  gDirectory->mkdir( "Towers" );

  float numEvents = m_nTrackDist->Integral(0,11);
  if( numEvents <= 0.0 ) numEvents = 1.0;
  if( m_nEvents > 0.0 ){
    m_largeMulGTRC->Scale( 1.0/(g_nTower*2*numEvents) );
    numEvents = m_nEvents;
  }

  m_nTrackDist->Write(0, TObject::kOverwrite);
  m_maxHitDist->Write(0, TObject::kOverwrite);
  m_numHitGTRC->Write(0, TObject::kOverwrite);
  m_numHitGTRCHE->Write(0, TObject::kOverwrite);
  m_largeMulGTRC->Write(0, TObject::kOverwrite);
  m_rawTOT->Write(0, TObject::kOverwrite);
  m_totT0X7L->Write(0, TObject::kOverwrite);
  m_totT0X7H->Write(0, TObject::kOverwrite);
  m_totT0X7TL->Write(0, TObject::kOverwrite);
  m_totT0X7TH->Write(0, TObject::kOverwrite);
  m_trkRMS->Write(0, TObject::kOverwrite);
  m_trkRMS1TWR->Write(0, TObject::kOverwrite);
  m_trkRMS2TWR->Write(0, TObject::kOverwrite);
  m_rmsProf1TWR->Write(0, TObject::kOverwrite);
  m_rmsProf2TWR->Write(0, TObject::kOverwrite);
  m_tresProfX->Write(0, TObject::kOverwrite);
  m_tresProfY->Write(0, TObject::kOverwrite);
  m_numClsDist->Write(0, TObject::kOverwrite);
  m_dirzDist->Write(0, TObject::kOverwrite);
  m_armsDist->Write(0, TObject::kOverwrite);
  m_res->Write(0, TObject::kOverwrite);
  m_resSel->Write(0, TObject::kOverwrite);
  m_lrec->Write(0, TObject::kOverwrite);
  m_ldigi->Write(0, TObject::kOverwrite);
  m_lcls->Write(0, TObject::kOverwrite);

  TH1F* sigDist = new TH1F("sigDist","sigDist", 100, 0, 600.0);
  TH1F* sigRMS = new TH1F("sigRMS","sigRMS", 100, 0, 0.5);
  TH1F* sigTrad = new TH1F("sigTrad","sigTrad", 100, 0, 10.0);
  for( int i=1; i<101; i++ ){
    sigDist->SetBinContent( i, m_sigDist->GetBinError( i ) );
    sigRMS->SetBinContent( i, m_sigRMS->GetBinError( i ) );
    sigTrad->SetBinContent( i, m_sigTrad->GetBinError( i ) );
  }
  sigDist->Write(0, TObject::kOverwrite);
  sigRMS->Write(0, TObject::kOverwrite);
  sigTrad->Write(0, TObject::kOverwrite);

  //
  // save MIP filter related stuff
  //
  m_hAcdTileCount->Write(0, TObject::kOverwrite);
  m_hAcdTotalEnergy->Write(0, TObject::kOverwrite);
  m_hCalEnergyRaw->Write(0, TObject::kOverwrite);
  m_hCalTotalEnergyRaw->Write(0, TObject::kOverwrite);
  m_hNumCalXtal->Write(0, TObject::kOverwrite);
  m_HAcdTileCount->Write(0, TObject::kOverwrite);
  m_HAcdTotalEnergy->Write(0, TObject::kOverwrite);
  m_HCalEnergyRaw->Write(0, TObject::kOverwrite);
  m_HNumCalXtal->Write(0, TObject::kOverwrite);
  //
  // save trgger eff related stuff
  //
  m_sixInARow->Write(0, TObject::kOverwrite);
  m_sixInARowWithTrig->Write(0, TObject::kOverwrite);
  m_orphanTrig->Write(0, TObject::kOverwrite);
  m_sixInARowMIP->Write(0, TObject::kOverwrite);
  m_sixInARowWithTrigMIP->Write(0, TObject::kOverwrite);
  m_orphanTrigMIP->Write(0, TObject::kOverwrite);
  m_sixInARowAll->Write(0, TObject::kOverwrite);
  m_sixInARowWithTrigAll->Write(0, TObject::kOverwrite);
  m_sixInARowCut->Write(0, TObject::kOverwrite);
  m_sixInARowWithTrigCut->Write(0, TObject::kOverwrite);
  m_trigComb->Write(0, TObject::kOverwrite);
  m_totTrig->Write(0, TObject::kOverwrite);
  m_totNonTrig->Write(0, TObject::kOverwrite);
  m_nBadLayersNonTrig->Write(0, TObject::kOverwrite);
  m_nBadLayersTrig->Write(0, TObject::kOverwrite);
  m_deltaWindOpenTime->Write(0, TObject::kOverwrite);
  m_hitCut->Write(0, TObject::kOverwrite);
  m_trackCut->Write(0, TObject::kOverwrite);
  //
  // save timestamp TTree
  //
  TTree *tree = new TTree("timeStamps","time stamps");
  tree->Branch("revision",&m_revision,"revision/D");
  tree->Branch("startTime",&m_startTime,"startTime/D");
  tree->Branch("endTime",&m_endTime,"endTime/D");
  tree->Branch("firstRunId",&m_firstRunId,"firstRunId/i"); // unsigned int
  tree->Branch("lastRunId",&m_lastRunId,"lastRunId/i"); // unsigned int
  tree->Fill();
  tree->Write();
  

  if( runFitTot ) fitTot();
  m_fracBatTot->Write(0, TObject::kOverwrite);
  m_fracErrDist->Write(0, TObject::kOverwrite);
  m_chisqDist->Write(0, TObject::kOverwrite);
  m_chargeScale->Write(0, TObject::kOverwrite);
  m_entries->Write(0, TObject::kOverwrite);
  m_langauWidth->Write(0, TObject::kOverwrite);
  m_langauGSigma->Write(0, TObject::kOverwrite);
  m_dirProfile->Write(0, TObject::kOverwrite);
  for( int i=0; i!=7; i++)
    m_chist[i]->Write(0, TObject::kOverwrite);

  if( gDirectory->cd( "ChargeDist" ) )
    gDirectory->Write();
  else {
    gDirectory->mkdir( "ChargeDist" );
    gDirectory->cd( "ChargeDist" );
    UInt_t nhist = m_chargeHist.size() / g_nTower;
    char rdir[]= "T16";
    for( unsigned int tw=0; tw<m_towerVar.size(); tw++ ){
      int tower = m_towerVar[ tw ].towerId;
      sprintf( rdir, "T%d", tower );
      gDirectory->mkdir( rdir );
      gDirectory->cd( rdir );
      for( UInt_t i=0; i!=nhist; i++ )
        m_chargeHist[tower*nhist+i]->Write(0, TObject::kOverwrite);
      gDirectory->cd( ".." );
    }
  }
  gDirectory->cd( ".." );
  //
  // save tot TTree
  //
  UInt_t rawTOT, icharge;
  tree = new TTree("totInfo","TOT information");
  tree->Branch("rawTOT",&rawTOT,"rawTOT/i"); // unsigned int
  tree->Branch("charge",&icharge,"charge/i"); // unsigned int
  for( UInt_t itot=0; itot<m_totInfo.size(); itot++){
    rawTOT = m_totInfo[itot].rawTOT;
    icharge = m_totInfo[itot].charge;
    tree->Fill();
  }
  tree->Write();

  //
  // calculate trigger and hit efficiencies
  //
  float eff, error, offset, frac, entry;
  TH1F* tEff = new TH1F( "trigEff", "trigger efficiency", ncut, 0, ncut );
  TH1F* hEff = new TH1F( "hitEff", "hit efficiency", ncut, 0, ncut );
  for( int icut=0; icut<ncut; icut++ ){
    entry = m_trackCut->GetBinContent( icut+1 );
    if( entry > 0 ){
      eff = m_hitCut->GetBinContent( icut+ 1 ) / entry;
      error = eff*(1-eff)/entry;
      if( error > 0.0 ) error = sqrt( error );
      else error = 1.0E-5;
      std::cout.setf( std::ios::fixed );
      std::cout.precision( 5 );
      std::cout << "LAT hit efficiency " << icut << ": " 
                << eff << " +- " << error << std::endl;
      hEff->SetBinContent( icut+1, eff );
      hEff->SetBinError( icut+1, error );
    }
    entry = m_sixInARowCut->GetBinContent( icut+1 );
    if( entry > 0 ){
      eff = m_sixInARowWithTrigCut->GetBinContent( icut+ 1 ) / entry;
      error = eff*(1-eff)/entry;
      if( error > 0.0 ) error = sqrt( error );
      else error = 1.0E-5;
      std::cout.setf( std::ios::fixed );
      std::cout.precision( 5 );
      std::cout << "LAT trigger efficiency " << icut << ": " 
                << eff << " +- " << error << std::endl;
      tEff->SetBinContent( icut+1, eff );
      tEff->SetBinError( icut+1, error );
    }
  }
  hEff->Write(0, TObject::kOverwrite);
  tEff->Write(0, TObject::kOverwrite);

  if( !m_badStrips ) return;

  saveOccHists();

  char name[] = "leffT00";
  TH1F* hist;
  TH1F* thits = new TH1F("thits", "thits", g_nTower, 0, g_nTower);
  TH1F* ttrks = new TH1F("ttrks", "ttrks", g_nTower, 0, g_nTower);
  TH1F* theff = new TH1F("theff", "Tower Hit Efficinecies", g_nTower, 0, g_nTower);
  TH1F* tteff = new TH1F("tteff", "Tower Trigger Efficiencies", g_nTower, 0, g_nTower);
  float bins[24] = {0.0001,0.01,0.014,0.02,0.028,0.04,0.056,0.08,0.12,0.17,0.25,0.35,0.5,0.7,1.0,1.4,2.0,2.8,4.0,5.6,8.0,12.0,16.0,24.0};
  TH1F* ineffDist = new TH1F("ineffDist", "ineffDist", 23, bins);
  for( UInt_t tw = 0; tw != m_towerVar.size(); ++tw){
    m_towerVar[tw].saveHists( saveWaferOcc );
    for( int unp = 0; unp != g_nUniPlane; unp++){
      ttrks->Fill( m_towerVar[tw].towerId, m_towerVar[tw].bsVar[unp].tLayer );
      thits->Fill( m_towerVar[tw].towerId, m_towerVar[tw].bsVar[unp].hLayer );
    }
    //
    // check for outliers
    //
    for( int unp=0; unp<g_nUniPlane; unp++){
      //
      // efficiency
      sprintf( name, "leffT%d", m_towerVar[tw].towerId );
      hist = (TH1F*)m_rootFile->FindObjectAny( name );
      eff = hist->GetBinContent( unp+1 );
      error = hist->GetBinError( unp+1 );
      if( error == 0.0 ) error = 1.0;
      ineffDist->Fill( 100*(1-eff)+0.0001 );
      if( (eff-0.98)/error<-3.0 && m_log.is_open() ){
        layerId lid( unp );
        std::string lname = lid.getLayerName();
        std::cout << m_towerVar[tw].hwserial << " T" << m_towerVar[tw].towerId 
                  << " " << lname
                  << " low efficiency: " << eff << " +- " << error << std::endl;
        m_log << m_towerVar[tw].hwserial << " T" << m_towerVar[tw].towerId 
              << " " << lname
              << " low efficiency: " << eff << " +- " << error << std::endl;
      }
      //
      // check for alignment outliers
      offset = m_towerVar[tw].resProf->GetBinContent( unp+1 );
      error = m_towerVar[tw].resProf->GetBinError( unp+1 );
      if( error == 0.0 ) error = 1.0;
      if( (fabs(offset)-0.1)/error > 3.0 && m_log.is_open() ){
        layerId lid( unp );
        std::string lname = lid.getLayerName();
        std::cout << m_towerVar[tw].hwserial << " T" << m_towerVar[tw].towerId 
                  << " " << lname
                  << " large offset: " << offset << " +- " << error << std::endl;
        m_log << m_towerVar[tw].hwserial << " T" << m_towerVar[tw].towerId 
              << " " << lname
              << " large offset: " << offset << " +- " << error << std::endl;
      }
    }
    //
    // check noise flare (sum 14-64)
    frac = m_towerVar[tw].numHitGTRC->Integral( 15, 65 ) 
      / ( numEvents * g_nUniPlane * 2 );
    if( frac > 1.0E-4 && m_log.is_open() ){
         std::cout << m_towerVar[tw].hwserial << " T" << m_towerVar[tw].towerId 
                   << " high fraction of large GTRC hits: " << frac << std::endl;
         m_log << m_towerVar[tw].hwserial << " T" << m_towerVar[tw].towerId 
               << " high fraction of large GTRC hits: " << frac << std::endl;
    }
  }

  entry = ttrks->Integral();
  if( entry > 0 ){
    eff = thits->Integral() / entry;
    error = eff*(1-eff)/entry;
    if( error > 0.0 ) error = sqrt( error );
    else error = 1.0E-5;
    std::cout.setf( std::ios::fixed );
    std::cout.precision( 5 );
    std::cout << "LAT hit efficiency: " << eff << " +- " << error << std::endl;
  }

  entry = m_sixInARowMIP->Integral();
  if( entry > 0 ){
    eff = m_sixInARowWithTrigMIP->Integral() / entry;
    error = eff*(1-eff)/entry;
    if( error > 0.0 ) error = sqrt( error );
    else error = 1.0E-5;
    std::cout.setf( std::ios::fixed );
    std::cout.precision( 5 );
    std::cout << "LAT trigger efficiency: " << eff << " +- " << error 
              << std::endl;
  }
  
  tteff->Add( m_sixInARowWithTrigMIP );
  tteff->Divide( m_sixInARowMIP );
  for( UInt_t tw = 0; tw != m_towerVar.size(); ++tw){
    int twr = m_towerVar[tw].towerId;
    // calculate binomial error
    entry = ttrks->GetBinContent( twr+1 );
    eff = thits->GetBinContent( twr+1 ) / entry;
    error = 1.0E-5;
    if( entry > 0 ){
      error = eff*(1-eff)/entry;
      if( error > 0.0 ) error = sqrt( error );
      else error = 1.0E-5;
    }
    theff->SetBinContent( twr+1, eff );
    theff->SetBinError( twr+1, error );
    if( m_log.is_open() ){
      std::cout.setf( std::ios::fixed );
      std::cout.precision( 5 );
      std::cout << m_towerVar[tw].hwserial << " T" << twr
                << " tower efficiency: " << eff << " +- " << error 
                << std::endl;
      m_log.setf( std::ios::fixed );
      m_log.precision( 5 );
      m_log << m_towerVar[tw].hwserial << " T" << twr
            << " tower efficiency: " << eff << " +- " << error 
            << std::endl;
    }
    eff = tteff->GetBinContent( twr+1 );
    entry = m_sixInARowMIP->GetBinContent( twr+1 );
    error = eff*(1-eff)/entry;
    if( error > 0.0 ) error = sqrt( error );
    else error = 1.0E-5;
    tteff->SetBinError( twr+1, error );
  }
  thits->Write(0, TObject::kOverwrite);
  ttrks->Write(0, TObject::kOverwrite);
  theff->Write(0, TObject::kOverwrite);
  tteff->Write(0, TObject::kOverwrite);
  ineffDist->Write(0, TObject::kOverwrite);
  //
  //
  //
  TH1F* tteffT0X7 = new TH1F("tteffT0X7", "Tower Trigger Efficiencies", g_nTower, 0, g_nTower);
  tteffT0X7->Add( m_sixInARowWithTrigT0X7 );
  tteffT0X7->Divide( m_sixInARowT0X7 );
  for( UInt_t tw = 0; tw != m_towerVar.size(); ++tw){
    int twr = m_towerVar[tw].towerId;
    eff = tteffT0X7->GetBinContent( twr+1 );
    entry = m_sixInARowMIP->GetBinContent( twr+1 );
    error = eff*(1-eff)/entry;
    if( error > 0.0 ) error = sqrt( error );
    else error = 1.0E-5;
    tteffT0X7->SetBinError( twr+1, error );
    std::cout.setf( std::ios::fixed );
    std::cout.precision( 5 );
    std::cout << m_towerVar[tw].hwserial << " T" << twr
              << " T0X7 trigger efficiency: " << eff << " +- " << error 
              << std::endl;
  }
  tteffT0X7->Write(0, TObject::kOverwrite);
  m_sixInARowT0X7->Write(0, TObject::kOverwrite);
  m_sixInARowWithTrigT0X7->Write(0, TObject::kOverwrite);
  //
  TH1F* tteffT3X7 = new TH1F("tteffT3X7", "Tower Trigger Efficiencies", g_nTower, 0, g_nTower);
  tteffT3X7->Add( m_sixInARowWithTrigT3X7 );
  tteffT3X7->Divide( m_sixInARowT3X7 );
  for( UInt_t tw = 0; tw != m_towerVar.size(); ++tw){
    int twr = m_towerVar[tw].towerId;
    eff = tteffT3X7->GetBinContent( twr+1 );
    entry = m_sixInARowT3X7->GetBinContent( twr+1 );
    error = eff*(1-eff)/entry;
    if( error > 0.0 ) error = sqrt( error );
    else error = 1.0E-5;
    tteffT3X7->SetBinError( twr+1, error );
    std::cout.setf( std::ios::fixed );
    std::cout.precision( 5 );
    std::cout << m_towerVar[tw].hwserial << " T" << twr
              << " T3X7 trigger efficiency: " << eff << " +- " << error 
              << std::endl;
  }
  m_sixInARowT3X7->Write(0, TObject::kOverwrite);
  m_sixInARowWithTrigT3X7->Write(0, TObject::kOverwrite);
  tteffT3X7->Write(0, TObject::kOverwrite);

  //
  // print out tower offset
  if( m_log.is_open() ){
    for( UInt_t tw = 0; tw != m_towerVar.size(); ++tw){
      int twr = m_towerVar[tw].towerId;
      std::cout.precision( 3 );
      std::cout << m_towerVar[tw].hwserial << " T" << twr
                << " xshift: " << m_tresProfX->GetBinContent( twr+1 )
                << " +- " << m_tresProfX->GetBinError( twr+1 )
                << " yshift: " << m_tresProfY->GetBinContent( twr+1 )
                << " +- " << m_tresProfY->GetBinError( twr+1 ) << std::endl;
      m_log.precision( 3 );
      m_log << m_towerVar[tw].hwserial << " T" << twr
            << " xshift: " << m_tresProfX->GetBinContent( twr+1 )
            << " +- " << m_tresProfX->GetBinError( twr+1 )
            << " yshift: " << m_tresProfY->GetBinContent( twr+1 )
            << " +- " << m_tresProfY->GetBinError( twr+1 ) << std::endl;
    }
  }

  //
  // Error counts
  //
  if( m_log.is_open() && m_numErrors > 0 ){
    std::cout << "Number of data errors: " << m_numErrors << std::endl;
    m_log << "Number of data errors: " << m_numErrors << std::endl;
  }
}


void TkrHits::saveOccHists(){
  std::cout << "save occupancy histograms" << std::endl;
  for( int i=0; i<g_nLayer/3; i++) 
    m_brmsDist[i]->Write(0, TObject::kOverwrite);
  m_locc->Write(0, TObject::kOverwrite);
  m_ltrk->Write(0, TObject::kOverwrite);
  m_dist->Write(0, TObject::kOverwrite);
  m_occDist->Write(0, TObject::kOverwrite);
  m_poissonDist->Write(0, TObject::kOverwrite);
  m_aPos[0]->Write(0, TObject::kOverwrite);
  m_aPos[1]->Write(0, TObject::kOverwrite);
  m_aPos[2]->Write(0, TObject::kOverwrite);
  m_aPos[3]->Write(0, TObject::kOverwrite);
  m_aPos[4]->Write(0, TObject::kOverwrite);
}


void TkrHits::printEvent( int iEvent )
{
  std::cout << "Event#: " << iEvent << std::endl;
}


void TkrHits::analyzeEvent()
{
#ifdef PRINT_DEBUG
    std::cout << "start analyzeEvent" << std::endl;
#endif
    if( ! m_towerInfoDefined ) setTowerInfo();
#ifdef PRINT_DEBUG
    std::cout << "setTowerInfo done " << std::endl;
#endif
    if( ! MIPfilter() ) return;
#ifdef PRINT_DEBUG
    std::cout << "MIPfilter done " << std::endl;
#endif
    //MIPfilter(); // result of MIP is used in passCut()
    monitorTKR();
#ifdef PRINT_DEBUG
    std::cout << "monitorTKR done " << std::endl;
#endif

    if(! passCut()) return;
#ifdef PRINT_DEBUG
    std::cout << "passCut done " << std::endl;
#endif

    getReconClusters();
#ifdef PRINT_DEBUG
    std::cout << "getReconCluster done " << std::endl;
#endif
    getDigiClusters();
#ifdef PRINT_DEBUG
    std::cout << "getDigiCluster done " << std::endl;
#endif
    selectGoodClusters();
#ifdef PRINT_DEBUG
    std::cout << "selectGoodCluster done " << std::endl;
#endif
    
    fillOccupancy( 0 );
#ifdef PRINT_DEBUG
    std::cout << "fillOcc done " << std::endl;
#endif
    fillTot();
#ifdef PRINT_DEBUG
    std::cout << "fillTot done " << std::endl;
#endif

}


bool TkrHits::MIPfilter()
{
  //
  // get timestamp and runId
  //
  Double_t ts = m_digiEvent->getTimeStamp(); 
  if( m_startTime < 0 ) m_startTime = ts;
  if( ts > m_endTime ) m_endTime = ts;

  //RunInfo runInfo = m_digiEvent->getMetaEvent().run();
  //UInt_t runId = runInfo.id();
  UInt_t runId = m_digiEvent->getRunId();
  if( m_firstRunId == 0 ){
    m_firstRunId = runId;
    m_lastRunId = runId;
    //std::cout << "run# " << runId << " " << m_digiEvent->getRunId() << std::endl;
  }
  else if( runId > m_lastRunId ) m_lastRunId = runId;
  else if( runId < m_firstRunId ) m_firstRunId = runId;

  //
  // MIP filter
  // pre-selection based on GEM/# of tracks
  //
  m_MIPtot = true;
  m_MIPeff = true;
  m_HighEnergy = false;
  for( int icut=0; icut<ncut; icut++) m_cut[icut] = true;
  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();
  assert(tkrRecon != 0);
  TObjArray* tracks = tkrRecon->getTrackCol();
  // select only 1-track event
  if( tracks->GetEntries() != 1 ) m_cut[1] = false;
#ifdef PRINT_DEBUG
  std::cout << "# of tracks: " << tracks->GetEntries() << std::endl;
#endif
  TkrTrack* track = dynamic_cast<TkrTrack*>(tracks->At(0));
  if( track ) 
    if( track->getKalEnergy() < 700 ) m_cut[2] = false;
  if( (m_digiEvent->getGem().getConditionSummary()&0x3) != 3 ) 
    m_cut[3] = false;
  //if( (m_digiEvent->getGem().getConditionSummary()&0x7) != 7 ) 
  if( m_digiEvent->getGem().getCnoVector() != 0 ) 
    m_cut[11] = false;
  if( m_digiEvent->getGem().getDeltaWindowOpenTime() < 150 )
    m_cut[4] = false;

  //
  // check ACD variables
  //
  AcdReconV2* acdRecon = m_reconEvent->getAcdReconV2();
  assert(acdRecon!= 0);
  AcdEventTopology& eventTopology = acdRecon->getEventTopology();
  AcdEventTopology* pTopo = &eventTopology;
  assert(pTopo!=0);
  m_acdTileCount = eventTopology.getTileCount();
  m_acdTotalEnergy = eventTopology.getTotalTileEnergy();
  m_hAcdTileCount->Fill( m_acdTileCount );
  m_hAcdTotalEnergy->Fill( m_acdTotalEnergy );

  //
  // check CAL variables
  //
  CalRecon* calRecon = m_reconEvent->getCalRecon();
  assert(calRecon != 0);
  TObjArray* clusters = calRecon->getCalClusterCol();
  int numClusters = clusters->GetEntries();
#ifdef PRINT_DEBUG
  std::cout << "# of clusters: " << clusters->GetEntries() << std::endl;
#endif
  if(numClusters == 2 or numClusters == 3){
    std::cout << "Invalid # of clusters: " << clusters->GetEntries() << std::endl;
    numClusters = 0;
  }
  else if(numClusters > 3){
    numClusters -= 2; // uber and uber2
  }
  //Event::CalCluster* calCluster = pCals->front();
  m_calEnergyRaw = 0;
  int num = 0;
  if( numClusters > 0 ){
    for( int cl=0; cl!=numClusters; cl++){
      CalCluster* calCluster = dynamic_cast<CalCluster*>(clusters->At(cl));
      if(calCluster) {
	m_calEnergyRaw += calCluster->getMomParams().getEnergy();
	num++;
      }
    }
  }
  if( m_digiEvent->getGem().getCnoVector()==0 ){
    m_hCalTotalEnergyRaw->Fill( m_calEnergyRaw );
    if( m_calEnergyRaw >1E4 ) m_HighEnergy = true;
  }
  if( num>0 ) m_calEnergyRaw /= num;
  m_hCalEnergyRaw->Fill( m_calEnergyRaw );
  TObjArray* xtals = calRecon->getCalXtalRecCol();
  int numXtals = xtals->GetEntries();
#ifdef PRINT_DEBUG
  std::cout << "# of Xtals: " << xtals->GetEntries() << std::endl;
#endif
  const int maxXtal = 7;
  UInt_t numXtal[maxXtal];
  for( int i=0; i<maxXtal; i++) numXtal[i] = 0;
  for( int xtal=0; xtal!=numXtals; xtal++){
    CalXtalRecData* calXtal = dynamic_cast<CalXtalRecData*>(xtals->At(xtal));
    if( calXtal->getEnergy() <= 0 ) continue;
    CalXtalId xtalId = calXtal->getPackedId();
    Short_t layerId = xtalId.getLayer();
    if( layerId > 7 ) std::cout << "Xtal layerId: " << layerId << std::endl;
    numXtal[layerId]++;
  }
  m_numCalXtal = 0;
  for( int i=0; i<maxXtal; i++){
    if( numXtal[i] > m_numCalXtal ) m_numCalXtal = numXtal[i];
  }
  m_hNumCalXtal->Fill( m_numCalXtal );
  
  //if( m_acdTileCount < 1 || m_acdTileCount > 2 ) m_MIPfilter = false;
  //if( abs(m_acdTotalEnergy - 1.8) > 1.0 ) m_MIPfilter = false;
  //if( abs(m_calEnergyRaw - 110) > 50 ) m_MIPfilter = false;
  //if( m_numCalXtal > 2 ) m_MIPfilter = false;
  if( abs(m_calEnergyRaw - 110) > 90 ) m_cut[7] = false;
  if( abs(m_calEnergyRaw - 110) > 50 ) m_cut[8] = false;
  if( m_acdTileCount < 1 || m_acdTileCount > 2 ) m_cut[9] = false;
  if( m_numCalXtal > 2 ) m_cut[10] = false;


  //
  // find 6-in-a-row
  //
  const TObjArray* tkrDigiCol = m_digiEvent->getTkrDigiCol();
#ifdef PRINT_DEBUG
  std::cout << "# of tkrGigi: " << tkrDigiCol->GetEntries() << std::endl;
#endif
  if (tkrDigiCol){
    bool sixInARow[g_nTower], triggered[g_nTower], badTower[g_nTower];
    int nHitTC[g_nTower][g_nTkrLayer],
      nBadTC[g_nTower][g_nTkrLayer];
    UShort_t tkrVector = m_digiEvent->getGem().getTkrVector();
    for(UInt_t tower=0; tower<g_nTower; tower++){
      UShort_t tvector = (1<<tower);
      triggered[tower] = tkrVector & tvector; // triggered tower
      badTower[tower] = false;
      sixInARow[tower] = false;
      m_goodTrackTowerFlag[ tower ] = triggered[tower];
      for (UInt_t bilayer=0; bilayer<g_nTkrLayer; bilayer++){
        nHitTC[tower][bilayer] = 0;
        nBadTC[tower][bilayer] = 0;
      }
    }
    
    TIter tkrIter(tkrDigiCol);
    TkrDigi *tkrDigi = 0;
    Int_t numT0X7=0, numNonT0X7=0, numT3X7=0, numNonT3X7=0;
    while ( ( tkrDigi = (TkrDigi*)tkrIter.Next() ) ) { 
      Int_t tower = tkrDigi->getTower().id();
      Int_t bilayer = tkrDigi->getBilayer();
      UInt_t numHits = tkrDigi->getNumHits();
      if( numHits > 0 ){
        Int_t tot = tkrDigi->getToT(0);
        Int_t toth = tkrDigi->getToT(1);
        if( toth > tot ) tot = toth;
        //if( tot<5 || tot>252 ) m_cut[5] = true;
        if( tot<1 || tot>252 ) m_cut[5] = true;
        int cmin = bilayer - 2;
        if( cmin < 0 ) cmin = 0;
        int cmax = bilayer+1;
        if( cmax > g_nTkrLayer ) cmax = g_nTkrLayer;
        for(int tc=cmin; tc<cmax; tc++){ // loop trigger combinations
          nHitTC[tower][tc]++;
          if( tot<10 || tot>252 ) nBadTC[tower][tc]++;
          if( nHitTC[tower][tc] >= 6 ){
            sixInARow[tower] = true;
            m_trigComb->Fill( tc );
            if( abs(tc-6) < 2 ) // trigger combination invovling X7
              if( tower == 0 ) numT0X7++;
              else if( tower == 3 ) numT3X7++;
            else if( tower == 0 ) numNonT0X7++;
            else if( tower == 3 ) numNonT3X7++;
          }
        }
      }
    }

    tkrIter.Reset();
    while ( ( tkrDigi = (TkrDigi*)tkrIter.Next() ) ) { 
      Int_t tower = tkrDigi->getTower().id();
      if( !sixInARow[tower] ) continue;
      if( tkrDigi->getNumHits() < 1 ) continue;
      Int_t tot = tkrDigi->getToT(0);
      Int_t toth = tkrDigi->getToT(1);
      if( toth > tot ) tot = toth;
      if( triggered[tower] ) m_totTrig->Fill( tot );
      else m_totNonTrig->Fill( tot );
    }

    int nSixInARow = 0;
    for(UInt_t tower=0; tower<g_nTower; tower++){
      if( sixInARow[tower] )
        for (UInt_t bilayer=0; bilayer<g_nTkrLayer; bilayer++){
          if( triggered[tower] ) 
            m_nBadLayersTrig->Fill( nBadTC[tower][bilayer] );
          else m_nBadLayersNonTrig->Fill( nBadTC[tower][bilayer] );
          //if( nBadTC[tower][bilayer] > 1 ) badTower[tower] = true;
        }
      if( badTower[tower] ){
        sixInARow[tower] = false;
        m_goodTrackTowerFlag[ tower ] = false;
      }
      if( sixInARow[tower] ){
        nSixInARow++;
        // search for additional tower
        for(UInt_t tw=tower+1; tw<g_nTower; tw++){ 
          if( sixInARow[tw] ){
            int row1 = tower/4;
            int row2 = tw/4;
            int col1 = tower%4;
            int col2 = tw%4;
            bool adj = false;
            if( row1==row2 && abs(col1-col2)<2 ) adj = true;
            if( col1==col2 && abs(row1-row2)<2 ) adj = true;
            if( ! adj )  m_cut[6] = false;
            // any towers needs to be adjacent
          }
        }
      }
    }

    // combine cuts
    m_cut[12] = m_cut[1] && m_cut[5];
    m_cut[13] = m_cut[12] && m_cut[11];
    m_cut[14] = m_cut[13] && m_cut[3];
    m_cut[15] = m_cut[14] && m_cut[6];
    m_cut[16] = m_cut[15] && m_cut[4] && m_cut[2];
    m_cut[17] = m_cut[16] && m_cut[7];
    m_cut[18] = m_cut[17] && m_cut[8];
    m_cut[19] = m_cut[18] && m_cut[9] && m_cut[10];
    m_MIPtot = m_cut[16];
    m_MIPeff = m_cut[19];
    m_cut[20] = m_MIPeff && (numT0X7>0) && (numNonT0X7==0);
    m_cut[21] = m_MIPeff && (numT0X7>0) && (numNonT0X7>0);
    m_cut[22] = m_MIPeff && (numT0X7==0);
    m_cut[23] = m_MIPeff && (numT3X7>0) && (numNonT3X7==0);
    m_cut[24] = m_MIPeff && (numT3X7>0) && (numNonT3X7>0);
    m_cut[25] = m_MIPeff && (numT3X7==0);

    for(UInt_t tower=0;tower<g_nTower; tower++){
      if( sixInARow[tower] ){
        if( nSixInARow > 1 ){
          m_sixInARow->Fill( tower );
          if( m_MIPeff )
            m_sixInARowMIP->Fill( tower );
          if( m_MIPtot ){
            if( m_cut[20] ) m_sixInARowT0X7->Fill( tower );
            if( m_cut[23] ) m_sixInARowT3X7->Fill( tower );
          }
          // check various additional cuts
          for( int icut=0; icut<ncut; icut++ )
            if( m_cut[icut] ) m_sixInARowCut->Fill( icut+0.5 );
          if( triggered[tower] ){
            m_sixInARowWithTrig->Fill( tower );
            for( int icut=0; icut<ncut; icut++ )
              if( m_cut[icut] ) m_sixInARowWithTrigCut->Fill( icut+0.5 );
            if( m_MIPeff )
              m_sixInARowWithTrigMIP->Fill( tower );
            if( m_MIPtot ){
              if( m_cut[20] ) m_sixInARowWithTrigT0X7->Fill( tower );
              if( m_cut[23] ) m_sixInARowWithTrigT3X7->Fill( tower );
            }
          }
        }
        if( m_MIPeff ){
          m_sixInARowAll->Fill( tower );
          if( triggered[tower] )
            m_sixInARowWithTrigAll->Fill( tower );
        }
      }
      else if( triggered[tower] ){
        m_orphanTrig->Fill( tower );
        if( m_MIPeff ) m_orphanTrigMIP->Fill( tower );
      }
    }
  }

  return m_cut[11]; // veto CNO
}



float TkrHits::calcCharge( layerId lid, int iStrip, int tot)
{
  // convert TOT raw count to micro second
  float time = (tot << 2) * 0.05;

  int tw = m_towerPtr[ lid.tower ];
  int unp = lid.uniPlane;
  // TOT to charge conversion
  if( m_towerVar[tw].tcVar[unp].totThreshold[iStrip] < 0.0 ) return -1.0;
  float charge = m_towerVar[tw].tcVar[unp].totThreshold[iStrip] 
    + time*m_towerVar[tw].tcVar[unp].totGain[iStrip]
    + time*time*m_towerVar[tw].tcVar[unp].totQuad[iStrip];
  //std::cout << lid.tower << " " << unp << " " << iStrip << " " << tot << " " << time << " " << charge << std::endl;
  
  return charge;
}


void TkrHits::fillTot() 
{
  
#ifdef PRINT_DEBUG
  std::cout << "fillTot start" << std::endl;
#endif

  if( !m_MIPtot ) return;

  for( unsigned int cls=0; cls<m_clusters.size(); cls++){
    Cluster* cluster = m_clusters[cls];
    layerId lid = getLayerId( cluster );
    int tower = lid.tower;
    int unp = lid.uniPlane;
    int tw = m_towerPtr[ tower ];

    // require only a single strip
    if(cluster->getSize() != 1) continue;
    
    int iStrip = cluster->getFirstStrip();
    
    int tot = cluster->getRawToT();
    if( tot == 0 ) continue;
 
    float charge;
    if( m_correctedTot ) charge = cluster->getCorrectedTot();
    else charge  = calcCharge( lid, iStrip, tot);
    if( charge < 0.0 ) continue;
    // empirical correction factor
    charge /= ( 1 + m_totAngleCF * (1+m_dir.z()) ); 
    m_chist[5]->Fill( charge );
    if( m_cut[17] ) m_chist[4]->Fill( charge );
    int idirz = int( 20 + m_dir.z()*20 );
    if( idirz < 4 ) m_chist[idirz]->Fill( charge );
    m_dirProfile->Fill( m_dir.z(), m_dir.z() );

    // fill TOT per layer for most cases
    // or per TFE for charge scale calibration
    int ibin = int( charge * nTotHistBin / maxTot );
    int iDiv = 0;
    iDiv = iStrip * m_nDiv / g_nStrip;
    //std::cout << tw << " " << unp << " " << iDiv << " " << ibin << std::endl;
    if( ibin < nTotHistBin && ibin >=0 )
      m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin]++;

    // save raw and corrected TOT inforation for later processing.
    totInfo totinfo;
    totinfo.rawTOT = tot + maxRawTOT*iStrip + maxRawTOT*g_nStrip*unp 
      + maxRawTOT*g_nStrip*g_nUniPlane*tower;
    UInt_t cfac = int( 1000 /  ( 1 + m_totAngleCF * (1+m_dir.z()) ) );
    if( cfac > 1000 ) cfac = 1000;
    if( cfac < 0 ) cfac = 0;
    UInt_t icharge = int( 1000 * charge );
    if( icharge > 100000 ) icharge = 100000;
    if( icharge < 0 ) icharge = 0;
    totinfo.charge = cfac*100000 + icharge;
    m_totInfo.push_back( totinfo );
  }
  
#ifdef PRINT_DEBUG
  std::cout << "fillTot end" << std::endl;
#endif
}


void TkrHits::fitTot()
{  
  // define Gaussian convolved Laudau function.
  TF1 *ffit = new TF1( "langau2", langau2fun, 0, 30, 6 );
  ffit->SetParNames( "Width", "MP", "Area", "GSigma", "RSigma", "GFrac" );
  ffit->SetLineWidth( 2 );
  ffit->SetLineColor( 2 );
  std::cout << "Start fit." << std::endl;
  m_chargeHist.clear();
  
  gDirectory->mkdir( "ChargeDist" );
  gDirectory->cd( "ChargeDist" );
  char rdir[] = "T16";
  int totalCount = 0;
  
  const float meanChargeScale = 1.12, rangeChargeScale=0.3;
  for( unsigned int tw=0; tw<m_towerVar.size(); tw++ ){
    int tower = m_towerVar[ tw ].towerId;
    std::cout << "Tower " << tower << ": ";
    sprintf( rdir, "T%d", tower );
    gDirectory->mkdir( rdir );
    gDirectory->cd( rdir );

    for(int unp = 0; unp != g_nUniPlane; ++unp) {
      layerId lid( unp );
      std::string lname = lid.getLayerName();
      std::cout << " " << lname;
      bool alarm = false;
      for(int iDiv = 0; iDiv != m_nDiv; ++iDiv){
        m_towerVar[tw].tcVar[unp].chargeScale[iDiv] = meanChargeScale;
        
        float area, ave, rms;
        Double_t *par, *error;
        float peak, errPeak, width, errWidth;
          
        // fill charge for each FE
        char name[] = "chargeT00X00fe0000";
        if( m_nDiv != 1 )
          sprintf(name,"chargeT%d%sfe%d", tower, lname.c_str(), iDiv);
        else
          sprintf(name,"chargeT%d%s", tower, lname.c_str());
        TH1F* chargeHist = new TH1F(name, name, nTotHistBin, 0, maxTot);
        chargeHist->SetLineWidth( 2 );
        float binWidth = maxTot / nTotHistBin;
        for( int ibin=0; ibin!=nTotHistBin; ibin++){
          if( m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin] < 0 ){
            std::cout << "invalid count " << tower << " " << lname << " " 
                      << iDiv << " " << ibin << " " 
                      << m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin]
                      << std::endl;
            continue;
          }
          totalCount +=  m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin];
          chargeHist->Fill( (ibin+0.5)*binWidth, 
                            m_towerVar[tw].tcVar[unp].chargeDist[iDiv][ibin] );
        }
        m_chargeHist.push_back( chargeHist );
        //std::cout << "T" << tower << " " << lname << " " << iDiv << std::endl;
        //if( m_histMode ) continue; // do not fit during hist mode.

        // fit charge for each FE
        area = chargeHist->Integral();
        ave = chargeHist->GetMean();
        rms = chargeHist->GetRMS();
        m_entries->Fill( area );
        //std::cout << area << " " << ave << " " << rms << std::endl;
        if( area<200 || ave==0.0 || rms==0.0 ){ 
          m_log << "T" << tower << " " << lname
                << " " << iDiv << ", Entries: " << area
                << ", Mean: " << ave << ", RMS: " << rms 
                << " skipped." << std::endl;
          continue;
        }
        
        //
        // check the presence of bad TOT values
        // try not to include these TOT in the fit.
        //
        int bin = int(ave*0.5/binWidth) + 1;
        float fracBadTot = chargeHist->Integral(1,bin)
          / chargeHist->Integral();
        m_fracBatTot->Fill( fracBadTot );

        float lowLim = ave - 1.4 * rms;
        if( fracBadTot > 0.05 && lowLim < ave*0.5 ){
          lowLim = ave*0.5;
          alarm = true;
          std::cout << std::endl << "WARNING, large bad TOT fraction: " 
                    << fracBadTot << ", T" << tower << " " << lname
                    << " " << iDiv;
          m_log << "WARNING, large bad TOT fraction: " 
                << fracBadTot << ", T" << tower << " " << lname
                << " " << iDiv << std::endl;
        }

        ffit->SetParLimits( 0, 0.10, 0.5 );
        ffit->SetParLimits( 1, 0.0, ave*2 );
        ffit->SetParLimits( 2, 0.0, area*0.4 );
        ffit->SetParLimits( 3, 0.35, 0.85 );
        ffit->SetRange( lowLim, ave+2.5*rms );
        ffit->SetParameters( rms*0.2, ave*0.75, area*0.1, rms*0.4 );
        ffit->FixParameter( 4, m_RSigma );
        ffit->FixParameter( 5, m_GFrac );
        chargeHist->Fit( "langau2", "RBQ" );
        
        //0:width(scale) 1:peak 2:total area 3:width(sigma)
        par = ffit->GetParameters();
        error = ffit->GetParErrors();
        
        peak = float( *(par+1) );
        errPeak = float( *(error+1) );
        if( peak > 0 ) m_fracErrDist->Fill( errPeak*sqrt(area/1000)/peak );
        
        width = float( *(par+3) );
        errWidth = float( *(error+3) );

        float chisq = ffit->GetChisquare();
        float ndf = ffit->GetNDF();
        if( ndf > 0 ) m_chisqDist->Fill( chisq/ndf );
        if( peak > 0.0 ){
          float chargeScale =  m_peakMIP / peak;
          addScaledHist( chargeHist, m_chist[6], chargeScale );
          m_chargeScale->Fill( chargeScale );
          m_langauWidth->Fill( *(par+0) ); //  width (scale)
          m_langauGSigma->Fill( *(par+3) ); // width (sigma)
          if( fabs(chargeScale-meanChargeScale) > rangeChargeScale ){
            alarm = true;
            std::cout << std::endl << "WARNING, Abnormal charge scale: " 
                      << chargeScale << ", T" << tower << " " << lname
                      << " " << iDiv;
            m_log << "WARNING, Abnormal charge scale: "
                  << chargeScale << ", T" << tower << " " << lname 
                  << " " << iDiv << std::endl;
            if( chargeScale > meanChargeScale ) 
              chargeScale = meanChargeScale + rangeChargeScale;
            if( chargeScale < meanChargeScale ) 
              chargeScale = meanChargeScale - rangeChargeScale;
          }
          if( chisq/ndf > maxChisq ){ // large chisq/ndf
            alarm = true;
            std::cout << std::endl << "WARNING, large chisq/ndf: "
                      << chisq/ndf << ", T" << tower << " " << lname
                      << " " << iDiv;
            m_log << "WARNING, large chisq/ndf: "
                  << chisq/ndf << ", T" << tower << " " << lname
                  << " " << iDiv << std::endl;
          }
          // large peak fit error
          if( errPeak*sqrt(area/1000)/peak > maxFracErr 
              || errPeak*sqrt(area/1000)/peak < minFracErr ){ 
            alarm = true;
            std::cout << std::endl << "WARNING, abnormal peak fit error: "
                      << errPeak*sqrt(area/1000)/peak << ", T" << tower 
                      << " " << lname << " " << iDiv;
            m_log << "WARNING, abnormal peak fit error: "
                  << errPeak*sqrt(area/1000)/peak << ", T" << tower << " " 
                  << lname << " " << iDiv << std::endl;
          }
          m_towerVar[tw].tcVar[unp].chargeScale[iDiv] = chargeScale;
        }
        else{
          alarm = true;
          std::cout << std::endl << "WARNING, negative peak value: "
                    << peak << ", T" << tower << " " << lname 
                    << " " << iDiv;
          m_log << "WARNING, negative peak value: "
                << peak << ", T" << tower << " " << lname 
                << " " << iDiv << std::endl;
        }
        
        m_log << "Fit T" << tower << " " << lname << " " 
              << iDiv << ' ';
        m_log.precision(3);
        m_log << area << ' ' << ave << ' ' << rms << ", " << *(par+0) << ' '
              << *(par+1) << ' ' << *(par+2) << ' ' << *(par+3) << ", "
              << errPeak/peak << " " << int(chisq+0.5) << "/" << ndf
              << std::endl;
        
      }
      if( alarm ) std::cout << std::endl;
    }
    std::cout << std::endl;
    gDirectory->cd( ".." );
  }
  gDirectory->cd( ".." );
  //gDirectory->ls();
  //gDirectory->pwd();
  std::cout << "total TOT count: " << totalCount << std::endl;
  if( m_log.is_open() )
    m_log << "total TOT count: " << totalCount << std::endl;

  float GSigma, RSigma = m_RSigma, GFrac = m_GFrac;
  for( int i=6; i>=0; i--){
    float area = m_chist[i]->Integral();
    float ave = m_chist[i]->GetMean();
    float rms = m_chist[i]->GetRMS();
    ffit->SetParLimits( 0, 0.0, rms );
    ffit->SetParLimits( 1, 0.0, ave*2 );
    ffit->SetParLimits( 2, 0.0, area*0.4 );
    ffit->SetParLimits( 3, 0.0, rms );
    ffit->SetRange( ave-1.25*rms, ave+2*rms );
    if( i >= 5 ){
      GSigma = rms * 0.4;
      ffit->ReleaseParameter( 3 );
      ffit->ReleaseParameter( 4 );
      ffit->ReleaseParameter( 5 );
      ffit->SetParameters(rms*0.5, ave*0.75, area*0.1, GSigma, RSigma, GFrac);
    }
    else{
      ffit->SetParameters( rms*0.5, ave*0.75, area*0.1 );
      ffit->FixParameter( 3, GSigma );
      ffit->FixParameter( 4, RSigma );
      ffit->FixParameter( 5, GFrac );
    }
    m_chist[i]->Fit( "langau2", "RBQ" );
    m_chist[i]->SetLineWidth( 2 );
    if( i == 5 ){
      GSigma = ffit->GetParameter( 3 );
      RSigma = ffit->GetParameter( 4 );
      GFrac = ffit->GetParameter( 5 );
    }    
  }
}


void TkrHits::addScaledHist( TH1F* hist1, TH1F* hist2, float xscale ){

  for( int ibin=0; ibin<hist1->GetNbinsX(); ibin++){ 
    float count = hist1->GetBinContent( ibin+1 );
    if( count == 0.0 ) continue;
    Double_t low = hist1->GetBinLowEdge( ibin+1 ) * xscale;
    Double_t size = hist1->GetBinWidth( ibin+1 ) * xscale;
    Double_t high = low + size;
    Int_t ibin2 = hist2->FindBin( low );
    float sum = 0.0, cvalue = 0.0;
    while( true ){
      float count2 = hist2->GetBinContent( ibin2 );
      Double_t low2 = hist2->GetBinLowEdge( ibin2 );
      Double_t size2 = hist2->GetBinWidth( ibin2 );
      Double_t high2 = low2 + size2;
      if( ibin2 >= hist2->GetNbinsX()+1 ){
        hist2->SetBinContent( ibin2, count-sum+count2 );
        sum = count;
        break;
      }
      if( low2 >= high || high2 <= low ) break;
      else if( low>=low2 && high<=high2 ) cvalue = count;
      else if( low<=low2 && high>=high2 ) cvalue = count*size2/size;
      else if( low<=low2 && high<=high2 ) cvalue = count*(high-low2)/size;
      else if( low>=low2 && high>=high2 ) cvalue = count*(high2-low)/size;
      else
        std::cout << "orphan case: " << low << " " << high << " " << low2 << " " << high2 << std::endl;
      hist2->SetBinContent( ibin2, cvalue+count2 );
      sum += cvalue;
      ibin2++;
    }
    if( fabs(sum/count-1) > 1.0E-4 )
      std::cout << "invalid sum: " << sum << " " << count << std::endl;
  }

}


void TkrHits::setTowerInfo(){

  if( m_towerVar.size() == g_nTower ){
    m_towerInfoDefined = true;
    return;
  }

  const TObjArray* tkrDigiCol = m_digiEvent->getTkrDigiCol();
  if (!tkrDigiCol) return;
  
  std::vector<int> strips;
  // Loop over all TkrDigis
  TIter tkrIter(tkrDigiCol);
  TkrDigi *tkrDigi = 0;
  while ( ( tkrDigi = (TkrDigi*)tkrIter.Next() ) ) {
    Int_t tower = tkrDigi->getTower().id();
    // register tower id pointer if this is the first encounter
    //    if(tower!=2 && tower!=3)
    //      std::cout<<tower<<" "<<m_towerVar.size()<<std::endl;
    if( m_towerPtr[ tower ] < 0 ){
      m_towerPtr[ tower ] = m_towerVar.size();
      m_towerVar.push_back( towerVar( tower, m_badStrips ) );
    }
  }
}


layerId TkrHits::getLayerId( Cluster* cluster ){
  
  int tower = cluster->getTowerId();
  int unp = cluster->getUniPlane();

  if( tower < 0 || tower >= g_nTower || unp < 0 || unp >= g_nUniPlane ){
    m_numErrors++;
    layerId lid( 0 );
    lid.setTower( -1 );
    return lid;
  }
  layerId lid( unp );
  lid.setTower ( tower );
  return lid;
}


layerId TkrHits::getLayerId( const TkrCluster* cluster )
{

  commonRootData::TkrId id = cluster->getTkrId();
  int tower = TowerId( id.getTowerX(), id.getTowerY() ).id();
  int view = id.getView();
  int layer = cluster->getLayer();
  if( tower < 0 || tower >= g_nTower || layer < 0 || layer >= g_nLayer || view < 0 || view >= g_nView ){
    m_numErrors++;
    layerId lid( 0 );
    lid.setTower( -1 );
    return lid;
  }
  layerId lid( layer, view, tower);
  return lid;
}

void TkrHits::monitorTKR(){

  // Loop over all TkrDigis
  const TObjArray* tkrDigiCol = m_digiEvent->getTkrDigiCol();
  if (!tkrDigiCol) return;
  TIter tkrIter(tkrDigiCol);
  TkrDigi *tkrDigi = 0;
  while ( ( tkrDigi = (TkrDigi*)tkrIter.Next() ) ) { 
    Int_t totl = tkrDigi->getToT(0);
    if( totl > 0 ) m_rawTOT->Fill( totl );
    Int_t toth = tkrDigi->getToT(1);
    if( toth > 0 ) m_rawTOT->Fill( toth );
    Int_t lastRC0Strip = tkrDigi->getLastController0Strip();

    UInt_t numHits = tkrDigi->getNumHits();
    // Loop through collection of hit strips for this TkrDigi
    UInt_t ihit;
    Int_t numl=0, numh=0;
    for (ihit = 0; ihit < numHits; ihit++) {
      // Retrieve the strip number
      if( abs(tkrDigi->getStrip(ihit)-819) < 12 ) continue;
      if( tkrDigi->getStrip(ihit) > lastRC0Strip ) numh++;
      else numl++;
    }
    Int_t tower = tkrDigi->getTower().id();
    Int_t tw = -1;
    if( tower >= 0 && tower < g_nTower ) tw = m_towerPtr[ tower ];
    if( tw < 0 ){
      if( m_log.is_open() && m_numErrors == 0 ){
        std::cout << "Data Error, invalid tower ID: " << tower << std::endl;
        m_log << "Data Error, invalid tower ID: " << tower << std::endl;
        std::cout << "Data error is issued only once and total number of " 
                  << "occurnaces will be reported at the end of analysis." 
                  << std::endl;
        m_log << "Data error is issued only once and total number of " 
              << "occurnaces will be reported at the end of analysis." 
              << std::endl;
      }
      m_numErrors++;
      continue;
    }
    Int_t layer = tkrDigi->getBilayer(); 
    GlastAxis::axis viewId = tkrDigi->getView(); 
    int view = (viewId == GlastAxis::X) ? 0 : 1; 
    if( numl > 0 ){
      m_numHitGTRC->Fill( numl );
      if( m_HighEnergy )m_numHitGTRCHE->Fill( numl );
      m_towerVar[tw].numHitGTRC->Fill( numl );
    }
    if( numh > 0 ){
      m_numHitGTRC->Fill( numh );
      if( m_HighEnergy )m_numHitGTRCHE->Fill( numh );
      m_towerVar[tw].numHitGTRC->Fill( numh );
    }
    if( numl > maxHitGTRC || numh > maxHitGTRC ){
      if( layer < 0 || layer >= g_nLayer || view < 0 || view >= g_nView ){
        m_numErrors++;
        continue;
      }
      layerId lid( layer, view );
      if( numl > maxHitGTRC ) m_largeMulGTRC->Fill( lid.uniPlane );        
      if( numh > maxHitGTRC ) m_largeMulGTRC->Fill( lid.uniPlane );
    }
    if( m_log.is_open() && (numh>bufferSizeGTRC||numl>bufferSizeGTRC) ){
      char vw[] = "XY";
      std::cout << m_towerVar[tw].hwserial << "" << tower
                << vw[view] << layer << " invalid GTRC multiplicity: " 
                << numh << " " << numl << std::endl;
      m_log << m_towerVar[tw].hwserial << " T" << tower
            << vw[view] << layer << " invalid GTRC multiplicity: " 
            << numh << " " << numl << std::endl;
    }
  } 
}


void TkrHits::getDigiClusters()
{
#ifdef PRINT_DEBUG
  std::cout << "getDigiClusters start" << std::endl;
#endif
  //
  // check data errors
  if( m_numErrors > 1E4 ){
    std::cout << "Too many data errrors are observed: " << m_numErrors 
              << std::endl << "EXIT." << std::endl;
    exit( EXIT_FAILURE );
  }
  //
  //
  // clear cluster information.
  for( UInt_t tw=0; tw!=m_towerVar.size(); tw++)
    for( int unp=0; unp!=g_nUniPlane; unp++)
      m_towerVar[tw].digiClusters[unp].clear();
  
  // The full collection of TkrDigis for this event
  const TObjArray* tkrDigiCol = m_digiEvent->getTkrDigiCol();
  if (!tkrDigiCol) return;
  
  std::vector<int> strips;
  // Loop over all TkrDigis
  TIter tkrIter(tkrDigiCol);
  TkrDigi *tkrDigi = 0;
  while ( ( tkrDigi = (TkrDigi*)tkrIter.Next() ) ) {
    // Identify the tower and layer
    Int_t tower = tkrDigi->getTower().id();
    Int_t tw = m_towerPtr[ tower ];
    if( tw < 0 ) continue; // invalid trower ID
    Int_t totl = tkrDigi->getToT(0);
    Int_t toth = tkrDigi->getToT(1);
    Int_t lastRC0Strip = tkrDigi->getLastController0Strip();

    Int_t layer = tkrDigi->getBilayer();
    if( layer < 0 || layer >= g_nLayer ){
      m_numErrors++;
      continue;
    }
  
    // Returns the orientation of the strips
    GlastAxis::axis viewId = tkrDigi->getView();
    int view = (viewId == GlastAxis::X) ? 0 : 1;
    if( view < 0 || view >= g_nView ){
      m_numErrors++;
      continue;
    }

    layerId lid( layer, view );
    int uniPlane = lid.uniPlane;
    
    strips.clear();
    UInt_t numHits = tkrDigi->getNumHits();
    // Loop through collection of hit strips for this TkrDigi
    UInt_t ihit;
    for (ihit = 0; ihit < numHits; ihit++) {
      // Retrieve the strip number
      Int_t iStrip = tkrDigi->getStrip(ihit);
      if( iStrip < 0 || iStrip >= g_nStrip ){
        m_numErrors++;
        continue;
      }
      strips.push_back( iStrip );
      m_towerVar[tw].dHits[uniPlane][iStrip]++;      
    }
    std::sort( strips.begin(), strips.end() ); // sort strip# for clustering.
    for( UInt_t i=0; i!=strips.size(); i++){
      bool newCls = false;
      if( m_towerVar[tw].digiClusters[uniPlane].empty() ) newCls = true;
      else if( !m_towerVar[tw].digiClusters[uniPlane].back().addStrip( strips[i] ) )
        newCls = true;
      
      if( newCls ){
        m_ldigi->Fill( uniPlane );
        if( strips[i] <= lastRC0Strip )
          m_towerVar[tw].digiClusters[uniPlane].push_back( Cluster( strips[i], totl ) );
        else
          m_towerVar[tw].digiClusters[uniPlane].push_back( Cluster( strips[i], toth ) );
      }
    }
  }
#ifdef PRINT_DEBUG
  std::cout << "getDigiClusters end" << std::endl;
#endif
}

void TkrHits::getReconClusters()
{
#ifdef PRINT_DEBUG
  std::cout << "getReconClusters start" << std::endl;
#endif
  
  // initialize recon cluster info
  for( unsigned int tw=0; tw<m_towerVar.size(); tw++){
    int twr = m_towerPtr[ m_towerVar[tw].towerId ];
    if( twr != int(tw) ) {
      std::cout << "Invalid tower id: " << twr << " != " << tw << std::endl;
      if( m_log.is_open() )
        m_log << "Invalid tower id: " << twr << " != " << tw << std::endl;
      exit( EXIT_FAILURE );
    }
    for( int unp=0; unp<g_nUniPlane; unp++)
      m_towerVar[tw].reconClusters[unp] = 0;
  }
  m_trackTowerList.clear();
   
  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon();
  assert(tkrRecon != 0);
 
  int lastTower = -1;
  int numRecCls = 0;
  
  TkrTrack* tkrTrack = m_track;
  TIter trk1HitsItr(tkrTrack);
  TkrTrackHit* pTrk1Hit = 0;
  while( (pTrk1Hit = (TkrTrackHit*)trk1HitsItr.Next()) ) {    
    const TkrCluster* cluster = (pTrk1Hit->getClusterPtr());
    if(!cluster) continue;
    numRecCls++;
    layerId lid = getLayerId( cluster );
    //std::cout << lid.tower << " " << lid.uniPlane << " " << lid.layer << " " << lid.view << std::endl;
    // check if this cluster belong to good track segment.
    if( lid.tower < 0 ){
      m_numErrors++;
      continue;
    }
    if( !m_goodTrackTowerFlag[lid.tower] ) continue;
    int tw = m_towerPtr[ lid.tower ];
    if( tw < 0 ){
      m_numErrors++;
      continue;
    }
    m_towerVar[tw].reconClusters[lid.uniPlane] = cluster;
    if( lid.tower != lastTower ){
      lastTower = lid.tower;
      m_trackTowerList.push_back( lastTower );
      if( lid.view == 0 )
        m_towerVar[tw].center[1] = (cluster->getPosition()).Y();
      else
        m_towerVar[tw].center[0] = (cluster->getPosition()).X();
    }
    m_lrec->Fill( lid.uniPlane );
  }
#ifdef PRINT_DEBUG
  std::cout << "getReconClusters end" << std::endl;
#endif
}

void TkrHits::selectGoodClusters(){
  
#ifdef PRINT_DEBUG
  std::cout << "selectGoodClusters start" << std::endl;
#endif

  bool display = false;
  m_clusters.clear();
  //
  // register new raw clusters if it is close to the track position
  //
  Cluster* cluster;
  for( UInt_t tw=0; tw<m_trackTowerList.size(); tw++ ){
    int tower = m_trackTowerList[tw]; // order of towers for a track
    int twr = m_towerPtr[tower];
    for (int unp=g_nUniPlane-1; unp>=0; unp--){
      layerId lid( unp ); 
      int layer = lid.layer;
      int view = lid.view;
      float zpos = posZ[view][layer];
      //std::cout << layer << " " << view << " " 
      //        << m_towerVar[twr].digiClusters[unp].size() << std::endl;

      int numCls = 0;
      for( UInt_t i=0; i!= m_towerVar[twr].digiClusters[unp].size(); i++){
        cluster = &( m_towerVar[twr].digiClusters[unp].at(i) );

        if( tower == 0 && layer == 7 && view == 0 )
          if( cluster->getFirstStrip() > 768 )
            m_totT0X7H->Fill( cluster->getRawToT() );
          else m_totT0X7L->Fill( cluster->getRawToT() );
        
        // calculate position
        float pos = ( cluster->getLastStrip() + cluster->getFirstStrip() )/2;
        int ladder = int( pos * 4 / g_nStrip );
        pos = pos - g_nStrip * 0.5 + 0.5;
        pos = pos * stripPitch + ladderGap * (ladder - 1.5 ) + m_towerVar[twr].center[view];
        //std::cout <<  cluster->getLastStrip() << " " 
        //<< cluster->getFirstStrip() << " " << pos << std::endl;
        if( view == 0 )
          cluster->setXYZ( pos, m_towerVar[twr].center[1], zpos );
        else
          cluster->setXYZ( m_towerVar[twr].center[0], pos, zpos );
        cluster->setId( tower, unp );
        if( m_towerVar[twr].reconClusters[unp] )
          cluster->setCorrectedTot( m_peakMIP * m_towerVar[twr].reconClusters[unp]->getMips() );
        
        // check if this cluster is close to the track position
        layerId tlid;

        // find closest recon clusters
        float dzmin=10000, dzmin2=10000, dz;
        int numSkip=0, umin=0, umin2=0, tl, tunp;
        const TkrCluster* tcls;
        for( int dl=1; dl<g_nLayer; dl++){
          for( int dir=-1; dir<2; dir+=2){
            tl = layer + dl*dir;
            if( tl < g_nLayer && tl >=0 ){
              tlid.setLayer( tl, view );
              tunp = tlid.uniPlane;
              tcls = m_towerVar[twr].reconClusters[tunp];
              if( tcls ){
                dz = fabs( zpos - tcls->getPosition().Z() );
                if( dz < dzmin ){
                  dzmin2 = dzmin;
                  umin2 = umin;
                  dzmin = dz;
                  umin = tunp;
                }
                else if( dz < dzmin2 ){
                  dzmin2 = dz;
                  umin2 = tunp;        
                }
                else numSkip++;
              }
            }
            if( numSkip > 1 ) break;
          }
          if( numSkip > 1 ) break;
        }

        if( dzmin2>1000 || dzmin>1000 ) continue;

        TVector3 pos1 = m_towerVar[twr].reconClusters[umin]->getPosition();
        TVector3 pos2 = m_towerVar[twr].reconClusters[umin2]->getPosition();

        float delta;
        if( lid.view == 0 ){
          pos = pos1.X() + ( pos2.X()-pos1.X() ) * ( zpos-pos1.Z() ) / ( pos2.Z() - pos1.Z() );
          delta = cluster->getPosition().X() - pos;
        }
        else{
          pos = pos1.Y() + ( pos2.Y()-pos1.Y() ) * ( zpos-pos1.Z() ) / ( pos2.Z() - pos1.Z() );
          delta = cluster->getPosition().Y() - pos;
        }
        
        if( display ) std::cout << tower << " " << lid.layer << " " << lid.view << ", " << delta << " " << " " << pos << " " << zpos;
        m_armsDist->Fill( delta );
        if( fabs(delta) > 3.0 ){
          if( display ) std::cout << " **************" << std::endl;
          continue;
        }
        if( display ) std::cout << std::endl;
        
        //
        // good cluster, multiple cluster per layer allowed.
        //
        numCls++;
        m_clusters.push_back( cluster ); 
        for(int iStrip = cluster->getFirstStrip(); 
            iStrip != int(cluster->getLastStrip()+1); ++iStrip){
          m_towerVar[twr].rHits[unp][iStrip]++;
          m_lcls->Fill( unp );
        }
        if( tower == 0 && layer == 7 && view == 0 )
          if( cluster->getFirstStrip() > 768 )
            m_totT0X7TH->Fill( cluster->getRawToT() );
          else m_totT0X7TL->Fill( cluster->getRawToT() );
      }
      m_numClsDist->Fill( numCls );
      //std::cout << tower << " " << uniPlane << std::endl;
    }
  }

#ifdef PRINT_DEBUG
  std::cout << "selectGoodClusters end" << std::endl;
#endif

}

bool TkrHits::closeToTrack( const TkrCluster* cluster, TkrCluster* clusters[g_nTower][g_nUniPlane] )
{
  layerId lid = getLayerId( cluster ), tlid;
  int tower = lid.tower;
  int layer = lid.layer;
  int view = lid.view;
  float zpos = cluster->getPosition().Z();

  bool display = false;
  //if( layer == 4 && view == 0 ) display = true;
  if( display ) std::cout << "X4 ";

  // find closest hits
  float dzmin=10000, dzmin2=10000, dz;
  int numSkip=0, umin=0, umin2=0, tl, tunp;
  TkrCluster* tcls;
  for( int dl=1; dl<g_nLayer; dl++){
    for( int dir=-1; dir<2; dir+=2){
      tl = layer + dl*dir;
      if( tl < g_nLayer && tl >=0 ){
        tlid.setLayer( tl, view );
        tunp = tlid.uniPlane;
        tcls = clusters[tower][tunp];
        if( tcls ){
          dz = fabs( zpos - tcls->getPosition().Z() );
          if( dz < dzmin ){
            dzmin2 = dzmin;
            umin2 = umin;
            dzmin = dz;
            umin = tunp;
          }
          else if( dz < dzmin2 ){
            dzmin2 = dz;
            umin2 = tunp;        
          }
          else numSkip++;
        }
      }
      if( numSkip > 1 ) break;
    }
    if( numSkip > 1 ) break;
  }

  if( display ) std::cout << numSkip << " " << lid.uniPlane << " " 
                          << umin << " " << dzmin << " " 
                          << umin2 << " " << dzmin2 << ", ";
  
  if( dzmin2>1000 || dzmin>1000 ){
    if( display ) std::cout << std::endl;
    return false;
  }

  TVector3 pos1 = clusters[tower][umin]->getPosition();
  TVector3 pos2 = clusters[tower][umin2]->getPosition();

  float delta, pos;
  if( lid.view == 0 ){
    pos = pos1.X() + ( pos2.X()-pos1.X() ) * ( zpos-pos1.Z() ) / ( pos2.Z() - pos1.Z() );
    delta = cluster->getPosition().X() - pos;
  }
  else{
    pos = pos1.Y() + ( pos2.Y()-pos1.Y() ) * ( zpos-pos1.Z() ) / ( pos2.Z() - pos1.Z() );
    delta = cluster->getPosition().Y() - pos;
  }
  
  if( display ) std::cout << tower << " " << lid.layer << " " << lid.view << ", " << delta << " " << " " << pos << " " << zpos;
  m_armsDist->Fill( delta );
  if( fabs(delta) > 3.0 ){
    if( display ) std::cout << " **************" << std::endl;
    return false;
  }
  if( display ) std::cout << std::endl;

  return true;

}


bool TkrHits::passCut() 
{
  TkrRecon* tkrRecon = m_reconEvent->getTkrRecon(); 
  assert(tkrRecon != 0);
  
  TObjArray* tracks = tkrRecon->getTrackCol();
  int numTracks = tracks->GetEntries();
  m_nTrackDist->Fill( numTracks );
  m_trackRMS = -1.0;
  
  // select only 1 or 2 track event
  if( numTracks > 2 ) return false;
#ifdef DEBUG_PRINT
  std::cout << "# of tracks: " << tracks->GetEntries() << std::endl;
#endif
  
  // find a track with maximum number of hits.
  int maxHits = 0, nHits;
  for( int tk=0; tk!=numTracks; tk++){
    TkrTrack* track = dynamic_cast<TkrTrack*>(tracks->At(tk));
    if(track) {
      nHits = track->getNumFitHits();
      if( nHits > maxHits ){
        maxHits = nHits;
        m_track = track;
        m_pos = track->getInitialPosition();
        m_dir = track->getInitialDirection();
      }
    }
  }
  m_maxHitDist->Fill( maxHits );
  if( maxHits < 6 ) return false;
  //m_trkRMS->Fill( m_track->getScatter() );
  m_trkRMS->Fill( getTrackRMS() );
  //std::cout << m_trackRMS << " " << m_track->getScatter()/m_trackRMS << std::endl;
  m_dirzDist->Fill( m_dir.Z() );
  float maxDirZ = m_maxDirZ;
  if( m_badStrips ) maxDirZ = -0.7;
  if( m_dir.Z() > maxDirZ ) return false;
  //if( m_track->getScatter() > 2.0E-4 ) return false;
  if( m_trackRMS > m_maxTrackRMS ) return false;

  //
  // check MIP filter variable after track cuts
  //
  m_HAcdTileCount->Fill( m_acdTileCount );
  m_HAcdTotalEnergy->Fill( m_acdTotalEnergy );
  m_HCalEnergyRaw->Fill( m_calEnergyRaw );
  m_HNumCalXtal->Fill( m_numCalXtal );
  m_deltaWindOpenTime->Fill( m_digiEvent->getGem().getDeltaWindowOpenTime() );

  return true;
}


float TkrHits::getTrackRMS(){
  if( m_trackRMS >= 0.0 ) return m_trackRMS;

#ifdef DEBUG_PRINT
  std::cout << "getTrackRMS started." << std::endl;
#endif
  //
  // register hits and perform linear fit
  //
  Int_t towerBits = 0;
  std::vector<Double_t> vx[g_nTower], vzx[g_nTower], 
    vy[g_nTower], vzy[g_nTower];
  std::vector<Int_t> upx[g_nTower], upy[g_nTower];
  float offsetx[g_nTower], offsety[g_nTower], trms[g_nTower];
  TkrTrack* tkrTrack = m_track;
  TIter trk1HitsItr(tkrTrack);
  TkrTrackHit* pTrk1Hit = 0;
  while( (pTrk1Hit = (TkrTrackHit*)trk1HitsItr.Next()) ) {    
    const TkrCluster* cluster = (pTrk1Hit->getClusterPtr());
    if(!cluster) continue;
    layerId lid = getLayerId( cluster );
    int twr = lid.tower;
    if( twr < 0 ){
      m_numErrors++;
      continue;
    }
    if( m_towerPtr[lid.tower] < 0 ){
      m_numErrors++;
      continue;
    }
    towerBits |= (1<<lid.tower);
    Double_t z = (cluster->getPosition()).Z();
    if( lid.view == 0 ){
      offsety[twr] = (cluster->getPosition()).Y();
      vx[twr].push_back( (cluster->getPosition()).X() );
      vzx[twr].push_back( z );
      upx[twr].push_back( lid.uniPlane );
    }
    else{
      offsetx[twr] = (cluster->getPosition()).X();
      vy[twr].push_back( (cluster->getPosition()).Y() );
      vzy[twr].push_back( z );
      upy[twr].push_back( lid.uniPlane );
    }
  }
  //
  // check tower with hits and do linear fit
  //
  float sum=0.0, sumsq=0.0, num=0.0;
  Double_t x0[g_nTower], dxdz[g_nTower], y0[g_nTower], dydz[g_nTower];
  Int_t ref1twr=-1, ref2twr=-1;
  UInt_t ref1num=0, ref2num=0, numTracks=0;
  m_trackTowerList.clear();
  for( int twr=0; twr<g_nTower; twr++){
    //if( towerBits&(1<<twr) ) 
    //only tower with sufficient number of hits is used.
    if( ! m_goodTrackTowerFlag[ twr ] ) continue;
    if( vx[twr].size() < 3 || vy[twr].size() < 3 ) continue;
    m_trackTowerList.push_back( twr );
    numTracks++;
    UInt_t numhits = vx[twr].size() + vy[twr].size();
    if( numhits >  ref1num ){
      ref2num = ref1num;
      ref2twr = ref1twr;
      ref1num = numhits;
      ref1twr = twr;      
    }
    else if( numhits > ref2num ){
      ref2num = numhits;
      ref2twr = twr;      
    }

    // fit xz and yz
    leastSquareLinearFit( vx[twr], vzx[twr], x0[twr], dxdz[twr] );
    leastSquareLinearFit( vy[twr], vzy[twr], y0[twr], dydz[twr] );

    // calculare combined rms in both x and y
    float numt=0, sumt=0, sumsqt=0;
    for( UInt_t i=0; i<vx[twr].size(); i++){
      float del = vx[twr][i] - ( x0[twr]+dxdz[twr]*vzx[twr][i] );
      num++;
      sum += del;
      sumsq += del*del;
      numt++;
      sumt += del;
      sumsqt += del*del;
    }
    for( UInt_t i=0; i<vy[twr].size(); i++){
      float del = vy[twr][i] - ( y0[twr]+dydz[twr]*vzy[twr][i] );
      num++;
      sum += del;
      sumsq += del*del;
      numt++;
      sumt += del;
      sumsqt += del*del;
    }
    float mean = sumt / numt;
    float rmsq = sumsqt / numt - mean*mean;
    if( rmsq > 0.0 ) trms[twr] = sqrt( rmsq );
    else trms[twr] = 0.0;
    if( trms[twr] > m_maxTrackRMS ) 
      m_goodTrackTowerFlag[ twr ] = false;
  }
  float mean = sum / num;
  float rmsq = sumsq / num - mean*mean;
  if( rmsq > 0.0 ) m_trackRMS = sqrt( rmsq );
  else m_trackRMS = 0.0;
  if( numTracks == 0 ) return m_trackRMS;

#ifdef DEBUG_PRINT
  std::cout << "getTrackRMS, linear fit finished." << std::endl;
#endif

  for( UInt_t itwr=0; itwr<m_trackTowerList.size(); itwr++){
    int twr = m_trackTowerList[itwr];
    if( m_trackTowerList.size() == 1 ){
      m_trkRMS1TWR->Fill( m_trackRMS );
      if( m_trackRMS < 2.0 ) 
        m_rmsProf1TWR->Fill( twr, m_trackRMS );
    }
    else if( m_trackTowerList.size() == 2 ){
      m_trkRMS2TWR->Fill( m_trackRMS );
      if( trms[twr] < 2.0 )
        m_rmsProf2TWR->Fill( twr, trms[twr] );
    }
  }

  // reset mpos
  Double_t z = m_pos.Z();
  m_pos.SetXYZ( x0[ref1twr] + dxdz[ref1twr]*z, y0[ref1twr] + dydz[ref1twr]*z, z );
  // reset dir
  Double_t dirz = -1 / sqrt( dxdz[ref1twr]*dxdz[ref1twr] + dydz[ref1twr]*dydz[ref1twr] + 1.0 );
  m_dir.SetXYZ( dxdz[ref1twr]*dirz, dydz[ref1twr]*dirz, dirz );

  if( m_trackRMS > 0.4 ) return m_trackRMS;
  if( numTracks > 2 ) return m_trackRMS;
  if( numTracks==2 && ( ref1twr<0 || ref2twr<0 ) ){ // invalid
    std::cout << "invalid reference tower numbers: " << ref1twr 
              << " " << ref2twr << std::endl;
    if( m_log.is_open() )
      m_log << "invalid reference tower numbers: " << ref1twr 
            << " " << ref2twr << std::endl;
    return m_trackRMS;  
  }
  //
  // look at residual for each tower and layer.
  //
#ifdef DEBUG_PRINT
  std::cout << "getTrackRMS, start residual calculation." << std::endl;
#endif
  float trkpos, residual, boundx, boundy, boundz, dldz, dist, trad;
  const float maxrms = 0.3, maxtrad=6.0, mindist=10.0, maxdist=450;
  for( UInt_t itwr=0; itwr!=m_trackTowerList.size(); itwr++){
    int twr = m_trackTowerList[itwr];
    if( trms[twr] > maxrms ) continue;
    int reftwr = -1;
    if( numTracks==2 ){ 
      // determine reference tower
      if( twr!=ref1twr ) reftwr = ref1twr;
      else if( twr!=ref2twr ) reftwr = ref2twr;
      else{
        std::cout << "invalid tower numbers: " << twr << " " << ref1twr 
                  << " " << ref2twr << std::endl;
        if( m_log.is_open() )
          m_log << "invalid tower numbers: " << twr << " " << ref1twr 
                << " " << ref2twr << std::endl;
        continue;
      }
#ifdef DEBUG_PRINT
      std::cout << "getTrackRMS: " << twr << " " << reftwr
                << ", " << offsetx[twr] << " " << offsety[twr] 
                << ", " << offsetx[reftwr] << " " << offsety[reftwr] 
                << ", " << x0[reftwr] << " " << y0[reftwr]
                << ", " << dxdz[reftwr] << " " << dydz[reftwr]
                << std::endl;
#endif
      //
      // calculate track intercept at tower boudary
      if( offsetx[twr] == offsetx[reftwr] ){ // boundary at y
        boundy = 0.5 * ( offsety[twr] + offsety[reftwr] );
        boundz = (boundy-y0[reftwr]) / dydz[reftwr];
        boundx = boundz * dxdz[reftwr] + x0[reftwr];
        // total thickness of wall
        trad = sqrt( 1 + (1+dxdz[reftwr]*dxdz[reftwr])/(dydz[reftwr]*dydz[reftwr]) );
      }
      else if( offsety[twr] == offsety[reftwr] ){ // boundary at x
        boundx = 0.5 * ( offsetx[twr] + offsetx[reftwr] );
        boundz = (boundx-x0[reftwr]) / dxdz[reftwr];
        boundy = boundz * dydz[reftwr] + y0[reftwr];
        trad = sqrt( 1 + (1+dydz[reftwr]*dydz[reftwr])/(dxdz[reftwr]*dxdz[reftwr]) );
      }
      else numTracks = 3;
      dldz = sqrt( dxdz[reftwr]*dxdz[reftwr] + dydz[reftwr]*dydz[reftwr] + 1 );
    }
    //
    // residual for x
    for( UInt_t i=0; i<vx[twr].size(); i++){
      residual = vx[twr][i] - ( x0[twr]+dxdz[twr]*vzx[twr][i] );
      if( fabs( residual ) > 1.0 ) continue;
      if( m_MIPtot && dirz<m_maxDirZ ) 
        m_towerVar[m_towerPtr[twr]].resProf->Fill( upx[twr][i], residual );
      if( numTracks==2 ){
        trkpos = x0[reftwr] + dxdz[reftwr]*vzx[twr][i];
        if( fabs( trkpos-boundx ) < 5.0 ) continue;
        residual = vx[twr][i] - trkpos;
        dist = dldz * fabs(vzx[twr][i]-boundz);
        if( trms[reftwr]<0.3 && trad<maxtrad )
          m_sigDist->Fill( dist, residual );
        if( trad<maxtrad && dist>mindist && dist<maxdist)
          m_sigRMS->Fill( trms[reftwr], residual );
        if( trms[reftwr]<0.3 && dist>mindist && dist<maxdist)
          m_sigTrad->Fill( trad, residual );
        if( fabs(residual)<3.0 && trms[reftwr]<0.3 
            && trad<maxtrad && dist>mindist && dist<maxdist ){
          m_tresProfX->Fill( twr, residual );
          m_towerVar[m_towerPtr[twr]].resX->Fill( residual );
        }
      }
    }
    //
    // residual for y
    for( UInt_t i=0; i<vy[twr].size(); i++){
      residual = vy[twr][i] - ( y0[twr]+dydz[twr]*vzy[twr][i] );
      if( twr==0 && upy[twr][i]==4 && m_MIPtot ){
        m_res->Fill( residual );
        if( dirz<m_maxDirZ ) m_resSel->Fill( residual );
      }
      if( fabs( residual ) > 1.0 ) continue;
      if( m_MIPtot && dirz<m_maxDirZ ) 
        m_towerVar[m_towerPtr[twr]].resProf->Fill( upy[twr][i], residual );
      if( numTracks==2 ){
        trkpos = y0[reftwr] + dydz[reftwr]*vzy[twr][i];
        if( fabs( trkpos-boundy ) < 5.0 ) continue;
        residual = vy[twr][i] - trkpos;
        dist = dldz * fabs(vzy[twr][i]-boundz);
        if( trms[reftwr]<maxrms && trad<maxtrad )
          m_sigDist->Fill( dist, residual );
        if( trad<maxtrad && dist>mindist && dist<maxdist)
          m_sigRMS->Fill( trms[reftwr], residual );
        if( trms[reftwr]<maxrms && dist>mindist && dist<maxdist)
          m_sigTrad->Fill( trad, residual );
        if( fabs(residual)<3.0 && trms[reftwr]<maxrms 
            && trad<maxtrad && dist>mindist && dist<maxdist ){
          m_tresProfY->Fill( twr, residual );
          m_towerVar[m_towerPtr[twr]].resY->Fill( residual );
        }
      }
    }
  }

#ifdef DEBUG_PRINT
  std::cout << "getTrackRMS finished." << std::endl;
#endif
  return m_trackRMS;

}


Double_t TkrHits::leastSquareLinearFit( std::vector<Double_t> &vy, 
                                        std::vector<Double_t> &vx, 
                                        Double_t &y0, Double_t &dydx ){

  Double_t sumX=0.0, sumXX=0.0, sumXY=0.0, sumY=0.0, sumYY=0.0;
  UInt_t num = vy.size();
  if( num < 3 ){
    y0 = 0;
    dydx = 0;
    return -1.0;
  }

  for( UInt_t i=0; i<num; i++){
    sumX += vx[i];
    sumXX += vx[i]*vx[i];
    sumXY += vx[i]*vy[i];
    sumYY += vy[i]*vy[i];
    sumY += vy[i];
  }
  dydx = ( sumXY*num - sumX*sumY ) / ( sumXX*num - sumX*sumX );
  y0 = ( sumY - dydx*sumX ) / num;
  Double_t rms = sumYY - 2*dydx*sumXY - 2*y0*sumY 
    + dydx*dydx*sumXX + 2*dydx*y0*sumX + y0*y0*num;
  if( rms > 0.0 ) rms = sqrt( rms / num );
  else rms = 0.0;
  //std::cout << rms << std::endl;
  return rms;
}


void TkrHits::fillOccupancy( int tDiv ) 
{
#ifdef PRINT_DEBUG
  std::cout << "fillOccupancy start" << std::endl;
#endif
  
  int lStrip = g_nStrip/g_nLadder;

  //initialize container
  int nHits[g_nTower][g_nLayer][g_nView][g_nWafer+1];
  for( unsigned int tw=0; tw<m_towerVar.size(); tw++){
    int tower = m_towerVar[tw].towerId;
    for( int layer=0; layer<g_nLayer; layer++)
      for( int view=0; view<g_nView; view++)
        for( int i=0; i<g_nWafer+1; i++) nHits[tower][layer][view][i] = 0;
  }
  
  //
  // first loop to register hits and get tower offset
  //
  int hitLayers[g_nLayer];
  int hitPlanes[g_nLayer][g_nView];
  for( int layer=0; layer!=g_nLayer; layer++){
    hitLayers[layer]=0;
    for( int vw=0; vw!=g_nView; vw++) hitPlanes[layer][vw]=0;
  }
  
  for( unsigned int cls=0; cls<m_clusters.size(); cls++){
    Cluster* cluster = m_clusters[cls];
    layerId lid = getLayerId( cluster );
    int tower = lid.tower;
    int view = lid.view;
    int layer = lid.layer;
    
    for(int iStrip = cluster->getFirstStrip(); 
        iStrip != int(cluster->getLastStrip()+1); ++iStrip){
      nHits[tower][layer][view][iStrip/lStrip]++;
      nHits[tower][layer][view][g_nWafer]++;
    }
    
    hitLayers[layer]++;
    hitPlanes[layer][view]++;
  }
  
  //
  // main loop to fill occupancy and track position
  //
  float pos, apos, tpos[g_nView], lpos, posz;
  float dist, dxz, dyz, dz, delta, deltax, deltay;
  float dirX=m_dir.X()/m_dir.Z(), dirY=m_dir.Y()/m_dir.Z(), 
    preX=m_pos.X(), preY=m_pos.Y(), preXZ=m_pos.Z(), preYZ=m_pos.Z();
  int aview, preLayer=g_nLayer;
  int lastTower=-1, nTowers=0, towers[2], strips[g_nView];

  int numCls = m_clusters.size();
  for( int cls=0; cls<numCls; cls++){
    Cluster* cluster = m_clusters[cls];
    layerId lid = getLayerId( cluster );
    int tower = lid.tower;
    int view = lid.view;
    int layer = lid.layer;
    int unp = lid.uniPlane;

    // check if track moves across towers.
    if( tower != lastTower ){
      if( lastTower < 0 ){ 
        lastTower = tower;
        towers[0] = tower;
        nTowers = 1;
      }
      else{
        towers[0] = lastTower;
        towers[1] = tower;
        nTowers = 2;
        lastTower = -1;
      }
    }
    
    // fill track positions in all layer between previous and current layer
    // carefull for moving across towers.
    int elyr = layer;
    if( cls == numCls-1 ) elyr = 0;
    for( int lyr=preLayer-1; lyr>= elyr; lyr--){ 
      // disable following statement since it reduces the tracks in the layer edge.
      //if( nTowers > 1 ) continue; // avoid tower transition region
      // layers where hits are expected.
      // hit in the same layer or hits in both layers below and above.
      if( hitLayers[lyr] != 0
          || ( lyr!=0 && lyr!=g_nLayer-1 
               && hitLayers[lyr+1]>0 && hitLayers[lyr-1]>0 ) ){
        posz =  posZ[view][lyr];
        dxz = posz - preXZ;
        tpos[0] = preX + dirX*dxz;
        dyz = posz - preYZ;
        tpos[1] = preY + dirY*dyz;
        
        int margin = 20;
        for( int tw=0; tw<nTowers; tw++){ // check both tower
          int twr = towers[tw];
          int vtw = m_towerPtr[twr];
          int numActive=0;
          int numVG=0;
          for( int vw=0; vw!=g_nView; vw++){
            strips[vw]=-1;
            lpos = tpos[vw] - m_towerVar[vtw].center[vw];
            for( int iw=0; iw!=g_nWafer; ++iw){
              float stp = ( lpos-ladderGap*(iw-1.5) ) / stripPitch 
                + g_nStrip/2;
              if( stp>=iw*lStrip && stp<(iw+1)*lStrip ) strips[vw] = int(stp);
            }
            if( strips[vw] > 0 ){ 
              numActive++;
              int stp = strips[vw] % lStrip;
              if( stp>margin && stp<lStrip-margin ) numVG++;
            }
          }
          // check if track go thorugh active region in all views
          if( numActive == g_nView ){
            for( int vw=0; vw!=g_nView; vw++){
              layerId lid( lyr, vw );
              if( m_MIPeff ) 
                m_towerVar[vtw].bsVar[lid.uniPlane].tHits[strips[vw]]++;
              if( numVG == g_nView ){
                if( m_MIPeff ) m_towerVar[vtw].bsVar[lid.uniPlane].tLayer++;
                for(int icut=0; icut<ncut; icut++)
                  if( m_cut[icut] ) m_trackCut->Fill( icut+0.5 );
              }
              m_ltrk->Fill( lid.uniPlane );
              // layer with associated hits
              if( hitPlanes[layer][vw]>0 ){
                if( m_MIPeff ) 
                  m_towerVar[vtw].bsVar[lid.uniPlane].eHits[strips[vw]]++;
                if( numVG == g_nView ){
                  if( m_MIPeff ) 
                    m_towerVar[vtw].bsVar[lid.uniPlane].hLayer++;
                  for(int icut=0; icut<ncut; icut++)
                    if( m_cut[icut] ) m_hitCut->Fill( icut+0.5 );
                }
              }
            }
          }
        } // tower loop
      } // valid layer
    } // layer loop
    preLayer = layer;

    TVector3 position = cluster->getPosition();
    posz = position.Z();
    if( fabs(posz-posZ[view][layer]) > 0.01 )
      std::cout << "Incosistent z position: " << layer << " " << view << ", " 
                << posZ[view][layer] << " != " << posz << std::endl;
    dyz = posz - preYZ;
    dxz = posz - preXZ;
    deltax = preX + dirX*dxz - position.X();
    deltay = preY + dirY*dyz - position.Y();
    
    if( view == 0 ){
      dz = dxz;
      delta = deltax;
    }
    else{
      dz = dyz;
      delta = deltay;
    }
    
    float dx = dirX*dz;
    float dy = dirY*dz;
    dist =sqrt( dz*dz+dx*dx+dy*dy );
    m_dist->Fill( dist );
    if( dist < 30 ) dist = 30;
    delta *= (35.0/dist);
    m_brmsDist[layer/3]->Fill( delta );
    //if( layer==4 && view==0 ) m_brmsDist[layer/3]->Fill( delta );
    
    // select good clusters
    if( (fabs(delta) > 3.0) || !m_MIPtot ) continue;
    
    if( view == 0 ){
      aview = 1;
      pos = deltax;
      apos = deltay;
      if( dxz > 10.0 ) dirX = ( position.X() - preX ) / dxz;
      preX = position.X();
      preXZ = position.Z();
    }
    else{
      aview = 0;
      pos = deltay;
      apos = deltax;
      if( dyz > 10.0 ) dirY = ( position.Y() - preY ) / dyz;
      preY = position.Y();
      preYZ = position.Z();
    }
    
    //std::cout << layer << " " << view << ", " << pos << " " << apos
    //      << std::endl;
    
    int twr = m_towerPtr[tower];
    m_locc->Fill( lid.uniPlane );
    for(int iStrip = cluster->getFirstStrip(); 
        iStrip != int(cluster->getLastStrip()+1); ++iStrip){
      m_towerVar[twr].bsVar[unp].lHits[iStrip]++;
      if( nHits[tower][layer][aview][g_nWafer] > 0 ){
        for( int iw=0; iw<g_nWafer; iw++ )
          if( nHits[tower][layer][aview][iw] > 0 ){
            m_towerVar[twr].bsVar[unp].nHits[iStrip][iw][tDiv]++;
            m_aPos[iw]->Fill( apos-89.5*(iw-1.5) );
          }
      }
      else{
        for( int iw=0; iw<g_nWafer; iw++ )
          if( fabs( apos-89.5*(iw-1.5) ) < 44.4 ){
            m_towerVar[twr].bsVar[unp].nHits[iStrip][iw][tDiv]++;
            m_aPos[iw]->Fill( apos-89.5*(iw-1.5) );
          }
      }
    }
  }
  
#ifdef PRINT_DEBUG
  std::cout << "fillOccupancy end" << std::endl;
#endif
}
