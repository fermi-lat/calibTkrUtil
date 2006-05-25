#include "TkrNoiseOcc.h"

TkrNoiseOcc::TkrNoiseOcc() 
{}

TkrNoiseOcc::~TkrNoiseOcc() 
{}

void 
TkrNoiseOcc::setCritStripRate(float crit_strip_rate){
  m_crit_strip_rate = crit_strip_rate;
}

void 
TkrNoiseOcc::setCritLayerRate(float crit_layer_rate){
  m_crit_layer_rate = crit_layer_rate;
}

void 
TkrNoiseOcc::setTrigCut(int trig_cut){
  m_trig_cut = trig_cut;
}

/*
void
TkrNoiseOcc::initPar(int coincidence_cut, int multi_ld, int multi_hd) 
{
  m_coincidence_cut = coincidence_cut;
  m_multi_ld = multi_ld;
  m_multi_hd = multi_hd;
}
*/

void 
TkrNoiseOcc::initAnalysis(int nEvent, int evt_interval, int coincidence_cut,
			  int multi_ld, int multi_hd, int periodic_trig) {

  int tower, bilayer, xyview;
  int ix;

  m_nEvent          = nEvent;	   
  m_evt_interval    = evt_interval;   
  m_coincidence_cut = coincidence_cut;
  m_multi_ld        = multi_ld;       
  m_multi_hd        = multi_hd;
  m_periodic_trig   = periodic_trig;


  m_crit_strip_rate = 5.0e-5;
  m_crit_layer_rate = 5.0e-2;
  m_trig_cut = 0;

  m_nx = (int)(m_nEvent/m_evt_interval)+1;
  m_event_counter   = 0;
  
  for(tower=0;tower<g_nTower; tower++){
    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      vTkrExposure[tower][bilayer] = new float[m_nx];
      for(ix=0; ix<m_nx; ix++) vTkrExposure[tower][bilayer][ix]=0.0;
      
      for(xyview=0; xyview<g_nView; xyview++) {
	vTkrStripOcc[tower][bilayer][xyview]   = new float[m_nx];
	vTkrLayerOcc[tower][bilayer][xyview]   = new float[m_nx];
	vTkrHitMap[tower][bilayer][xyview]     = new float[g_nStripsPerLayer];
	vTkrNoiseMul[tower][bilayer][xyview]   = new float[MUL_MAX];
	vTkrNoiseTot0[tower][bilayer][xyview]  = new float[TOT_MAX];
	vTkrNoiseTot1[tower][bilayer][xyview]  = new float[TOT_MAX];

	for(ix=0; ix<m_nx; ix++)   vTkrStripOcc[tower][bilayer][xyview][ix]  =0.0;
	for(ix=0; ix<m_nx; ix++)   vTkrLayerOcc[tower][bilayer][xyview][ix]  =0.0;
	for(ix=0; ix<g_nStripsPerLayer; ix++) vTkrHitMap[tower][bilayer][xyview][ix] =0.0;
	for(ix=0; ix<150; ix++)  vTkrNoiseMul[tower][bilayer][xyview][ix]  =0.0;
	for(ix=0; ix<300; ix++)  vTkrNoiseTot0[tower][bilayer][xyview][ix] =0.0;
	for(ix=0; ix<300; ix++)  vTkrNoiseTot1[tower][bilayer][xyview][ix] =0.0;
      }
    }
  }
  
}

void
TkrNoiseOcc::setDigiEvtPtr(DigiEvent *digiEvt) {
  m_digiEvt = digiEvt;
}

void
TkrNoiseOcc::anaDigiEvt() {

  int tower, bilayer, xyview;
  int err_cout = 0;
  int ievent, ihit, val;
  int numHitLayer[g_nTower][g_nTkrLayer][g_nView];
  int tot0[g_nTower][g_nTkrLayer][g_nView], tot1[g_nTower][g_nTkrLayer][g_nView];
  int buf_stripId[g_nTower][g_nTkrLayer][g_nView][128];

  UShort_t gemCondSummary, gemTkrVector;


  const Gem& gem = m_digiEvt->getGem();
  gemCondSummary = gem.getConditionSummary();
  gemTkrVector   = gem.getTkrVector();

  if (m_periodic_trig==1) {
    //Gem
    if ((gemCondSummary&0x20)==0x0) {
      m_event_counter +=1;
      return;
    }
    //std::cout << "m_event_counter = " << m_event_counter << std::endl;
  }

  for(tower=0; tower<g_nTower; tower++) {
    for(bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      for(xyview=0; xyview<g_nView; xyview++) {
	numHitLayer[tower][bilayer][xyview]=0;
	tot0[tower][bilayer][xyview]=TOT_INI_VAL;
	tot1[tower][bilayer][xyview]=TOT_INI_VAL;
      }
    }
  }

  ievent = m_digiEvt->getEventId(); 

  ///////////////////////
  // Tkr Digi Analysis
  ///////////////////////
  
  // The full collection of TkrDigis for this event
  const TObjArray* tkrDigiCol = m_digiEvt->getTkrDigiCol();
  if (!tkrDigiCol) return;
      
  // Total number of TkrDigis for this event
  //Int_t numTkrDigi = tkrDigiCol->GetEntries();

  // Loop over all TkrDigis
  TIter tkrIter(tkrDigiCol);
  TkrDigi *tkr = 0;
  while ((tkr = (TkrDigi*)tkrIter.Next())) {

    // Identify the tower and layer
    tower = tkr->getTower().id();
    bilayer = tkr->getBilayer();
      
    // Returns the orientation of the strips
    GlastAxis::axis view = tkr->getView();
    
    numHitLayer[tower][bilayer][view] = tkr->getNumHits();
    tot0[tower][bilayer][view] = tkr->getToT(0);
    tot1[tower][bilayer][view] = tkr->getToT(1);
    
    // Loop through collection of hit strips for this TkrDigi
    for (ihit = 0; ihit < numHitLayer[tower][bilayer][view]; ihit++) {
      // Retrieve the strip number
      buf_stripId[tower][bilayer][view][ihit] = tkr->getStrip(ihit);
    }
  }  //while(tkr)


  for(tower=0; tower<g_nTower; tower++){
    for(bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      
      /// check TKR trigger
      if ( (m_trig_cut==1) && ( ((gemTkrVector>>tower)&0x1)==1 ) ) {
	continue;
      }
      /// check track hits on the adjacent layers
      if (m_coincidence_cut>0) {
	if ( (numHitLayer[tower][bilayer][0]>0) && (numHitLayer[tower][bilayer][1]>0) ) {
	  continue;
	}
	if (bilayer>0) {
	  if ( (numHitLayer[tower][bilayer-1][0]>0) || (numHitLayer[tower][bilayer-1][1]>0) ) {
	    continue;
	  }
	} 
	if ( bilayer<(g_nTkrLayer-1) ){
	  if ( (numHitLayer[tower][bilayer+1][0]>0) || (numHitLayer[tower][bilayer+1][1]>0) ) {
	    continue;
	  }
	}
      }
      
      // Fill Exposure
      vTkrExposure[tower][bilayer][(int)(m_event_counter/m_evt_interval)] +=1.0;
      
      for(xyview=0; xyview<g_nView; xyview++){ 
	
	// Cut by Multiplicity
	if( numHitLayer[tower][bilayer][xyview]<m_multi_ld) {
	  continue;
	}
	if( numHitLayer[tower][bilayer][xyview]>m_multi_hd) {
	  continue;
	}
	
	// Fill Strip Occupancy
	vTkrStripOcc[tower][bilayer][xyview][(int)(m_event_counter/m_evt_interval)] += (float)numHitLayer[tower][bilayer][xyview];
	
	// Fill Event Occupancy
	if (numHitLayer[tower][bilayer][xyview]>0) {
	  vTkrLayerOcc[tower][bilayer][xyview][(int)(m_event_counter/m_evt_interval)] +=1.0;
	}
	// Fill Noise Multiplicity
	vTkrNoiseMul[tower][bilayer][xyview][numHitLayer[tower][bilayer][xyview]] +=1.0;
	// Fill HitMap
	for(ihit=0; ihit<numHitLayer[tower][bilayer][xyview]; ihit++) {
	  vTkrHitMap[tower][bilayer][xyview][buf_stripId[tower][bilayer][xyview][ihit]] +=1.0;
	}
	//Fill Tot
	val = tot0[tower][bilayer][xyview];
	if ( val<0 || 299<val ) val = 299;
	vTkrNoiseTot0[tower][bilayer][xyview][val] +=1.0;
	
	if ((err_cout==1)&&(val==299)) {
	  std::cout << " ievent = " << ievent;
	  std::cout << " tower = " << tower;
	  std::cout << " layer = " << bilayer*g_nView+xyview;
	  std::cout << " tot0 = "  << tot0[tower][bilayer][xyview] << std::endl;
	}	      

	val = tot1[tower][bilayer][xyview];
	if ( val<0 || 299<val ) val = 299;
	vTkrNoiseTot1[tower][bilayer][xyview][val] +=1.0;
	
	if ((err_cout==1)&&(val==299)) {
	  std::cout << " ievent = " << ievent;
	  std::cout << " tower = " << tower;
	  std::cout << " layer = " << bilayer*g_nView+xyview;
	  std::cout << " tot1 = "  << tot1[tower][bilayer][xyview] << std::endl;
	}	      
	
      }
    }
  }
  m_event_counter +=1;
}  


void 
TkrNoiseOcc::clearAnalysis() {

  int tower, bilayer, xyview;
  //int m_nx;

  m_event_counter = 0;
  /// delete histgram array
  for(tower=0;tower<g_nTower; tower++){
    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      delete[] vTkrExposure[tower][bilayer];
      vTkrExposure[tower][bilayer] = NULL;
      
      for(xyview=0; xyview<g_nView; xyview++) {
	delete[] vTkrStripOcc[tower][bilayer][xyview]; 
	delete[] vTkrLayerOcc[tower][bilayer][xyview]; 
	delete[] vTkrHitMap[tower][bilayer][xyview];   
	delete[] vTkrNoiseMul[tower][bilayer][xyview]; 
	delete[] vTkrNoiseTot0[tower][bilayer][xyview];
	delete[] vTkrNoiseTot1[tower][bilayer][xyview];
	
	vTkrStripOcc[tower][bilayer][xyview]   = NULL;
	vTkrLayerOcc[tower][bilayer][xyview]   = NULL;
	vTkrHitMap[tower][bilayer][xyview]     = NULL;
	vTkrNoiseMul[tower][bilayer][xyview]   = NULL;
	vTkrNoiseTot0[tower][bilayer][xyview]  = NULL;
	vTkrNoiseTot1[tower][bilayer][xyview]  = NULL;
      }
    }
  }
}


void TkrNoiseOcc::openHistFile(char* histFileName){
  m_histFile = new TFile(histFileName, "RECREATE");
}

void TkrNoiseOcc::closeHistFile(){
  if(m_histFile==NULL) {
    std::cout << "No output file is opened" << std::endl;
    return;
  }
  m_histFile->Close(); 
}

void
TkrNoiseOcc::saveAnaToHis(char* histFileName){

  openHistFile(histFileName);
  if(m_histFile==NULL) {
    std::cout << "No output file is opened" << std::endl;
    return;
  }

  TDirectory* tkrNoiseOcc_dir = m_histFile->mkdir("TkrNoiseOcc");
  writeAnaToHis(tkrNoiseOcc_dir);
  closeHistFile();
}

void 
TkrNoiseOcc::writeAnaToHis(TDirectory* dirTkrNoise){

  int tower, bilayer, xyview;
  int ix;

  float max_rate, ave_rate;
  float xmi, xma;

  if(dirTkrNoise==NULL) {
    std::cout << "No output directory is opened" << std::endl;
    return;
  }

  dirTkrNoise->cd();

  xmi = 0.5;
  xma = (float)(m_nx*m_evt_interval)+0.5;

  char hname[80], htitle[80], dirname[80];

  ///  histgram definitions
  //TDirectory *dirTower;
  TDirectory *dirExposure, *dirStripOcc, *dirLayerOcc, *dirStripRaw, *dirLayerRaw, *dirHitmap, *dirMulti, *dirTOT;
  TDirectory *dirExposureRt, *dirStripOccRt, *dirLayerOccRt, *dirStripRawRt, *dirLayerRawRt, *dirHitmapRt, *dirMultiRt, *dirTOTRt;

  dirExposureRt  = dirTkrNoise->mkdir("Exposure");
  dirStripRawRt  = dirTkrNoise->mkdir("StripRaw");
  dirLayerRawRt  = dirTkrNoise->mkdir("LayerRaw");
  dirStripOccRt  = dirTkrNoise->mkdir("StripOcc");
  dirLayerOccRt  = dirTkrNoise->mkdir("LayerOcc");
  dirHitmapRt    = dirTkrNoise->mkdir("Hitmap");
  dirMultiRt     = dirTkrNoise->mkdir("Multi");
  dirTOTRt       = dirTkrNoise->mkdir("TOT");

  TH1F *hTkrExposure[g_nTower][g_nTkrLayer];
  TH1F *hTkrStripOcc[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrLayerOcc[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrHitMap[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrNoiseMul[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrNoiseTot0[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrNoiseTot1[g_nTower][g_nTkrLayer][g_nView];


  for(tower=0;tower<g_nTower; tower++){

    sprintf(dirname, "Tower%d", tower);
    //dirTower = dirTkrNoise->mkdir(dirname);
    dirExposure  = dirExposureRt->mkdir(dirname);
    dirStripRaw  = dirStripRawRt->mkdir(dirname);
    dirLayerRaw  = dirLayerRawRt->mkdir(dirname);
    dirHitmap    = dirHitmapRt->mkdir(dirname);
    dirMulti     = dirMultiRt->mkdir(dirname);
    dirTOT       = dirTOTRt->mkdir(dirname);
    
    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {

      dirExposure->cd();
      sprintf(hname, "hTkrExposTwr%dbiLayer%d", tower, bilayer);
      sprintf(htitle, " TKR  Tower%d biLayer%d Exposure", tower, bilayer);
      hTkrExposure[tower][bilayer] = new TH1F(hname, htitle, m_nx, xmi, xma);
      for (ix=0; ix<m_nx; ix++) hTkrExposure[tower][bilayer]->SetBinContent(ix+1, vTkrExposure[tower][bilayer][ix]);
      
      for(xyview=0; xyview<g_nView; xyview++) {

	// StripOcc
	dirStripRaw->cd();
	sprintf(hname, "hTkrStripOccTwr%dLayer%d", tower, bilayer*g_nView+xyview);
	sprintf(htitle, "Tower%d-TKR Averaged Hit-Strip Occupancy per Layer%d", tower, bilayer*g_nView+xyview);
	hTkrStripOcc[tower][bilayer][xyview]  = new TH1F(hname, htitle, m_nx, xmi, xma);
	for (ix=0; ix<m_nx; ix++) hTkrStripOcc[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrStripOcc[tower][bilayer][xyview][ix]);
	hTkrStripOcc[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	// LayerOcc
	dirLayerRaw->cd();
	sprintf(hname, "hTkrLayerOccTwr%dLayer%d", tower, bilayer*g_nView+xyview);
	sprintf(htitle, "Tower%d-TKR Layer Occupancy per Layer%d", tower, bilayer*g_nView+xyview);
	hTkrLayerOcc[tower][bilayer][xyview]  = new TH1F(hname, htitle, m_nx, xmi, xma);
	for (ix=0; ix<m_nx; ix++) hTkrLayerOcc[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrLayerOcc[tower][bilayer][xyview][ix]);
	hTkrLayerOcc[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	// Hitmap
	dirHitmap->cd();
	sprintf(hname, "hTkrHitMapTwr%dLayer%d", tower, bilayer*g_nView+xyview);
	sprintf(htitle, "Tower%d-TKR Layer%d HitMap", tower, bilayer*g_nView+xyview);
	hTkrHitMap[tower][bilayer][xyview]  = new TH1F(hname, htitle, g_nStripsPerLayer, -0.5, g_nStripsPerLayer-0.5);
	for (ix=0; ix<g_nStripsPerLayer; ix++) hTkrHitMap[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrHitMap[tower][bilayer][xyview][ix]);
	hTkrHitMap[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	// Multiplicity
	dirMulti->cd();
	sprintf(hname, "hTkrNoiseMulTwr%dLayer%d", tower, bilayer*g_nView+xyview);
	sprintf(htitle, "Tower%d-TKR Layer%d Noise Multiplicity", tower, bilayer*g_nView+xyview);
	hTkrNoiseMul[tower][bilayer][xyview]  = new TH1F(hname, htitle, 150, -0.5, 149.5);
	for (ix=0; ix<150; ix++) hTkrNoiseMul[tower][bilayer][xyview]-> SetBinContent(ix+1, vTkrNoiseMul[tower][bilayer][xyview][ix]);
	hTkrNoiseMul[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	// TOT
	dirTOT->cd();
	sprintf(hname, "hTkrNoiseTot0Twr%dLayer%d", tower, bilayer*g_nView+xyview);
	sprintf(htitle, "Tower%d-TKR Layer%d ToT0", tower, bilayer*g_nView+xyview);
	hTkrNoiseTot0[tower][bilayer][xyview]  = new TH1F(hname, htitle, 300, -0.5, 299.5);
	for (ix=0; ix<300; ix++) hTkrNoiseTot0[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrNoiseTot0[tower][bilayer][xyview][ix]); 
	hTkrNoiseTot0[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	sprintf(hname, "hTkrNoiseTot1Twr%dLayer%d", tower, bilayer*g_nView+xyview);
	sprintf(htitle, "Tower%d-TKR Layer%d ToT1", tower, bilayer*g_nView+xyview);
	hTkrNoiseTot1[tower][bilayer][xyview]  = new TH1F(hname, htitle, 300, -0.5, 299.5);
	for (ix=0; ix<300; ix++) hTkrNoiseTot1[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrNoiseTot1[tower][bilayer][xyview][ix]);
	hTkrNoiseTot1[tower][bilayer][xyview]->Write("",TObject::kOverwrite);
      }
    }
  }
  
  /// scale histgram
  for (tower=0;tower<g_nTower; tower++){
    
    sprintf(dirname, "Tower%d", tower);
    dirStripOcc = dirStripOccRt->mkdir(dirname);
    dirLayerOcc = dirLayerOccRt->mkdir(dirname);

    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      for (xyview=0; xyview<g_nView; xyview++) { 

	dirStripOcc->cd();
	// scale of Strip Occupancy
	hTkrStripOcc[tower][bilayer][xyview]->Divide(hTkrExposure[tower][bilayer]);
	hTkrStripOcc[tower][bilayer][xyview]->Scale(1.0/g_nStripsPerLayer);
	hTkrStripOcc[tower][bilayer][xyview]->Write(0,TObject::kOverwrite);

	dirLayerOcc->cd();
	// scale of Layer Occupancy
	hTkrLayerOcc[tower][bilayer][xyview]->Divide(hTkrExposure[tower][bilayer]);
	hTkrLayerOcc[tower][bilayer][xyview]->Write(0,TObject::kOverwrite);

	// check maximum strip rate
	max_rate = hTkrStripOcc[tower][bilayer][xyview]->GetMaximum();
	if (max_rate>m_crit_strip_rate) {
	  std::cout << " Tower" << tower;
	  std::cout << " Layer" << bilayer*g_nView+xyview;
	  std::cout << std::scientific << " Max Strip Occupancy=" << max_rate << std::endl;
	}
	// check maximum layer rate
	max_rate = hTkrLayerOcc[tower][bilayer][xyview]->GetMaximum();
	if (max_rate>m_crit_layer_rate) {
	  std::cout << " Tower" << tower;
	  std::cout << " Layer" << bilayer*g_nView+xyview;
	  std::cout << std::scientific << " Max Layer ccupancy=" << max_rate << std::endl;
	}
	
      }
    }
  }
  
  
  /// Summary plot
  dirTkrNoise->cd();
  TH2F *hMaxStripOcc = new TH2F("hMaxStripOcc", "Max Noise Occupancy averaged over layer", g_nTower, -0.5, g_nTower-0.5, 36, -0.5, 35.5);
  TH2F *hMaxLayerOcc = new TH2F("hMaxLayerOcc", "Max Layer Occupancy averaged over layer", g_nTower, -0.5, g_nTower-0.5, 36, -0.5, 35.5);
  TH2F *hAveStripOcc = new TH2F("hAveStripOcc", "Average Noise Occupancy averaged over layer", g_nTower, -0.5, g_nTower-0.5, 36, -0.5, 35.5);
  TH2F *hAveLayerOcc = new TH2F("hAveLayerOcc", "Average Layer Occupancy averaged over layer", g_nTower, -0.5, g_nTower-0.5, 36, -0.5, 35.5);

  //TH1F *hTowerAveOcc = new TH1F("hTowerAveOcc", "Tower Average Noise Occupancy", 16, -0.5, 15.5);

  for(tower=0; tower<g_nTower; tower++){
    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      for (xyview=0; xyview<g_nView; xyview++) { 

	max_rate = hTkrStripOcc[tower][bilayer][xyview]->GetMaximum();
	hMaxStripOcc->Fill(tower, bilayer*g_nView+xyview, max_rate);

	max_rate = hTkrLayerOcc[tower][bilayer][xyview]->GetMaximum();
	hMaxLayerOcc->Fill(tower, bilayer*g_nView+xyview, max_rate);

	ave_rate = hTkrStripOcc[tower][bilayer][xyview]->GetSum()/hTkrStripOcc[tower][bilayer][xyview]->GetNbinsX();
	hAveStripOcc->Fill(tower, bilayer*g_nView+xyview, ave_rate);

	ave_rate = hTkrLayerOcc[tower][bilayer][xyview]->GetSum()/hTkrLayerOcc[tower][bilayer][xyview]->GetNbinsX();
	hAveLayerOcc->Fill(tower, bilayer*g_nView+xyview, ave_rate);

      }
    }  
  }
  hMaxStripOcc->Write("",TObject::kOverwrite);
  hMaxLayerOcc->Write("",TObject::kOverwrite);
  hAveStripOcc->Write("",TObject::kOverwrite);
  hAveLayerOcc->Write("",TObject::kOverwrite);


  /*
  //Overwite histgram
  for (tower=0;tower<g_nTower; tower++){
    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      hTkrExposure [tower][bilayer]->Write("",TObject::kOverwrite);
      for (xyview=0; xyview<g_nView; xyview++) { 
	hTkrStripOcc [tower][bilayer][xyview]->Write("",TObject::kOverwrite);
	hTkrLayerOcc [tower][bilayer][xyview]->Write("",TObject::kOverwrite);
	hTkrHitMap   [tower][bilayer][xyview]->Write("",TObject::kOverwrite);
	hTkrNoiseMul [tower][bilayer][xyview]->Write("",TObject::kOverwrite);
	hTkrNoiseTot0[tower][bilayer][xyview]->Write("",TObject::kOverwrite);
	hTkrNoiseTot1[tower][bilayer][xyview]->Write("",TObject::kOverwrite);
      }
    }
  }
  */


  dirTkrNoise->Write("",TObject::kOverwrite);
}
