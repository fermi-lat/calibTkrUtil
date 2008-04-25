#include "calibTkrUtil/TkrNoiseOcc.h"

int getUnpFromBilayer( int xyview, int bilayer ){
  int unp = -1;
  if(xyview==0){
    unp = 4 * ( bilayer/2 ) + 2;
    if(bilayer%2==0) unp--;
    //std::cout << "X" << bilayer << ": " << unp << std::endl;
  }
  else{
    unp = 4 * ( (bilayer+1)/2 );
    if(bilayer%2!=0) unp--;
    //std::cout << "Y" << bilayer << ": " << unp << std::endl;
  }	
  return unp;
}


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

void 
TkrNoiseOcc::initAnalysis(int nEvent, int evt_interval){

  int tower, bilayer, xyview;
  int ix;

  m_nEvent          = nEvent;	   
  m_evt_interval    = evt_interval;   

  //Default values
  m_coincidence_cut = 1;  //coincidence_cut;
  m_multi_ld        = -1; //multi_ld;       
  m_multi_hd        = 65; //multi_hd;
  m_periodic_trig   = 0;  //periodic_trig;
  m_trig_cut = 0;

  m_crit_strip_rate = 5.0e-5;
  m_crit_layer_rate = 5.0e-2;

  m_nx = (int)(m_nEvent/m_evt_interval)+1;
  m_event_counter   = 0;
  
  vEvtTime = new double[m_nx];

  for(tower=0;tower<g_nTower; tower++){
    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      //vTkrExposure[tower][bilayer] = new float[m_nx];
      //for(ix=0; ix<m_nx; ix++) vTkrExposure[tower][bilayer][ix]=0.0;
      for(xyview=0; xyview<g_nView; xyview++) {

	vTkrExposure[tower][bilayer][xyview]   = new float[m_nx];
	vTkrStripOcc[tower][bilayer][xyview]   = new float[m_nx];
	vTkrLayerOcc[tower][bilayer][xyview]   = new float[m_nx];
	vTkrHitMap[tower][bilayer][xyview]     = new float[g_nStripsPerLayer];
	vTkrWHitMap[tower][bilayer][xyview]     = new float[g_nStripsPerLayer];
	vTkrNoiseMul[tower][bilayer][xyview]   = new float[MUL_MAX];
	vTkrNoiseTot0[tower][bilayer][xyview]  = new float[TOT_MAX];
	vTkrNoiseTot1[tower][bilayer][xyview]  = new float[TOT_MAX];

	for(ix=0; ix<m_nx; ix++)   vTkrExposure[tower][bilayer][xyview][ix]  =0.0;
	for(ix=0; ix<m_nx; ix++)   vTkrStripOcc[tower][bilayer][xyview][ix]  =0.0;
	for(ix=0; ix<m_nx; ix++)   vTkrLayerOcc[tower][bilayer][xyview][ix]  =0.0;
	for(ix=0; ix<g_nStripsPerLayer; ix++) vTkrHitMap[tower][bilayer][xyview][ix] =0.0;
	for(ix=0; ix<g_nStripsPerLayer; ix++) vTkrWHitMap[tower][bilayer][xyview][ix] =0.0;
	for(ix=0; ix<150; ix++)  vTkrNoiseMul[tower][bilayer][xyview][ix]  =0.0;
	for(ix=0; ix<300; ix++)  vTkrNoiseTot0[tower][bilayer][xyview][ix] =0.0;
	for(ix=0; ix<300; ix++)  vTkrNoiseTot1[tower][bilayer][xyview][ix] =0.0;
      }
    }
  }
  
}

void 
TkrNoiseOcc::setCoincidenceCut(int coincidence_cut){
  m_coincidence_cut = coincidence_cut;
}

void 
TkrNoiseOcc::setMultiRange(int multi_ld, int multi_hd){
  m_multi_ld = multi_ld;       
  m_multi_hd = multi_hd;
}

void 
TkrNoiseOcc::setPeriodicTrigCut(int periodic_trig){
  m_periodic_trig = periodic_trig;
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
  m_evtTime = m_digiEvt->getTimeStamp(); 

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

      /*
      /// check track hits on the adjacent layers
      if (m_coincidence_cut>0) {
	//if ( (numHitLayer[tower][bilayer][0]>0) && (numHitLayer[tower][bilayer][1]>0) ) {
	//  continue;
	//}
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
      */
      
      for(xyview=0; xyview<g_nView; xyview++){ 
	
	
	/// check track hits on the adjacent layers
	if (m_coincidence_cut>0) {
	  if (numHitLayer[tower][bilayer][(xyview+1)%2]>0) {
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
	} // coincidence_cut
	// Fill Exposure
	vTkrExposure[tower][bilayer][xyview][(int)(m_event_counter/m_evt_interval)] +=1.0;
	
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
	  vTkrWHitMap[tower][bilayer][xyview][buf_stripId[tower][bilayer][xyview][ihit]] +=1.0/numHitLayer[tower][bilayer][xyview];
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

  if ((m_event_counter%m_evt_interval)==0) {
     vEvtTime[(int)(m_event_counter/m_evt_interval)] = m_evtTime;
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
      for(xyview=0; xyview<g_nView; xyview++) {
	delete[] vTkrExposure[tower][bilayer][xyview];
	delete[] vTkrStripOcc[tower][bilayer][xyview]; 
	delete[] vTkrLayerOcc[tower][bilayer][xyview]; 
	delete[] vTkrHitMap[tower][bilayer][xyview];   
	delete[] vTkrWHitMap[tower][bilayer][xyview];   
	delete[] vTkrNoiseMul[tower][bilayer][xyview]; 
	delete[] vTkrNoiseTot0[tower][bilayer][xyview];
	delete[] vTkrNoiseTot1[tower][bilayer][xyview];
	
	vTkrExposure[tower][bilayer][xyview]   = NULL;
	vTkrStripOcc[tower][bilayer][xyview]   = NULL;
	vTkrLayerOcc[tower][bilayer][xyview]   = NULL;
	vTkrHitMap[tower][bilayer][xyview]     = NULL;
	vTkrWHitMap[tower][bilayer][xyview]     = NULL;
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

  std::cout << "save noise histograms." << std::endl;

  dirTkrNoise->cd();
  xmi = 0.5;
  xma = (float)(m_nx*m_evt_interval)+0.5;

  char hname[80], htitle[80], dirname[80];

  ///  histgram definitions
  //TDirectory *dirTower;
  TDirectory *dirExposure, *dirStripOcc, *dirLayerOcc, *dirHitmap, *dirMulti, *dirTOT;
  TDirectory *dirExposureRt, *dirStripOccRt, *dirLayerOccRt, *dirHitmapRt, *dirMultiRt, *dirTOTRt;

  TDirectory *dirTowerAveOcc;

  dirExposureRt  = dirTkrNoise->mkdir("Exposure");
  dirStripOccRt  = dirTkrNoise->mkdir("StripOcc");
  dirLayerOccRt  = dirTkrNoise->mkdir("LayerOcc");
  dirHitmapRt    = dirTkrNoise->mkdir("Hitmap");
  dirMultiRt     = dirTkrNoise->mkdir("Multi");
  dirTOTRt       = dirTkrNoise->mkdir("TOT");

  dirTowerAveOcc = dirTkrNoise->mkdir("TowerAveOcc");

  TH1D *hEvtTime;
  TH1F *hTkrTowerAveOcc[g_nTower];
  TH1F *hTkrTowerSumExp[g_nTower];


  TH1F *hTkrExposure[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrStripOcc[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrLayerOcc[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrHitMap[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrWHitMap[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrNoiseMul[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrNoiseTot0[g_nTower][g_nTkrLayer][g_nView];
  TH1F *hTkrNoiseTot1[g_nTower][g_nTkrLayer][g_nView];

  sprintf(hname, "hEvtTime");
  sprintf(htitle, "Event Time (UT from 2001/1/1)");
  hEvtTime  = new TH1D(hname, htitle, m_nx, xmi, xma);
  for (ix=0; ix<m_nx; ix++) hEvtTime->SetBinContent(ix+1, vEvtTime[ix]);
  hEvtTime->Write("",TObject::kOverwrite);
  Double_t exposure[g_nTower][g_nTkrLayer*g_nView];

  for(tower=0;tower<g_nTower; tower++){

    sprintf(dirname, "Tower%d", tower);
    //dirTower = dirTkrNoise->mkdir(dirname);
    dirExposure  = dirExposureRt->mkdir(dirname);
    dirStripOcc  = dirStripOccRt->mkdir(dirname);
    dirLayerOcc  = dirLayerOccRt->mkdir(dirname);
    dirHitmap    = dirHitmapRt->mkdir(dirname);
    dirMulti     = dirMultiRt->mkdir(dirname);
    dirTOT       = dirTOTRt->mkdir(dirname);

    //Tower Average Occupancy
    dirTowerAveOcc->cd();
    sprintf(hname, "hTkrTowerAveOccTwr%d", tower);
    sprintf(htitle, "Tower%d-TKR tower-average noise occupancy", tower);
    hTkrTowerAveOcc[tower] = new TH1F(hname, htitle, m_nx, xmi, xma);
     
     //Sum of esposures of entire layers in Tower
    dirTowerAveOcc->cd();
    sprintf(hname, "hTkrTowerSumExpTwr%d", tower);
    sprintf(htitle, "Tower%d-TKR exposure of entire layer in Tower", tower);
    hTkrTowerSumExp[tower]  = new TH1F(hname, htitle, m_nx, xmi, xma);

    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {

      for(xyview=0; xyview<g_nView; xyview++) {

	// convert layer/view to layer name
	char lname[] = "X15";
	if( xyview == 0 ) sprintf(lname,"X%d", bilayer);
	else sprintf(lname,"Y%d", bilayer);

	// Exposure
	dirExposure->cd();
	//sprintf(hname, "hTkrExposTwr%dLayer%d",  tower, bilayer*g_nView+xyview);
	//sprintf(htitle, " TKR  Tower%d Layer%d Exposure", tower, bilayer*g_nView+xyview);
	sprintf(hname, "hTkrExposT%d%s",  tower, lname);
	sprintf(htitle, " TKR  Tower%d Layer %s Exposure", tower, lname);
	hTkrExposure[tower][bilayer][xyview] = new TH1F(hname, htitle, m_nx, xmi, xma);
	for (ix=0; ix<m_nx; ix++) hTkrExposure[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrExposure[tower][bilayer][xyview][ix]);

	// convert layer/view to uni-plane
	int unp = getUnpFromBilayer( xyview, bilayer );
	exposure[tower][unp] = hTkrExposure[tower][bilayer][xyview]->Integral();
      
	// StripOcc
	dirStripOcc->cd();
	sprintf(hname, "hTkrStripOccT%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Averaged Hit-Strip Occupancy per Layer %s", tower, lname);
	hTkrStripOcc[tower][bilayer][xyview]  = new TH1F(hname, htitle, m_nx, xmi, xma);
	for (ix=0; ix<m_nx; ix++) hTkrStripOcc[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrStripOcc[tower][bilayer][xyview][ix]);
	hTkrStripOcc[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	// LayerOcc
	dirLayerOcc->cd();
	sprintf(hname, "hTkrLayerOccT%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer Occupancy per Layer %s", tower, lname);
	hTkrLayerOcc[tower][bilayer][xyview]  = new TH1F(hname, htitle, m_nx, xmi, xma);
	for (ix=0; ix<m_nx; ix++) hTkrLayerOcc[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrLayerOcc[tower][bilayer][xyview][ix]);
	hTkrLayerOcc[tower][bilayer][xyview]->Write("",TObject::kOverwrite);
	
	// For tower average occupancy
	hTkrTowerAveOcc[tower]->Add(hTkrStripOcc[tower][bilayer][xyview]);
	hTkrTowerSumExp[tower]->Add(hTkrExposure[tower][bilayer][xyview]);


	///// Tune histogram scale
	// scale of Strip Occupancy
	hTkrStripOcc[tower][bilayer][xyview]->Divide(hTkrExposure[tower][bilayer][xyview]);
	hTkrStripOcc[tower][bilayer][xyview]->Scale(1.0/g_nStripsPerLayer);
	hTkrStripOcc[tower][bilayer][xyview]->Write(0,TObject::kOverwrite);

	// scale of Layer Occupancy
	hTkrLayerOcc[tower][bilayer][xyview]->Divide(hTkrExposure[tower][bilayer][xyview]);
	hTkrLayerOcc[tower][bilayer][xyview]->Write(0,TObject::kOverwrite);

	// check maximum strip rate
	max_rate = hTkrStripOcc[tower][bilayer][xyview]->GetMaximum();
	if (max_rate>m_crit_strip_rate) {
	  std::cout << " Tower" << tower;
	  std::cout << " Layer" << lname;
	  std::cout << std::scientific << " Max Strip Occupancy=" << max_rate << std::endl;
	}
	// check maximum layer rate
	max_rate = hTkrLayerOcc[tower][bilayer][xyview]->GetMaximum();
	if (max_rate>m_crit_layer_rate) {
	  std::cout << " Tower" << tower;
	  std::cout << " Layer" << lname;
	  std::cout << std::scientific << " Max Layer ccupancy=" << max_rate << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////
	// Hitmap
	dirHitmap->cd();
	sprintf(hname, "hTkrHitMapT%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s Hit Map", tower, lname);
	hTkrHitMap[tower][bilayer][xyview]  = new TH1F(hname, htitle, g_nStripsPerLayer, -0.5, g_nStripsPerLayer-0.5);
	for (ix=0; ix<g_nStripsPerLayer; ix++) hTkrHitMap[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrHitMap[tower][bilayer][xyview][ix]);
	hTkrHitMap[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	// WeightedHitmap
	dirHitmap->cd();
	sprintf(hname, "hTkrWHitMapT%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s Weighted Hit Map", tower, lname);
	hTkrWHitMap[tower][bilayer][xyview]  = new TH1F(hname, htitle, g_nStripsPerLayer, -0.5, g_nStripsPerLayer-0.5);
	for (ix=0; ix<g_nStripsPerLayer; ix++) hTkrWHitMap[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrWHitMap[tower][bilayer][xyview][ix]);
	hTkrWHitMap[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	// Multiplicity
	dirMulti->cd();
	sprintf(hname, "hTkrNoiseMulT%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s Noise Multiplicity", tower, lname);
	hTkrNoiseMul[tower][bilayer][xyview]  = new TH1F(hname, htitle, 150, -0.5, 149.5);
	for (ix=0; ix<150; ix++) hTkrNoiseMul[tower][bilayer][xyview]-> SetBinContent(ix+1, vTkrNoiseMul[tower][bilayer][xyview][ix]);
	hTkrNoiseMul[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	// TOT
	dirTOT->cd();
	sprintf(hname, "hTkrNoiseTot0T%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s ToT0", tower, lname);
	hTkrNoiseTot0[tower][bilayer][xyview]  = new TH1F(hname, htitle, 300, -0.5, 299.5);
	for (ix=0; ix<300; ix++) hTkrNoiseTot0[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrNoiseTot0[tower][bilayer][xyview][ix]); 
	hTkrNoiseTot0[tower][bilayer][xyview]->Write("",TObject::kOverwrite);

	sprintf(hname, "hTkrNoiseTot1T%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s ToT1", tower, lname);
	hTkrNoiseTot1[tower][bilayer][xyview]  = new TH1F(hname, htitle, 300, -0.5, 299.5);
	for (ix=0; ix<300; ix++) hTkrNoiseTot1[tower][bilayer][xyview]->SetBinContent(ix+1, vTkrNoiseTot1[tower][bilayer][xyview][ix]);
	hTkrNoiseTot1[tower][bilayer][xyview]->Write("",TObject::kOverwrite);


      } //xyview
    } // bilayer

    // Tower Average Occupancy scale
    hTkrTowerAveOcc[tower]->Divide(hTkrTowerSumExp[tower]);
    hTkrTowerAveOcc[tower]->Scale(1.0/g_nStripsPerLayer);
    hTkrTowerAveOcc[tower]->Write("",TObject::kOverwrite);
    hTkrTowerSumExp[tower]->Write("",TObject::kOverwrite);

  }
  dirTkrNoise->cd();
  
  //
  // save timestamp TTree
  //
  char leaf[] = "exposures[16][36]/D";
  sprintf( leaf, "exposures[%d][%d]/D", g_nTower, g_nTkrLayer*g_nView);
  TTree *tree = new TTree("tkrNoiseTree","tkrNoiseTree");
  tree->Branch("exposures",&exposure,leaf);
  tree->Fill();
  tree->Write();

  /// Summary plot
  TH2F *hMaxStripOcc = new TH2F("hMaxStripOcc", "Max Noise Occupancy averaged over layer", g_nTower, -0.5, g_nTower-0.5, 36, -0.5, 35.5);
  TH2F *hMaxLayerOcc = new TH2F("hMaxLayerOcc", "Max Layer Occupancy averaged over layer", g_nTower, -0.5, g_nTower-0.5, 36, -0.5, 35.5);
  TH2F *hAveStripOcc = new TH2F("hAveStripOcc", "Average Noise Occupancy averaged over layer", g_nTower, -0.5, g_nTower-0.5, 36, -0.5, 35.5);
  TH2F *hAveLayerOcc = new TH2F("hAveLayerOcc", "Average Layer Occupancy averaged over layer", g_nTower, -0.5, g_nTower-0.5, 36, -0.5, 35.5);

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
  
  dirTkrNoise->Write("",TObject::kOverwrite);
}
