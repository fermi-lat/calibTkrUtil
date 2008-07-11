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
TkrNoiseOcc::initAnalysis(UInt_t duration, UInt_t nx){

  int tower, bilayer, xyview;
  int ix;

  //m_nEvent          = nEvent;	   
  //m_evt_interval    = evt_interval;   

  //Default values
  m_coincidence_cut = 1;  //coincidence_cut;
  m_multi_ld        = -1; //multi_ld;       
  m_multi_hd        = 65; //multi_hd;
  m_periodic_trig   = 0;  //periodic_trig;
  m_trig_cut = 0;

  m_crit_strip_rate = 5.0e-5;
  m_crit_layer_rate = 5.0e-2;

  //m_nx = (int)(m_nEvent/m_evt_interval)+1;
  //m_event_counter   = 0;
  
  m_duration = duration; // minutes to monitor per run
  m_nx = nx; // number of bins per run
  m_iP = 0;

  for(tower=0;tower<g_nTower; tower++){
    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      for(xyview=0; xyview<g_nView; xyview++) {
	vTkrHitMap[tower][bilayer][xyview]     = new float[g_nStripsPerLayer];
	vTkrWHitMap[tower][bilayer][xyview]     = new float[g_nStripsPerLayer];
	vTkrNoiseMul[tower][bilayer][xyview]   = new float[MUL_MAX];
	vTkrNoiseTot0[tower][bilayer][xyview]  = new float[TOT_MAX];
	vTkrNoiseTot1[tower][bilayer][xyview]  = new float[TOT_MAX];
	for(ix=0; ix<g_nStripsPerLayer; ix++) 
	  vTkrHitMap[tower][bilayer][xyview][ix] = 0.0;
	for(ix=0; ix<g_nStripsPerLayer; ix++) 
	  vTkrWHitMap[tower][bilayer][xyview][ix] = 0.0;
	for(ix=0; ix<MUL_MAX; ix++)
	  vTkrNoiseMul[tower][bilayer][xyview][ix]  = 0.0;
	for(ix=0; ix<TOT_MAX; ix++){
	  vTkrNoiseTot0[tower][bilayer][xyview][ix] = 0.0;
	  vTkrNoiseTot1[tower][bilayer][xyview][ix] = 0.0;
	}
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
      //m_event_counter +=1;
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
  RunInfo runInfo = m_digiEvt->getMetaEvent().run();

  // check if new run
  //UInt_t runid = runInfo.id();
  UInt_t runid = m_digiEvt->getRunId(); 
  if( vParamTimeDep.size() == 0 || 
      runid != vParamTimeDep[m_iP].id() ){
    // check if there is the same runid
    bool newrun = true;
    for(UInt_t ip=0; ip<vParamTimeDep.size(); ip++)
      if( runid == vParamTimeDep[ip].id() ){
	newrun = false;
	m_iP = ip;
	std::cout << "Old run: " << runid << ", event id: " << ievent
		  << ", start time: " << runInfo.startTime() << std::endl;
      }
    // new run id
    if( newrun ){
      m_iP = vParamTimeDep.size();
      vParamTimeDep.push_back( paramTimeDep( runid, 
					     runInfo.startTime(), 
					     m_duration, m_nx ) );
      std::cout << "New run: " << runid << ", event id: " << ievent
		<< ", start time: " << runInfo.startTime() << std::endl;
    }
  }


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
	
	// Cut by Multiplicity
	if( (numHitLayer[tower][bilayer][xyview]<m_multi_ld)
	    || (numHitLayer[tower][bilayer][xyview]>m_multi_hd) ){
	  continue;
	  vParamTimeDep[m_iP].fill( tower, bilayer, xyview, m_evtTime, 1, 0 );
	}

	// Fill time dependent variables
	vParamTimeDep[m_iP].fill( tower, bilayer, xyview, m_evtTime,
				  1, numHitLayer[tower][bilayer][xyview] );
	
	// Fill Noise Multiplicity
	vTkrNoiseMul[tower][bilayer][xyview][numHitLayer[tower][bilayer][xyview]] +=1.0;
	// Fill HitMap
	for(ihit=0; ihit<numHitLayer[tower][bilayer][xyview]; ihit++) {
	  vTkrHitMap[tower][bilayer][xyview][buf_stripId[tower][bilayer][xyview][ihit]] +=1.0;
	  vTkrWHitMap[tower][bilayer][xyview][buf_stripId[tower][bilayer][xyview][ihit]] +=1.0/numHitLayer[tower][bilayer][xyview];
	}
	//Fill Tot
	val = tot0[tower][bilayer][xyview];
	if( val>TOT_MAX-2 ) val = TOT_MAX-2;
	if( val<0 ) val = TOT_MAX-1;
	vTkrNoiseTot0[tower][bilayer][xyview][val] +=1.0;
	
	if ((err_cout==1)&&(val>=TOT_MAX-2)) {
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

  //if ((m_event_counter%m_evt_interval)==0) {
  //   vEvtTime[(int)(m_event_counter/m_evt_interval)] = m_evtTime;
  //}
  //m_event_counter +=1;

}  


void 
TkrNoiseOcc::clearAnalysis() {

  int tower, bilayer, xyview;
  //int m_nx;

  //m_event_counter = 0;
  /// delete histgram array
  for(tower=0;tower<g_nTower; tower++){
    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      for(xyview=0; xyview<g_nView; xyview++) {
	delete[] vTkrHitMap[tower][bilayer][xyview];   
	delete[] vTkrWHitMap[tower][bilayer][xyview];   
	delete[] vTkrNoiseMul[tower][bilayer][xyview]; 
	delete[] vTkrNoiseTot0[tower][bilayer][xyview];
	delete[] vTkrNoiseTot1[tower][bilayer][xyview];
	
	vTkrHitMap[tower][bilayer][xyview]     = NULL;
	vTkrWHitMap[tower][bilayer][xyview]     = NULL;
	vTkrNoiseMul[tower][bilayer][xyview]   = NULL;
	vTkrNoiseTot0[tower][bilayer][xyview]  = NULL;
	vTkrNoiseTot1[tower][bilayer][xyview]  = NULL;
      }
    }
  }
  for( UInt_t i=0; i<vParamTimeDep.size(); i++ )
    vParamTimeDep[i].clear();
  vParamTimeDep.clear();
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
  UInt_t ix;

  //float max_rate, ave_rate;
  float xmi, xma;

  if(dirTkrNoise==NULL) {
    std::cout << "No output directory is opened" << std::endl;
    return;
  }

  std::cout << "save noise histograms." << std::endl;

  dirTkrNoise->cd();
  //xmi = 0.5;
  //xma = (float)(m_nx*m_evt_interval)+0.5;
  xmi = 0.0;
  xma = (float)m_duration;

  char hname[80], htitle[120], dirname[80];

  ///  histgram definitions
  //TDirectory *dirTower;
  TDirectory *dirTimeDep, *dirHitmap, *dirMulti, *dirTOT;
  TDirectory *dirTimeDepRt, *dirHitmapRt, *dirMultiRt, *dirTOTRt;

  TDirectory *dirTowerOcc;

  dirHitmapRt    = dirTkrNoise->mkdir("Hitmap");
  dirMultiRt     = dirTkrNoise->mkdir("Multi");
  dirTOTRt       = dirTkrNoise->mkdir("TOT");

  TH1F *hTkrHitMap;
  TH1F *hTkrWHitMap;
  TH1F *hTkrNoiseMul;
  TH1F *hTkrNoiseTot0;
  TH1F *hTkrNoiseTot1;

  /*
  sprintf(hname, "hEvtTime");
  sprintf(htitle, "Event Time (UT from 2001/1/1)");
  hEvtTime  = new TH1D(hname, htitle, m_nx, xmi, xma);
  for (ix=0; ix<m_nx; ix++) hEvtTime->SetBinContent(ix+1, vEvtTime[ix]);
  hEvtTime->Write("",TObject::kOverwrite);
  */

  for(tower=0;tower<g_nTower; tower++){

    sprintf(dirname, "Tower%d", tower);
    //dirTower = dirTkrNoise->mkdir(dirname);
    dirHitmap    = dirHitmapRt->mkdir(dirname);
    dirMulti     = dirMultiRt->mkdir(dirname);
    dirTOT       = dirTOTRt->mkdir(dirname);

    for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      for(xyview=0; xyview<g_nView; xyview++) {

	// convert layer/view to layer name
	char lname[] = "X15";
	if( xyview == 0 ) sprintf(lname,"X%d", bilayer);
	else sprintf(lname,"Y%d", bilayer);
	// convert layer/view to uni-plane
	//int unp = getUnpFromBilayer( xyview, bilayer );

	// Hitmap
	dirHitmap->cd();
	sprintf(hname, "hTkrHitMapT%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s Hit Map", tower, lname);
	hTkrHitMap = new TH1F(hname, htitle, g_nStripsPerLayer, 
			      -0.5, g_nStripsPerLayer-0.5);
	// WeightedHitmap
	sprintf(hname, "hTkrWHitMapT%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s Weighted Hit Map", tower, lname);
	hTkrWHitMap = new TH1F(hname, htitle, g_nStripsPerLayer, -0.5, g_nStripsPerLayer-0.5);
	for (ix=0; ix<g_nStripsPerLayer; ix++){
	  hTkrHitMap->SetBinContent(ix+1, vTkrHitMap[tower][bilayer][xyview][ix]);
	  hTkrWHitMap->SetBinContent(ix+1, vTkrWHitMap[tower][bilayer][xyview][ix]);
	}
	hTkrHitMap->Write("",TObject::kOverwrite);
	hTkrWHitMap->Write("",TObject::kOverwrite);

	// Multiplicity
	dirMulti->cd();
	sprintf(hname, "hTkrNoiseMulT%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s Noise Multiplicity", tower, lname);
	hTkrNoiseMul = new TH1F(hname, htitle, MUL_MAX, -0.5, MUL_MAX-0.5);
	for (ix=0; ix<MUL_MAX; ix++) 
	  hTkrNoiseMul-> SetBinContent(ix+1, vTkrNoiseMul[tower][bilayer][xyview][ix]);
	hTkrNoiseMul->Write("",TObject::kOverwrite);

	// TOT
	dirTOT->cd();
	sprintf(hname, "hTkrNoiseTot0T%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s ToT0", tower, lname);
	hTkrNoiseTot0 = new TH1F(hname, htitle, TOT_MAX, -0.5, TOT_MAX-0.5);
	sprintf(hname, "hTkrNoiseTot1T%d%s", tower, lname);
	sprintf(htitle, "Tower%d-TKR Layer %s ToT1", tower, lname);
	hTkrNoiseTot1 = new TH1F(hname, htitle, TOT_MAX, -0.5, TOT_MAX-0.5);
	for (ix=0; ix<TOT_MAX; ix++){
	  hTkrNoiseTot0->SetBinContent(ix+1, vTkrNoiseTot0[tower][bilayer][xyview][ix]); 
	  hTkrNoiseTot1->SetBinContent(ix+1, vTkrNoiseTot1[tower][bilayer][xyview][ix]);
	}
	hTkrNoiseTot1->Write("",TObject::kOverwrite);
	hTkrNoiseTot0->Write("",TObject::kOverwrite);

      } //xyview
    } // bilayer
  } // tower

  //
  // Time dependent plots.
  //

  //TH1D *hEvtTime;
  TH1F *hTkrTowerOcc;
  TH1F *hTkrTowerSumExp;
  TH1F *hTkrExposure;
  TH1F *hTkrStripOcc;
  TH1F *hTkrLayerOcc;
  /// Summary plot
  UInt_t ybin = g_nView * g_nTkrLayer;
  float ymax = ybin - 0.5;
  /*
  dirTkrNoise->cd();
  TH2F *hMaxStripOcc = new TH2F("hMaxStripOcc", 
				"Max Noise Occupancy averaged over layer", 
				g_nTower, -0.5, g_nTower-0.5, 
				ybin, -0.5, ymax);
  TH2F *hMaxLayerOcc = new TH2F("hMaxLayerOcc", 
				"Max Layer Occupancy averaged over layer", 
				g_nTower, -0.5, g_nTower-0.5, 
				ybin, -0.5, ymax);
  */
  TH2F *hExposureMap = new TH2F("hExposureMap", "Exposure per layer",
				g_nTower, -0.5, g_nTower-0.5, 
				ybin, -0.5, ymax);
  TH2F *hStripOccMap = new TH2F("hStripOccMap", 
				"Average Noise Occupancy averaged over layer", 
				g_nTower, -0.5, g_nTower-0.5, 
				ybin, -0.5, ymax);
  TH2F *hLayerOccMap = new TH2F("hLayerOccMap", 
				"Average Layer Occupancy averaged over layer", 
				g_nTower, -0.5, g_nTower-0.5, 
				ybin, -0.5, ymax);
  
  Double_t exposure[g_nTower][g_nTkrLayer*g_nView];
  float *vExposure, *vStripOcc, *vLayerOcc;
  
  for(UInt_t ip=0; ip<vParamTimeDep.size(); ip++){
    sprintf(dirname, "Run%d", vParamTimeDep[ip].id());
    std::cout << "New directory: " << dirname << std::endl;
    dirTimeDepRt = dirTkrNoise->mkdir(dirname);
    dirTowerOcc = dirTimeDepRt->mkdir("TowerOcc");
    for(tower=0;tower<g_nTower; tower++){
      
      // runid in name and title of histograms
      char nid[] = "-R0123456789";
      char runid[] = " RunID:0123456789";
      sprintf( nid, "-R%d", vParamTimeDep[ip].id() );
      sprintf( runid, " RunID:%d", vParamTimeDep[ip].id() );

      //Tower Average Occupancy
      dirTowerOcc->cd();
      sprintf(hname, "hTkrTowerOccTwr%d%s", tower, nid);
      sprintf(htitle, "Tower%d-TKR tower noise occupancy%s", 
	      tower, runid);
      hTkrTowerOcc = new TH1F(hname, htitle, m_nx, xmi, xma);
     
      //Sum of esposures of entire layers in Tower
      sprintf(hname, "hTkrTowerSumExpTwr%d%s", tower, nid);
      sprintf(htitle, "Tower%d-TKR exposure of entire layer in Tower%s", 
	      tower, runid);
      hTkrTowerSumExp  = new TH1F(hname, htitle, m_nx, xmi, xma );
      sprintf(dirname, "Tower%d", tower);
      dirTimeDep = dirTimeDepRt->mkdir(dirname);
      dirTimeDep->cd();
      
      for (bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
	for(xyview=0; xyview<g_nView; xyview++) {
	  // convert layer/view to layer name
	  char lname[] = "X15";
	  if( xyview == 0 ) sprintf(lname,"X%d", bilayer);
	  else sprintf(lname,"Y%d", bilayer);
	  // convert layer/view to uni-plane
	  int unp = getUnpFromBilayer( xyview, bilayer );

	  // Exposure
	  sprintf(hname, "hTkrExposT%d%s%s",  tower, lname, nid);
	  sprintf(htitle, " TKR  Tower%d Layer %s Exposure%s", 
		  tower, lname, runid);
	  hTkrExposure = new TH1F(hname, htitle, m_nx, xmi, xma);
	  vExposure = vParamTimeDep[ip].getExposure( tower, bilayer, xyview );
      
	  // StripOcc
	  sprintf(hname, "hTkrStripOccT%d%s%s", tower, lname, nid);
	  sprintf(htitle, "Tower%d- Averaged Strip Occupancy per Layer %s%s", 
		  tower, lname, runid);
	  hTkrStripOcc  = new TH1F(hname, htitle, m_nx, xmi, xma);
	  vStripOcc = vParamTimeDep[ip].getStripOcc( tower, bilayer, xyview );
	  
	  // LayerOcc
	  sprintf(hname, "hTkrLayerOccT%d%s%s", tower, lname, nid);
	  sprintf(htitle, "Tower%d- Layer-OR Occupancy per Layer %s%s", 
		  tower, lname, runid);
	  hTkrLayerOcc  = new TH1F(hname, htitle, m_nx, xmi, xma);
	  vLayerOcc = vParamTimeDep[ip].getLayerOcc( tower, bilayer, xyview );

	  //
	  //
	  // fill and save histograms
	  for (ix=0; ix<m_nx; ix++){
	    hTkrExposure->SetBinContent(ix+1, vExposure[ix]);
	    hTkrStripOcc->SetBinContent(ix+1, vStripOcc[ix]);
	    hTkrLayerOcc->SetBinContent(ix+1, vLayerOcc[ix]);
	  }	
	  exposure[tower][unp] = hTkrExposure->Integral();
	  hTkrExposure->Write("",TObject::kOverwrite);
	  hTkrStripOcc->Write("",TObject::kOverwrite);
	  hTkrLayerOcc->Write("",TObject::kOverwrite);
	  
	  // For tower average occupancy
	  hTkrTowerOcc->Add(hTkrStripOcc);
	  hTkrTowerSumExp->Add(hTkrExposure);
	  	  
	  // check maximum strip rate
	  /*
	  max_rate = hTkrStripOcc->GetMaximum();
	  if (max_rate>m_crit_strip_rate) {
	    std::cout << " Tower" << tower;
	    std::cout << " Layer" << lname;
	    std::cout << std::scientific << " Max Strip Occupancy=" << max_rate << std::endl;
	  }
	  // check maximum layer rate
	  max_rate = hTkrLayerOcc->GetMaximum();
	  if (max_rate>m_crit_layer_rate) {
	    std::cout << " Tower" << tower;
	    std::cout << " Layer" << lname;
	    std::cout << std::scientific << " Max Layer ccupancy=" << max_rate << std::endl;
	  }
	  */

	  // fill summary maps
	  /*
	  max_rate = hTkrStripOcc->GetMaximum();
	  hMaxStripOcc->Fill(tower, unp, max_rate);
	  max_rate = hTkrLayerOcc->GetMaximum();
	  hMaxLayerOcc->Fill(tower, unp, max_rate);
	  */
	  hExposureMap->Fill(tower, unp, hTkrExposure->GetSum() );
	  hStripOccMap->Fill(tower, unp, hTkrStripOcc->GetSum() );
	  hLayerOccMap->Fill(tower, unp, hTkrLayerOcc->GetSum() );
	} //xyview
      } // bilayer

      // Tower Average Occupancy scale
      dirTowerOcc->cd();
      //hTkrTowerAveOcc->Divide(hTkrTowerSumExp);
      //hTkrTowerAveOcc->Scale(1.0/g_nStripsPerLayer);
      hTkrTowerOcc->Write("",TObject::kOverwrite);
      hTkrTowerSumExp->Write("",TObject::kOverwrite);

    } // tower
  } // runid

  dirTkrNoise->cd();
  //hMaxStripOcc->Write("",TObject::kOverwrite);
  //hMaxLayerOcc->Write("",TObject::kOverwrite);
  hExposureMap->Write("",TObject::kOverwrite);
  hStripOccMap->Write("",TObject::kOverwrite);
  hLayerOccMap->Write("",TObject::kOverwrite);
  
  //
  // save TTree
  //
  TTree *tree = new TTree("tkrNoiseTree","tkrNoiseTree");

  UInt_t runid, startTime;
  tree->Branch( "runid", &runid, "runid/i");
  tree->Branch( "startTime", &startTime, "startTime/i");
  for( ix=0; ix<vParamTimeDep.size(); ix++){
    runid = vParamTimeDep[ix].id();
    startTime = vParamTimeDep[ix].startTime();
    tree->Fill();
  }
  tree->Write();
  
  dirTkrNoise->Write("",TObject::kOverwrite);
}



paramTimeDep::paramTimeDep( UInt_t id, UInt_t startT, UInt_t duration, 
			    UInt_t nx ){
  m_startTime = startT;
  m_duration = duration;
  m_nx = nx;
  m_binSize = 60.0 * m_duration / m_nx; // convert from second to minute
  m_id = id;

  for( UInt_t tower=0;tower<g_nTower; tower++){
    for( UInt_t bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      //m_Exposure[tower][bilayer] = new float[m_nx];
      //for( UInt_t ix=0; ix<m_nx; ix++) m_Exposure[tower][bilayer][ix]=0.0;
      for( UInt_t xyview=0; xyview<g_nView; xyview++) {
	m_Exposure[tower][bilayer][xyview] = new float[m_nx];
	m_StripOcc[tower][bilayer][xyview] = new float[m_nx];
	m_LayerOcc[tower][bilayer][xyview] = new float[m_nx];
	for( UInt_t ix=0; ix<m_nx; ix++){
	  m_Exposure[tower][bilayer][xyview][ix] = 0.0;
	  m_StripOcc[tower][bilayer][xyview][ix] = 0.0;
	  m_LayerOcc[tower][bilayer][xyview][ix] = 0.0;
	}
      }
    }
  }
}

void paramTimeDep::fill( int tower, int bilayer, int view, Double_t timeStamp,
			 float nExp, float nOcc ){
  
  if( nExp <= 0.0 ) return;
  Int_t ix = (int) (( timeStamp - m_startTime ) / m_binSize);
  if( ix >= m_nx ) ix = m_nx - 1;
  if( ix < 0 ) ix = 0;

  if( nOcc > 0 ){
    // Fill Strip Occupancy
    m_StripOcc[tower][bilayer][view][ix] += nOcc;
    // Fill layer-OR Occupancy
    m_LayerOcc[tower][bilayer][view][ix] += 1.0;
  }
  // Fill Exposure
  m_Exposure[tower][bilayer][view][ix] += 1.0;
}

void paramTimeDep::clear(){
  for( UInt_t tower=0;tower<g_nTower; tower++){
    for ( UInt_t bilayer=0; bilayer<g_nTkrLayer; bilayer++) {
      for( UInt_t xyview=0; xyview<g_nView; xyview++) {
	delete[] m_Exposure[tower][bilayer][xyview];
	delete[] m_StripOcc[tower][bilayer][xyview]; 
	delete[] m_LayerOcc[tower][bilayer][xyview]; 
	m_Exposure[tower][bilayer][xyview]   = NULL;
	m_StripOcc[tower][bilayer][xyview]   = NULL;
	m_LayerOcc[tower][bilayer][xyview]   = NULL;
      }
    }
  }
}
