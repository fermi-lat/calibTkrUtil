#include <TROOT.h>
#include "TkrNoiseRep.h"

TkrNoiseRep::TkrNoiseRep() 
{
  m_fail_strip_occ     = 1.0e-4;
  m_critical_strip_occ = 5.0e-5;
  m_critical_layer_occ = 0.08;
  m_critical_multi = 13; //20
  m_critical_multi_ratio = 1.0e-4;
}

TkrNoiseRep::~TkrNoiseRep() 
{}

void
TkrNoiseRep::openSvacFile(const char *filename){
  m_svacFile = new TFile(filename, "READ");
}

void 
TkrNoiseRep::findNoisyLayer(){

  int tower, bilayer, view;
  int layer;
  //double maxStripOcc, maxLayerOcc, multiRatio;
  char outputFilename[250]; 

  // initialize
  for(tower=0; tower<16; tower++){
    for(bilayer=0; bilayer<18; bilayer++){
      for(view=0; view<2; view++){
        m_failStripOcc[tower][bilayer][view]=0;
        m_failLayerOcc[tower][bilayer][view]=0;
        m_failAveLayerOcc[tower][bilayer][view]=0;
        m_failMultiRatio[tower][bilayer][view]=0;
      }
    }
  }

  sprintf(outputFilename, "%s/%s_BadLayer.txt", m_reportDirName, m_prefix);
  
  m_ancRootFile->cd();
  TH2F *hStrip = (TH2F*) gROOT->FindObject("hMaxStripOcc");
  TH2F *hLayer = (TH2F*) gROOT->FindObject("hMaxLayerOcc");
  TH2F *hMulti = (TH2F*) gROOT->FindObject("hLargeMultiRatio");
  TH2F *hAveLayer = (TH2F*) gROOT->FindObject("hAveLayerOcc");

  for(tower=0; tower<16; tower++){
    for(bilayer=0; bilayer<18; bilayer++){
      for(view=0; view<2; view++){
        layer = bilayer*2+view;
        
        if ( m_critical_strip_occ < (hStrip->GetBinContent(tower+1, layer+1)) ) {
          m_failStripOcc[tower][bilayer][view]=1;
        }
        if ( m_critical_layer_occ < (hLayer->GetBinContent(tower+1, layer+1)) ) {
          m_failLayerOcc[tower][bilayer][view]=1;
        }
        if ( m_critical_layer_occ < (hAveLayer->GetBinContent(tower+1, layer+1)) ) {
          m_failAveLayerOcc[tower][bilayer][view]=1;
        }
        if ( m_critical_multi_ratio < (hMulti->GetBinContent(tower+1, layer+1)) ) {
          m_failMultiRatio[tower][bilayer][view]=1;
        }
        
      }
    }
  }
  
}

void 
TkrNoiseRep::writeSummaryXML(){

  int  tower, bilayer, view;
  int  layer;
  char xmlFilename[250];
  char statusStr[10];

  if (m_test_status==1) strcpy(statusStr,"Pass");
  if (m_test_status==2) strcpy(statusStr,"Pass*");
  else strcpy(statusStr,"Fail");

  m_ancRootFile->cd();
  TH2F *hExpos = (TH2F*) gROOT->FindObject("hExposure");
  TH2F *hStrip = (TH2F*) gROOT->FindObject("hMaxStripOcc");
  TH2F *hLayer = (TH2F*) gROOT->FindObject("hMaxLayerOcc");
  TH2F *hMulti = (TH2F*) gROOT->FindObject("hLargeMultiRatio");

  sprintf(xmlFilename, "%s/%s_summary.xml", m_reportDirName, m_prefix);
  FILE *outfile = fopen(xmlFilename, "w");
  
  fprintf(outfile, "<?xml version=\"1.0\"?>\n");
  fprintf(outfile, "<tkrNoiseData>\n");
  fprintf(outfile, "<testRun status=\"%s\" runId=\"%s\" startTime=\"%s\" duration=\"%s\" minExposure=\"%.0f\"/>\n", statusStr, m_prefix, m_startTimeStr, m_durationStr, m_minExpos);
  // latAve
  fprintf(outfile, "<latAve>\n");
  fprintf(outfile, "<stripOcc mean=\"%.4e\" max=\"%.4e\"/>\n", m_latAveOccMean, m_latAveOccMax);
  fprintf(outfile, "</latAve>\n");
  // towerAve
  fprintf(outfile, "<towerAve>\n");
  for(tower=0; tower<16; tower++){
    fprintf(outfile, "<stripOcc tower=\"%d\" mean=\"%.4e\" max=\"%.4e\"/>\n", tower, m_aveTwrOccList[tower], m_maxTwrOccList[tower]);
  }
  fprintf(outfile, "</towerAve>\n");
  // bad layers
  fprintf(outfile, "<badLayers>\n");
  for(tower=0; tower<16; tower++){
    for(bilayer=0; bilayer<18; bilayer++){
      for(view=0; view<2; view++){
        layer = bilayer*2+view;
        if (m_failStripOcc[tower][bilayer][view]==1) {
          fprintf(outfile, "<stripOcc tower=\"%d\" biLayer=\"%d\" view=\"%d\" occupancy=\"%.6e\"/>\n", tower, bilayer, view, hStrip->GetBinContent(tower+1, layer+1));
        }
        if (m_failLayerOcc[tower][bilayer][view]==1) {
          fprintf(outfile, "<layerOcc tower=\"%d\" biLayer=\"%d\" view=\"%d\" occupancy=\"%.6e\"/>\n", tower, bilayer, view, hLayer->GetBinContent(tower+1, layer+1));
        }
        if (m_failMultiRatio[tower][bilayer][view]==1) {
          fprintf(outfile, "<multiRatio tower=\"%d\" biLayer=\"%d\" view=\"%d\" multiRatio=\"%.6e\"/>\n", tower, bilayer, view, hMulti->GetBinContent(tower+1, layer+1));
        }
      }
    }
  }

  fprintf(outfile, "</badLayers>\n");
  fprintf(outfile, "</tkrNoiseData>\n");

}

void 
TkrNoiseRep::makeAncRootFile(){

  char ancFilename[250];  
  sprintf(ancFilename, "%s/%s_anc.root", m_reportDirName, m_prefix); 
  m_ancRootFile = new TFile(ancFilename, "RECREATE");
  m_ancRootFile->cd();
  ////////////////////////////////
  ////// from makeNewHist.C //////
  ////////////////////////////////
  int tower, layer;
  int ix, nx;
  char dirname[80], hname[80], htitle[160];
  double exposure, num_strip, num_layerhit, ave_strip_occ, ave_layer_occ;
  //double ave_strip_occ_max[16];
  double tsum_exposure, tsum_num_strip, tsum_num_layerhit;
  double max_rate;
  double expos_limit = 100.0;
 
  
  sprintf(htitle, "Average Noise Strip Occupancy per Tower: RunID=%s", m_prefix);
  TH1F *hTwrAveStripOcc = new TH1F("hTwrAveStripOcc", htitle, 16, -0.5, 15.5);
  sprintf(htitle, "Layer Occupancy averaged over Tower: RunID=%s", m_prefix);
  TH1F *hTwrAveLayerOcc = new TH1F("hTwrAveLayerOcc", htitle, 16, -0.5, 15.5);
  
  sprintf(htitle, "Noise Strip Occupancy averaged per layer: RunID=%s", m_prefix);
  TH2F *hAveStripOcc = new TH2F("hAveStripOcc", htitle, 16, -0.5, 15.5, 36, -0.5, 35.5);
  sprintf(htitle, "Layer Occupancy: RunID=%s", m_prefix);
  TH2F *hAveLayerOcc = new TH2F("hAveLayerOcc", htitle, 16, -0.5, 15.5, 36, -0.5, 35.5);
  
  sprintf(htitle, "Transient Peak of Noise Strip Occupancy averaged per layer: RunID=%s", m_prefix);
  TH2F *hMaxStripOcc = new TH2F("hMaxStripOcc", htitle, 16, -0.5, 15.5, 36, -0.5, 35.5);
  sprintf(htitle, "Transient peak of Layer Occupancy: RunID=%s", m_prefix);
  TH2F *hMaxLayerOcc = new TH2F("hMaxLayerOcc", htitle, 16, -0.5, 15.5, 36, -0.5, 35.5);
  
  sprintf(htitle, "Number of Noise Events: RunID=%s", m_prefix);
  TH2F *hExposure = new TH2F("hExposure", htitle, 16, -0.5, 15.5, 36, -0.5, 35.5);


  for(tower=0; tower<16; tower++){

    tsum_exposure = 0.0;
    tsum_num_strip = 0.0;
    tsum_num_layerhit = 0.0;

    for(layer=0; layer<36; layer++){

      /// get histgram pointers
      
      // exposure
      sprintf(dirname, "TkrNoiseOcc/Exposure/Tower%d", tower);
      m_svacFile->cd(dirname);

      sprintf(hname, "hTkrExposTwr%dLayer%d", tower, layer);
      TH1F *hexpos = (TH1F*) gROOT->FindObject(hname);      
      if(hexpos==NULL) {
        sprintf(hname, "hTkrExposTwr%dbiLayer%d", tower, int(layer/2));
        hexpos = (TH1F*) gROOT->FindObject(hname);      
      }

      // Hitmap
      sprintf(dirname, "TkrNoiseOcc/Hitmap/Tower%d", tower);
      m_svacFile->cd(dirname);
      sprintf(hname, "hTkrHitMapTwr%dLayer%d", tower, layer);
      TH1F *hHitMap = (TH1F*) gROOT->FindObject(hname);

      // StripOcc
      sprintf(dirname, "TkrNoiseOcc/StripOcc/Tower%d", tower);
      m_svacFile->cd(dirname);
      sprintf(hname, "hTkrStripOccTwr%dLayer%d", tower, layer);
      TH1F *hStripOcc = (TH1F*) gROOT->FindObject(hname);

      // LayerOcc
      sprintf(dirname, "TkrNoiseOcc/LayerOcc/Tower%d", tower);
      m_svacFile->cd(dirname);
      sprintf(hname, "hTkrLayerOccTwr%dLayer%d", tower, layer);
      TH1F *hLayerOcc = (TH1F*) gROOT->FindObject(hname);
      

      /// hExposure ///
      exposure = hexpos->GetSum();
      hExposure->Fill(tower, layer, exposure);

      /// hAveStrip ///
      num_strip = hHitMap->GetSum();
      if (exposure>0.0) ave_strip_occ = num_strip/1536.0/exposure;
      else              ave_strip_occ = 0.0;
      hAveStripOcc->Fill(tower, layer, ave_strip_occ);
      
      tsum_exposure += exposure;
      tsum_num_strip += num_strip; 

      /// hAveLayer ///
      TH1F *hLayerHit = (TH1F*) hLayerOcc->Clone();
      hLayerHit->Multiply(hexpos);
      num_layerhit = hLayerHit->GetSum();
      if (exposure>0.0) ave_layer_occ = num_layerhit/exposure;
      else              ave_layer_occ = 0.0;
      hAveLayerOcc->Fill(tower, layer, ave_layer_occ);
      
      tsum_num_layerhit += num_layerhit; 

      nx = hexpos->GetNbinsX();
      if (hexpos->GetBinContent(nx)<expos_limit){
        hStripOcc->SetBinContent(nx, 0.0);
        hLayerOcc->SetBinContent(nx, 0.0);
      }
      /// hMaxStrip ///
      max_rate = hStripOcc->GetMaximum();
      hMaxStripOcc->Fill(tower, layer, max_rate);
      /// hMaxLayer ///
      max_rate = hLayerOcc->GetMaximum();
      hMaxLayerOcc->Fill(tower, layer, max_rate);
    }

    if (tsum_exposure>0.0) {
      hTwrAveStripOcc->Fill(tower, tsum_num_strip/1536.0/tsum_exposure);
      hTwrAveLayerOcc->Fill(tower, tsum_num_layerhit/tsum_exposure);
    } else {
      hTwrAveStripOcc->Fill(tower, 0.0); 
      hTwrAveLayerOcc->Fill(tower, 0.0);
    }
  }

  /// number of event  
  m_minExpos = hExposure->GetMinimum();

  ////////////////////////////////////////////////////////////
  ////// LAT average occipancy and Tower Peak Occupancy //////
  ////////////////////////////////////////////////////////////

  sprintf(dirname, "TkrNoiseOcc/TowerAveOcc");
  m_svacFile->cd(dirname);
  sprintf(hname, "hTkrTowerAveOccTwr0");
  TH1F *hOcc = (TH1F*)gROOT->FindObject(hname);
  
  m_ancRootFile->cd();
  TH1F *hLatOcc = hOcc->Clone();
  TH1F *hLatExp = hOcc->Clone();

  hLatOcc->SetName("hLatAveOcc");
  sprintf(htitle,"LAT Average Noise Occupancy :RunId=%s", m_prefix);
  hLatOcc->SetTitle(htitle);
  for(ix=1; ix<=hLatOcc->GetNbinsX(); ix++) {
    hLatOcc->SetBinContent(ix,0.0);
  }
  hLatExp->SetName("hLatExp");
  sprintf(htitle,"Sum of number of noise sampling events in all the layers of the 16-Tower LAT :RunId=%s", m_prefix);
  hLatExp->SetTitle(htitle);
  for(ix=1; ix<=hLatOcc->GetNbinsX(); ix++) {
    hLatExp->SetBinContent(ix,0.0);
  }
  
  sprintf(htitle, "Transient Maximum of Average Noise Strip Occupancy per Tower: RunID=%s", m_prefix);
  TH1F *hTwrAveStripOccMax = new TH1F("hTwrAveStripOccMax", htitle, 16, -0.5, 15.5);

  m_svacFile->cd(dirname);
  for(tower=0; tower<16; tower++) {
    sprintf(hname, "hTkrTowerAveOccTwr%d", tower);
    TH1F *hOcc = (TH1F*)gROOT->FindObject(hname);
    TH1F *hNumHits = (TH1F*) hOcc->Clone();
    sprintf(hname, "hTkrTowerSumExpTwr%d", tower);
    TH1F *hExpos = (TH1F*)gROOT->FindObject(hname);

    hTwrAveStripOccMax->Fill(tower, hOcc->GetMaximum());
    
    hNumHits->Multiply(hExpos);
    hLatOcc->Add(hNumHits);
    hLatExp->Add(hExpos);
  }
  
  hLatOcc->Divide(hLatExp);
  m_latAveOccMax  = hLatOcc->GetMaximum();
  m_latAveOccMean = hLatOcc->GetSum()/hLatOcc->GetNbinsX();

  
  //hNumHits->Delete();
  //hExpos->Delete();

  //////////////////////////////////////
  ////// Large Multiplicity Ratio //////
  //////////////////////////////////////
  int  numEvt, numMultiEvt;
  float fract;
  m_ancRootFile->cd();
  sprintf(hname, "hLargeMultiRatio");
  sprintf(htitle, "Fraction of noise events with strip multiplicity > %d", m_critical_multi);
  TH2F *hMultiRatio = new TH2F(hname, htitle, 16, -0.5, 15.5, 36, -0.5, 35.5);

  for(tower=0; tower<16; tower++) {
    sprintf(dirname, "TkrNoiseOcc/Multi/Tower%d", tower);
    m_svacFile->cd(dirname);
    for(layer=0; layer<36; layer++) {
      sprintf(hname, "hTkrNoiseMulTwr%dLayer%d", tower, layer); 
      TH1F *hMulti = (TH1F*) gROOT->FindObject(hname);
      numEvt = (int)hMulti->Integral(1,129);
      numMultiEvt = (int)hMulti->Integral(m_critical_multi,129);
      if (numEvt>0) fract = (float)numMultiEvt/(float)numEvt;
      else          fract = 0.0;
      hMultiRatio->Fill(tower, layer, fract);
    }
  }

  m_ancRootFile->Write();

  ///// check pass/fail /////
  m_test_status=1; // default : pass
  if ( m_fail_strip_occ < m_latAveOccMean ) m_test_status=0; //Fail
  if ( m_critical_strip_occ < m_latAveOccMean ) m_test_status=2;            //Warning
  if ( m_critical_layer_occ < hAveLayerOcc->GetMaximum() ) m_test_status=2; //Warning

}


void
TkrNoiseRep::setReportDirName(const char *reportDirName){
  m_reportDirName = reportDirName;
}

void
TkrNoiseRep::setPrefix(const char *prefix){
  //strcpy(m_prefix, prefix);
  m_prefix = prefix;
}

void
TkrNoiseRep::drawAncRootHist(const char *hname){

  char epsFilePath[250], epsFileName[250], gifFileName[250];

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); 

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetRightMargin(0.15);
  c1->SetLogz();

  m_ancRootFile->cd();
  TH2F *h2 = (TH2F*) gROOT->FindObject(hname);

  h2->GetXaxis()->SetTitle("Tower No.");
  h2->GetYaxis()->SetTitle("Layer No.");
  h2->Draw("COLZ");
  sprintf(epsFilePath, "%s/%s_%s.eps", m_reportDirName, m_prefix, hname);
  c1->SaveAs(epsFilePath);
  
  sprintf(epsFileName, "%s_%s.eps", m_prefix, hname);
  sprintf(gifFileName, "%s_%s.gif", m_prefix, hname);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_s.gif", m_prefix, hname);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);
}


void
TkrNoiseRep::drawLayerOccMax_sum(){

  const char *hname="hAveLayerOcc";
  char epsFilePath[250], epsFileName[250], gifFileName[250];
  
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); 

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetRightMargin(0.15);
  c1->SetLogz();

  m_ancRootFile->cd();
  TH2F *h2 = (TH2F*) gROOT->FindObject(hname);

  h2->SetContour(1);
  h2->SetMinimum(m_critical_layer_occ);
  h2->SetMaximum(1.0);
  h2->GetXaxis()->SetTitle("Tower No.");
  h2->GetYaxis()->SetTitle("Layer No.");
  h2->Draw("COLZ");
  sprintf(epsFilePath, "%s/%s_%s_sum.eps", m_reportDirName, m_prefix, hname);
  c1->SaveAs(epsFilePath);
  
  sprintf(epsFileName, "%s_%s_sum.eps", m_prefix, hname);
  sprintf(gifFileName, "%s_%s_sum.gif", m_prefix, hname);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_sum_s.gif", m_prefix, hname);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);
  
  h2->SetContour(0);
}


void
TkrNoiseRep::drawStripOccAve(){
  drawAncRootHist("hAveStripOcc");
}

void
TkrNoiseRep::drawStripOccMax(){
  drawAncRootHist("hMaxStripOcc");
}

void
TkrNoiseRep::drawLayerOccAve(){
  drawAncRootHist("hAveLayerOcc");
}

void
TkrNoiseRep::drawLayerOccMax(){
  drawAncRootHist("hMaxLayerOcc");
}

void
TkrNoiseRep::drawLargeMultiRatio(){
  drawAncRootHist("hLargeMultiRatio");
}

void
TkrNoiseRep::drawTowerAveStripOccHist(){

  char hname[80]; //dirname[80], htitle[80];
  char epsFilePath[250], epsFileName[250], gifFileName[250];
  char txtstr[20];

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  gStyle->SetOptStat(0);

  m_ancRootFile->cd();

  sprintf(hname, "hTwrAveStripOcc"); 
  TH1F *hTwrAveStripOcc = (TH1F*) gROOT->FindObject(hname);
  hTwrAveStripOcc->SetMinimum(5.0e-8);
  hTwrAveStripOcc->SetMaximum(2.0e-4);
  hTwrAveStripOcc->GetXaxis()->SetTitle("Tower No.");
  hTwrAveStripOcc->GetYaxis()->SetTitle("Strip Occupancy");
  hTwrAveStripOcc->Draw();
  
  sprintf(hname, "hTwrAveStripOccMax"); 
  TH1F *hTwrAveStripOccMax = (TH1F*) gROOT->FindObject(hname);
  hTwrAveStripOccMax->SetLineColor(2);
  hTwrAveStripOccMax->Draw("same");

  TLine *line1 = new TLine(-0.5, m_critical_strip_occ, 15.5, m_critical_strip_occ);
  line1->SetLineStyle(2);
  line1->SetLineColor(4);
  line1->Draw();

  sprintf(txtstr, "%.1e", m_critical_strip_occ);
  TText *text1 = new TText(15.5, m_critical_strip_occ, txtstr);
  text1->SetTextSize(0.04);
  text1->SetTextColor(4);
  text1->Draw();

  TLegend *leg = new TLegend(0.75,0.82, 0.9,0.9);
  leg->AddEntry(hTwrAveStripOcc,"Average","l");
  leg->AddEntry(hTwrAveStripOccMax,"Maximum","l");
  leg->SetBorderSize(1);
  leg->Draw();

  sprintf(epsFilePath, "%s/%s_TowerAveStripOcc.eps", m_reportDirName, m_prefix);
  c1->SaveAs(epsFilePath);

  sprintf(epsFileName, "%s_TowerAveStripOcc.eps", m_prefix);
  sprintf(gifFileName, "%s_TowerAveStripOcc.gif", m_prefix);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_TowerAveStripOcc_s.gif", m_prefix);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);

  hTwrAveStripOcc->SetMinimum();
  hTwrAveStripOcc->SetMaximum();
}

void
TkrNoiseRep::drawStripOcc_TowerAve(int tower){

  char hname[80]; //dirname[80], htitle[80];
  char epsFilePath[250], epsFileName[250], gifFileName[250];

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  gStyle->SetOptStat(0);

  m_svacFile->cd("TkrNoiseOcc/TowerAveOcc");
  sprintf(hname, "hTkrTowerAveOccTwr%d", tower); 
  TH1F *h_towerAveOcc = (TH1F*) gROOT->FindObject(hname);
   
  h_towerAveOcc->SetMinimum(5.0e-8);
  h_towerAveOcc->SetMaximum(2.0e-4);
  h_towerAveOcc->Draw();
  sprintf(epsFilePath, "%s/%s_TowerAveOcc_Tower%d.eps", m_reportDirName, m_prefix, tower);
  c1->SaveAs(epsFilePath);

  sprintf(epsFileName, "%s_TowerAveOcc_Tower%d.eps", m_prefix, tower);
  sprintf(gifFileName, "%s_TowerAveOcc_Tower%d.gif", m_prefix, tower);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_TowerAveOcc_Tower%d_s.gif", m_prefix, tower);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);

  h_towerAveOcc->SetMinimum();
  h_towerAveOcc->SetMaximum();
}

void
TkrNoiseRep::drawStripOccGr_TowerAve(int tower){
  
  char   dirname[80], hname[80]; //, t_title[80];
  int    ix, nx;
  double *vx; 
  double *vy;
  //double start_time;
  //time_t gmt_start_time;
  //int    offset_time;

  char   epsFilePath[250], epsFileName[250], gifFileName[250];

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  gStyle->SetOptStat(0);

  m_svacFile->cd("TkrNoiseOcc");
  TH1F *hEvtTime = (TH1F*) gROOT->FindObject("hEvtTime");
  nx = (hEvtTime->GetNbinsX());
  if ( (hEvtTime->GetBinContent(nx))==0.0 ) nx = nx-1;

  sprintf(dirname, "TkrNoiseOcc/TowerAveOcc");
  m_svacFile->cd(dirname);
  sprintf(hname, "hTkrTowerAveOccTwr%d", tower);
  TH1F *hStripOcc = (TH1F*) gROOT->FindObject(hname);

  vx = new double[nx];
  vy = new double[nx];

  for (ix=0; ix<nx; ix++) {
    vx[ix] = (double)hEvtTime->GetBinContent(ix+1)-m_start_time;
    vy[ix] = (double)hStripOcc->GetBinContent(ix+1);
  }
  
  TGraph *gr = new TGraph(nx, vx, vy);
  gr->SetMinimum(0.3e-7);
  gr->SetMaximum(0.3e-3);
  
  gr->SetMarkerStyle(7);
  gr->SetMarkerColor(2);
  gr->SetLineColor(2);

  gr->SetTitle(hStripOcc->GetTitle());
  gr->GetXaxis()->SetTitle(m_timeAxTitle);
  gr->Draw("AP");
  sprintf(epsFilePath, "%s/%s_TowerAveOccGr_Tower%d.eps", m_reportDirName, m_prefix, tower);
  c1->SaveAs(epsFilePath);

  sprintf(epsFileName, "%s_TowerAveOccGr_Tower%d.eps", m_prefix, tower);
  sprintf(gifFileName, "%s_TowerAveOccGr_Tower%d.gif", m_prefix, tower);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_TowerAveOccGr_Tower%d_s.gif", m_prefix, tower);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);

}

void
TkrNoiseRep::drawStripOccGr_LatAve(){
  
  //char   dirname[80], hname[80], htitle[120]; //, t_title[80];
  int    ix, nx; //, tower;
  double *vx; 
  double *vy;
  double xmax;
  char   epsFilePath[250], epsFileName[250], gifFileName[250];
  char txtstr[20];

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  gStyle->SetOptStat(0);

  m_svacFile->cd("TkrNoiseOcc");

  TH1F *hEvtTime = (TH1F*) gROOT->FindObject("hEvtTime");
  nx = (hEvtTime->GetNbinsX());
  if ( (hEvtTime->GetBinContent(nx))==0.0 ) nx = nx-1;

  m_ancRootFile->cd();
  TH1F *hLatOcc = (TH1F*) gROOT->FindObject("hLatAveOcc");

  vx = new double[nx];
  vy = new double[nx];

  for (ix=0; ix<nx; ix++) {
    vx[ix] = (double)hEvtTime->GetBinContent(ix+1)-m_start_time;
    vy[ix] = (double)hLatOcc->GetBinContent(ix+1);
  }
  
  TGraph *gr = new TGraph(nx, vx, vy);
  gr->SetMinimum(0.3e-7);
  gr->SetMaximum(0.3e-3);
  
  gr->SetMarkerStyle(7);
  gr->SetMarkerColor(2);
  gr->SetLineColor(2);

  gr->SetTitle(hLatOcc->GetTitle());
  gr->GetXaxis()->SetTitle(m_timeAxTitle);
  gr->Draw("AP");
  xmax = gr->GetXaxis()->GetXmax();

  TLine *line1 = new TLine(-0.5, m_critical_strip_occ, xmax, m_critical_strip_occ);
  line1->SetLineStyle(2);
  line1->SetLineColor(4);
  line1->Draw();

  sprintf(txtstr, "%.1e", m_critical_strip_occ);
  TText *text1 = new TText(xmax, m_critical_strip_occ, txtstr);
  text1->SetTextSize(0.04);
  text1->SetTextColor(4);
  text1->Draw();

  sprintf(epsFilePath, "%s/%s_LatAveOccGr.eps", m_reportDirName, m_prefix);
  c1->SaveAs(epsFilePath);

  sprintf(epsFileName, "%s_LatAveOccGr.eps", m_prefix);
  sprintf(gifFileName, "%s_LatAveOccGr.gif", m_prefix);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_LatAveOccGr_s.gif", m_prefix);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);

}

void
TkrNoiseRep::drawStripOcc_perLayer(int tower, int layer){

  char layerViewId[10], towerId[80];
  char dirname[80], hname[80], htitle[80];
  char epsFilePath[250], epsFileName[250], gifFileName[250];

  sprintf(towerId, "LAT_Tower%d", tower);
  if (layer%2==0) sprintf(layerViewId, "X%d", layer/2);
  else sprintf(layerViewId, "Y%d", (layer-1)/2);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  
  sprintf(dirname, "TkrNoiseOcc/StripOcc/Tower%d", tower);
  m_svacFile->cd(dirname);

  sprintf(hname, "hTkrStripOccTwr%dLayer%d", tower, layer);
  TH1F *h1 = (TH1F*) gROOT->FindObject(hname);

  
  sprintf(htitle, "Noise Occupancy monitor: %s %s", towerId, layerViewId);
  h1->SetTitle(htitle);
  h1->GetYaxis()->SetTitle("Hit-Strip Occupancy averaged over Layer");
  h1->SetMinimum(1.0e-7);
  h1->SetMaximum(1.0e-2);
  h1->SetStats(0);
  h1->SetMarkerStyle(6);
  h1->Draw("p9");

  sprintf(epsFilePath, "%s/%s_%s_%s_Occ.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  c1->SaveAs(epsFilePath);
  
  sprintf(epsFileName, "%s_%s_%s_Occ.eps", m_prefix, towerId, layerViewId);
  sprintf(gifFileName, "%s_%s_%s_Occ.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_%s_Occ_s.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);


  h1->SetMinimum();
  h1->SetMaximum();
}


void
TkrNoiseRep::drawLayerOcc_perLayer(int tower, int layer){

  char layerViewId[10], towerId[80];
  char dirname[80], hname[80], htitle[80];
  char epsFilePath[250], epsFileName[250], gifFileName[250];
  //char figFileName[200];

  sprintf(towerId, "LAT_Tower%d", tower);
  if (layer%2==0) sprintf(layerViewId, "X%d", layer/2);
  else sprintf(layerViewId, "Y%d", (layer-1)/2);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  
  sprintf(dirname, "TkrNoiseOcc/LayerOcc/Tower%d", tower);
  m_svacFile->cd(dirname);

  sprintf(hname, "hTkrLayerOccTwr%dLayer%d", tower, layer);
  TH1F *h1 = (TH1F*) gROOT->FindObject(hname);

  
  sprintf(htitle, "Noise Occupancy monitor: %s %s", towerId, layerViewId);
  h1->SetTitle(htitle);
  h1->GetYaxis()->SetTitle("Layer Occupancy");
  h1->SetMinimum(1.0e-4);
  h1->SetMaximum(1.0e-0);
  h1->SetStats(0);
  h1->SetMarkerStyle(6);
  h1->Draw("p9");

  sprintf(epsFilePath, "%s/%s_%s_%s_LyrOcc.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  c1->SaveAs(epsFilePath);
  
  sprintf(epsFileName, "%s_%s_%s_LyrOcc.eps", m_prefix, towerId, layerViewId);
  sprintf(gifFileName, "%s_%s_%s_LyrOcc.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_%s_LyrOcc_s.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);
  
  h1->SetMinimum();
  h1->SetMaximum();
}

void
TkrNoiseRep::drawStripOccGr_perLayer(int tower, int layer){

  char layerViewId[10], towerId[80];
  char dirname[80], hname[80], htitle[80];
  int    ix, nx;
  double *vx; 
  double *vy;
  char epsFilePath[250], epsFileName[250], gifFileName[250];

  sprintf(towerId, "LAT_Tower%d", tower);
  if (layer%2==0) sprintf(layerViewId, "X%d", layer/2);
  else sprintf(layerViewId, "Y%d", (layer-1)/2);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  gStyle->SetOptStat(0);

  m_svacFile->cd("TkrNoiseOcc");
  TH1F *hEvtTime = (TH1F*) gROOT->FindObject("hEvtTime");
  nx = (hEvtTime->GetNbinsX());
  if ( (hEvtTime->GetBinContent(nx))==0.0 ) nx = nx-1;

  sprintf(dirname, "TkrNoiseOcc/StripOcc/Tower%d", tower);
  m_svacFile->cd(dirname);

  sprintf(hname, "hTkrStripOccTwr%dLayer%d", tower, layer);
  TH1F *hStripOcc = (TH1F*) gROOT->FindObject(hname);

  vx = new double[nx];
  vy = new double[nx];
  
  for (ix=0; ix<nx; ix++) {
    vx[ix] = (double)hEvtTime->GetBinContent(ix+1)-m_start_time;
    vy[ix] = (double)hStripOcc->GetBinContent(ix+1);
  }
  
  TGraph *gr = new TGraph(nx, vx, vy);
  gr->SetMinimum(1.0e-7);
  gr->SetMaximum(1.0e-2);
  
  gr->SetMarkerStyle(7);
  gr->SetMarkerColor(2);
  gr->SetLineColor(2);

  sprintf(htitle, "Noise Occupancy monitor: %s %s", towerId, layerViewId);
  gr->SetTitle(htitle);
  gr->GetYaxis()->SetTitle("Hit-Strip Occupancy averaged over Layer");
  //gr->SetTitle(hStripOcc->GetTitle());
  gr->GetXaxis()->SetTitle(m_timeAxTitle);
  gr->Draw("AP");
  //sprintf(epsFilePath, "%s/%s_TowerAveOccGr_Tower%d.eps", m_reportDirName, m_prefix, tower);
  sprintf(epsFilePath, "%s/%s_%s_%s_OccGr.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  c1->SaveAs(epsFilePath);

  sprintf(epsFileName, "%s_%s_%s_OccGr.eps", m_prefix, towerId, layerViewId);
  sprintf(gifFileName, "%s_%s_%s_OccGr.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_%s_OccGr_s.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);

}



void
TkrNoiseRep::drawLayerOccGr_perLayer(int tower, int layer){

  char layerViewId[10], towerId[80];
  char dirname[80], hname[80], htitle[80];
  int    ix, nx;
  double *vx; 
  double *vy;
  char epsFilePath[250], epsFileName[250], gifFileName[250];

  sprintf(towerId, "LAT_Tower%d", tower);
  if (layer%2==0) sprintf(layerViewId, "X%d", layer/2);
  else sprintf(layerViewId, "Y%d", (layer-1)/2);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  gStyle->SetOptStat(0);

  m_svacFile->cd("TkrNoiseOcc");
  TH1F *hEvtTime = (TH1F*) gROOT->FindObject("hEvtTime");
  nx = (hEvtTime->GetNbinsX());
  if ( (hEvtTime->GetBinContent(nx))==0.0 ) nx = nx-1;

  sprintf(dirname, "TkrNoiseOcc/LayerOcc/Tower%d", tower);
  m_svacFile->cd(dirname);

  sprintf(hname, "hTkrLayerOccTwr%dLayer%d", tower, layer);
  TH1F *hLayerOcc = (TH1F*) gROOT->FindObject(hname);

  vx = new double[nx];
  vy = new double[nx];
  
  for (ix=0; ix<nx; ix++) {
    vx[ix] = (double)hEvtTime->GetBinContent(ix+1)-m_start_time;
    vy[ix] = (double)hLayerOcc->GetBinContent(ix+1);
  }
  
  TGraph *gr = new TGraph(nx, vx, vy);
  gr->SetMinimum(1.0e-4);
  gr->SetMaximum(1.0e-0);
  
  gr->SetMarkerStyle(7);
  gr->SetMarkerColor(2);
  gr->SetLineColor(2);

  sprintf(htitle, "Noise Occupancy monitor: %s %s", towerId, layerViewId);
  gr->SetTitle(htitle);
  gr->GetYaxis()->SetTitle("Layer Occupancy");
  //gr->SetTitle(hStripOcc->GetTitle());
  gr->GetXaxis()->SetTitle(m_timeAxTitle);
  gr->Draw("AP");
  //sprintf(epsFilePath, "%s/%s_TowerAveOccGr_Tower%d.eps", m_reportDirName, m_prefix, tower);
  sprintf(epsFilePath, "%s/%s_%s_%s_LyrOccGr.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  c1->SaveAs(epsFilePath);

  sprintf(epsFileName, "%s_%s_%s_LyrOccGr.eps", m_prefix, towerId, layerViewId);
  sprintf(gifFileName, "%s_%s_%s_LyrOccGr.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_%s_LyrOccGr_s.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);

}


void
TkrNoiseRep::drawStripHist(int tower, int layer){

  char layerViewId[10], towerId[80];
  char dirname[80], hname[80], htitle[80];
  char epsFilePath[250], epsFileName[250], gifFileName[250];
  double expos;

  sprintf(towerId, "LAT_Tower%d", tower);
  if (layer%2==0) sprintf(layerViewId, "X%d", layer/2);
  else sprintf(layerViewId, "Y%d", (layer-1)/2);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  
  //Expos
  sprintf(dirname, "TkrNoiseOcc/Exposure/Tower%d", tower);
  m_svacFile->cd(dirname);

  sprintf(hname, "hTkrExposTwr%dLayer%d", tower, layer);
  //TH1F *h_exp = (TH1F*) m_svacFile->FindObjectAny(hname);      
  TH1F *h_exp = (TH1F*) gROOT->FindObject(hname);      
  if(h_exp==NULL) {
    sprintf(hname, "hTkrExposTwr%dbiLayer%d", tower, int(layer/2));
    //h_exp = (TH1F*) m_svacFile->FindObjectAny(hname);      
    h_exp = (TH1F*) gROOT->FindObject(hname);      
  }
  expos = h_exp->GetSum();

  //Hitmap
  sprintf(dirname, "TkrNoiseOcc/Hitmap/Tower%d", tower);
  m_svacFile->cd(dirname);

  sprintf(hname, "hTkrHitMapTwr%dLayer%d", tower, layer);
  TH1F *h1 = (TH1F*) gROOT->FindObject(hname);
  h1->Scale(1.0/expos);
  sprintf(htitle, "Noise Occupancy per strip: %s %s %s", towerId, layerViewId, m_prefix);
  h1->SetTitle(htitle);
  h1->GetXaxis()->SetTitle("Strip ID");
  h1->GetYaxis()->SetTitle("Noise Occupancy");
  h1->SetStats(0);
  //h1->SetMarkerStyle(6);
  h1->Draw();

  sprintf(epsFilePath, "%s/%s_%s_%s_Strip.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  c1->SaveAs(epsFilePath);
  
  sprintf(epsFileName, "%s_%s_%s_Strip.eps", m_prefix, towerId, layerViewId);
  sprintf(gifFileName, "%s_%s_%s_Strip.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_%s_Strip_s.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);
  
  //sprintf(figFileName, "%s/%s_%s_%s_LyrOcc.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  //c1->SaveAs(figFileName);

}

void
TkrNoiseRep::drawMultiHist(int tower, int layer){

  char layerViewId[10], towerId[80];
  char dirname[80], hname[80], htitle[80];
  char epsFilePath[250], epsFileName[250], gifFileName[250];

  sprintf(towerId, "LAT_Tower%d", tower);
  if (layer%2==0) sprintf(layerViewId, "X%d", layer/2);
  else sprintf(layerViewId, "Y%d", (layer-1)/2);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  
  //Multi
  sprintf(dirname, "TkrNoiseOcc/Multi/Tower%d", tower);
  m_svacFile->cd(dirname);

  sprintf(hname, "hTkrNoiseMulTwr%dLayer%d", tower, layer);
  TH1F *h1 = (TH1F*) gROOT->FindObject(hname);
  sprintf(htitle, "Hit Multiplicity: %s %s %s", towerId, layerViewId, m_prefix);
  h1->SetTitle(htitle);
  h1->GetXaxis()->SetTitle("Number of hits per event");
  h1->GetYaxis()->SetTitle("Number of events");
  h1->SetStats(0);
  //h1->SetMarkerStyle(6);
  h1->Draw();
  
  sprintf(epsFilePath, "%s/%s_%s_%s_Multi.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  c1->SaveAs(epsFilePath);
  
  sprintf(epsFileName, "%s_%s_%s_Multi.eps", m_prefix, towerId, layerViewId);
  sprintf(gifFileName, "%s_%s_%s_Multi.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_%s_Multi_s.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);
  
}

void
TkrNoiseRep::drawTotHist(int tower, int layer){

  char layerViewId[10], towerId[80];
  char dirname[80], hname[80], htitle[80];
  char epsFilePath[250], epsFileName[250], gifFileName[250];

  sprintf(towerId, "LAT_Tower%d", tower);
  if (layer%2==0) sprintf(layerViewId, "X%d", layer/2);
  else sprintf(layerViewId, "Y%d", (layer-1)/2);

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 400);
  c1->SetLogy();
  
  //TOT
  sprintf(dirname, "TkrNoiseOcc/TOT/Tower%d", tower);
  m_svacFile->cd(dirname);

  //ToT0
  sprintf(hname, "hTkrNoiseTot0Twr%dLayer%d", tower, layer);
  TH1F *h1 = (TH1F*) gROOT->FindObject(hname);
  sprintf(htitle, "Noise Hit ToT0: %s %s %s", towerId, layerViewId, m_prefix);
  h1->SetTitle(htitle);
  h1->GetXaxis()->SetTitle("ToT");
  h1->GetYaxis()->SetTitle("Number of events");
  h1->SetStats(0);
  //h1->SetMarkerStyle(6);
  h1->Draw();
  
  sprintf(epsFilePath, "%s/%s_%s_%s_ToT0.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  c1->SaveAs(epsFilePath);
  
  sprintf(epsFileName, "%s_%s_%s_ToT0.eps", m_prefix, towerId, layerViewId);
  sprintf(gifFileName, "%s_%s_%s_ToT0.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_%s_ToT0_s.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);
  
  //ToT1
  sprintf(hname, "hTkrNoiseTot1Twr%dLayer%d", tower, layer);
  TH1F *h2 = (TH1F*) gROOT->FindObject(hname);
  sprintf(htitle, "Noise Hit ToT1: %s %s %s", towerId, layerViewId, m_prefix);
  h2->SetTitle(htitle);
  h2->GetXaxis()->SetTitle("ToT");
  h2->GetYaxis()->SetTitle("Number of events");
  h2->SetStats(0);
  //h1->SetMarkerStyle(6);
  h2->Draw();
  
  sprintf(epsFilePath, "%s/%s_%s_%s_ToT1.eps", m_reportDirName, m_prefix, towerId, layerViewId);
  c1->SaveAs(epsFilePath);
  
  sprintf(epsFileName, "%s_%s_%s_ToT1.eps", m_prefix, towerId, layerViewId);
  sprintf(gifFileName, "%s_%s_%s_ToT1.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 750, 500);
  sprintf(gifFileName, "%s_%s_%s_ToT1_s.gif", m_prefix, towerId, layerViewId);
  epstogif(m_reportDirName, epsFileName, gifFileName, 300, 200);
  
}

void
TkrNoiseRep::getTwrOccList(){

  int tower;
  char hname[80];

  for(tower=0; tower<16; tower++) {

    // Maximum
    m_svacFile->cd("TkrNoiseOcc/TowerAveOcc");
    sprintf(hname, "hTkrTowerAveOccTwr%d", tower); 
    TH1F *h_towerAveOcc = (TH1F*) gROOT->FindObject(hname);
    m_maxTwrOccList[tower] = h_towerAveOcc->GetMaximum();

    // Mean (not accurate)
    m_aveTwrOccList[tower] = h_towerAveOcc->GetSum()/h_towerAveOcc->GetNbinsX();
    //if (m_maxTwrOccList[tower]>m_critical_strip_occ) {
    //  m_test_status = 0;
    //}
    
  }
}

void
TkrNoiseRep::getTimeStr(){

  //double start_time, end_time;
  time_t gmt_start_time;
  int    offset_time;
  int    nx;

  m_svacFile->cd("TkrNoiseOcc");
  TH1F *hEvtTime = (TH1F*) gROOT->FindObject("hEvtTime");
  
  offset_time = (365*31+8)*(24*3600);

  m_start_time = (double)hEvtTime->GetBinContent(1);
  nx = hEvtTime->GetNbinsX();
  m_end_time = (double)hEvtTime->GetBinContent(nx);
  if (m_end_time==0.0) m_end_time = (double)hEvtTime->GetBinContent(nx-1);
  sprintf(m_durationStr, "%.0f", m_end_time-m_start_time);

  gmt_start_time = (int)m_start_time+offset_time;
  strftime(m_startTimeStr, 79, "%Y-%b-%d %H:%M:%S", gmtime(&gmt_start_time));
  strftime(m_timeAxTitle, 79, "Time [s] from %Y %b %d %H:%M:%S", gmtime(&gmt_start_time));
}


void
TkrNoiseRep::epstogif(const char *figFileDir, const char *epsFileName, const char *gifFileName, int xsize, int ysize){

  //int xsize, ysize;
  char cmdstr[300];

  std::string pwd = gSystem->WorkingDirectory();
  gSystem->cd(figFileDir);
  
  sprintf(cmdstr, "pstopnm -ppm -xborder 0 -yborder 0 -portrait -xsize %d -ysize %d %s", xsize, ysize, epsFileName);
  gSystem->Exec(cmdstr);
  sprintf(cmdstr, "ppmtogif %s001.ppm > %s", epsFileName, gifFileName);
  gSystem->Exec(cmdstr);
  sprintf(cmdstr, "rm -f %s001.ppm", epsFileName);
  gSystem->Exec(cmdstr);

  gSystem->cd(pwd.c_str());
}

void
TkrNoiseRep::generateSummary(){

  int tower;
  
  getTimeStr();
  makeAncRootFile();
  getTwrOccList();
  findNoisyLayer();

  for (tower=0; tower<16; tower++) {
    drawStripOccGr_TowerAve(tower);
  }
  drawStripOccGr_LatAve();

  drawStripOccAve();
  drawStripOccMax();
  drawLayerOccAve();
  drawLayerOccMax();
  drawLayerOccMax_sum();
  drawLargeMultiRatio();
  drawTowerAveStripOccHist();

  generateSummaryPage();
}

void
TkrNoiseRep::generateSummaryPage(){

  //int  status_pass;
  int  tower, bilayer, view;
  char xychar[2] = {'X', 'Y'};
  char font_color[20];
  char layerRep[100];
  char repFilename[250];

  getTwrOccList();

  sprintf(repFilename, "%s/%s_report.html", m_reportDirName, m_prefix);
  
  FILE *outfile = fopen(repFilename, "w");
  
  fprintf(outfile, "<html> <head>\n");
  fprintf(outfile, "<title>LAT TKR noise occupancy summary RunID=%s</title>\n", m_prefix);
  fprintf(outfile, "<style TYPE=\"text/css\">\n");
  fprintf(outfile, "<!--\n");
  //fprintf(outfile, "BODY,TH,TD {font-size:10pt;}\n");
  //fprintf(outfile, "DT,DD {font-size:12pt;}\n");
  fprintf(outfile, "DT {font-weight:bold;}\n");
  fprintf(outfile, "-->\n");
  fprintf(outfile, "</style>\n");
  fprintf(outfile, "</head>\n");
  fprintf(outfile, "\n");
  fprintf(outfile, "<body>\n");
  fprintf(outfile, "<h1>LAT TKR noise occupancy summary RunID=%s</h1>\n", m_prefix);
  fprintf(outfile, "\n");
  fprintf(outfile, "<dl>\n");
  fprintf(outfile, "<dt>RunID</dt><dd>%s</dd>\n", m_prefix);
  fprintf(outfile, "<dt>Start Time</dt><dd>%s</dd>\n", m_startTimeStr);
  fprintf(outfile, "<dt>Duraion [s]</dt><dd>%s</dd>\n", m_durationStr);
  fprintf(outfile, "<dt>Number of noise-sample triggers (minimum in the entire layers)</dt><dd>%.0f</dd>\n", m_minExpos);

  ///// Status
  fprintf(outfile, "<dt>Pass or Fail (LAT-average strip occupancy<%1.e)</dt>\n", m_fail_strip_occ);
  if (m_test_status == 1) {
    fprintf(outfile, "<dd><font color=blue>Pass</font></dd>\n");
  } else if (m_test_status == 2) {
    fprintf(outfile, "<dd><font color=green>Pass*</font></dd>\n");
  } else {
    fprintf(outfile, "<dd><font color=red>Fail</font></dd>\n");
  }

  ///// Warning Message
  // LAT average strip occupancy
  fprintf(outfile, "<dt>LAT-average strip occupancy </dt>\n");
  if ( m_fail_strip_occ < m_latAveOccMean ) {
    fprintf(outfile, "<dd><font color=red> %.3e </font></dd>\n", m_latAveOccMean);
  } else if ( m_critical_strip_occ < m_latAveOccMean ) {
    fprintf(outfile, "<dd><font color=green> %.3e </font></dd>\n", m_latAveOccMean);
  } else {
    fprintf(outfile, "<dd><font color=blue> %.3e </font></dd>\n", m_latAveOccMean);
  }
  
  // Time average layer occupancy
  fprintf(outfile, "<dt>Noisy Layers with average layer occupancy > %.1e </dt>\n", m_critical_layer_occ);
  fprintf(outfile, "<dd><font color=green>");
  for(tower=0; tower<16; tower++) {
    for(bilayer=0; bilayer<18; bilayer++) {
      for(view=0; view<2; view++) {
        if (m_failAveLayerOcc[tower][bilayer][view]==1) {
          fprintf(outfile, "Bay#%d-%c%d:", tower, xychar[view], bilayer); 
        }
      }
    }
  }
  fprintf(outfile, "</font></dd>\n");
  
  ///// Marginal
  fprintf(outfile, "<dt>(Marginal) Noisy Layers </dt>\n");
  fprintf(outfile, "<dd>\n");
  fprintf(outfile, "<dl>\n");
  fprintf(outfile, "<dt>Strip occuapncy > %.1e temporarily</dt>\n", m_critical_strip_occ);
  fprintf(outfile, "<dd>\n");
  for(tower=0; tower<16; tower++) {
    for(bilayer=0; bilayer<18; bilayer++) {
      for(view=0; view<2; view++) {
        if (m_failStripOcc[tower][bilayer][view]==1) {
          sprintf(layerRep, "%s_LAT_Tower%d_%c%d_report.html", m_prefix, tower, xychar[view], bilayer);
          fprintf(outfile, "<a href=\"%s\">Bay#%d-%c%d</a>\n", layerRep, tower, xychar[view], bilayer);
        }
      }
    }
  }
  fprintf(outfile, "</dd>\n");
  fprintf(outfile, "<dt>Layer occuapncy > %.1e temporarily</dt>\n", m_critical_layer_occ);
  fprintf(outfile, "<dd>\n");
  for(tower=0; tower<16; tower++) {
    for(bilayer=0; bilayer<18; bilayer++) {
      for(view=0; view<2; view++) {
        if (m_failLayerOcc[tower][bilayer][view]==1) {
          sprintf(layerRep, "%s_LAT_Tower%d_%c%d_report.html", m_prefix, tower, xychar[view], bilayer);
          fprintf(outfile, "<a href=\"%s\">Bay#%d-%c%d</a>\n", layerRep, tower, xychar[view], bilayer);
        }
      }
    }
  }
  fprintf(outfile, "</dd>\n");
  fprintf(outfile, "<dt>Ratio of large multiplicy(>%d) noise > %.1e </dt>\n", m_critical_multi, m_critical_multi_ratio);
  fprintf(outfile, "<dd>\n");
  for(tower=0; tower<16; tower++) {
    for(bilayer=0; bilayer<18; bilayer++) {
      for(view=0; view<2; view++) {
        if (m_failMultiRatio[tower][bilayer][view]==1) {
          sprintf(layerRep, "%s_LAT_Tower%d_%c%d_report.html", m_prefix, tower, xychar[view], bilayer);
          fprintf(outfile, "<a href=\"%s\">Bay#%d-%c%d</a>\n", layerRep, tower, xychar[view], bilayer);
        }
      }
    }
  }
  fprintf(outfile, "</dd>\n");

  fprintf(outfile, "</dl>\n");

  fprintf(outfile, "</ul>\n");
  fprintf(outfile, "</dd>\n");
  fprintf(outfile, "</dl>\n");
  fprintf(outfile, "<hr>\n");
  fprintf(outfile, "\n");

    
  fprintf(outfile, "<ul>\n");

  fprintf(outfile, "<li><h3> Summary </h3>\n");

  fprintf(outfile, "<ul>\n");
  fprintf(outfile, "<li><h4> Time Variation of the LAT-average noise occupancy </h4>\n");
  fprintf(outfile, "<table border=\"1\">\n");
  fprintf(outfile, "<tr>\n");
  fprintf(outfile, "<th>Time Variation</th>\n");
  fprintf(outfile, "<th>Mean</th>\n");
  fprintf(outfile, "<th>Maximum</th>\n");
  fprintf(outfile, "</tr>\n");
  fprintf(outfile, "<tr>\n");
  fprintf(outfile, "<td>\n");
  fprintf(outfile, "<a href=\"%s_LatAveOccGr.gif\"><img src=\"%s_LatAveOccGr_s.gif\"></a>\n", m_prefix, m_prefix);
  fprintf(outfile, "</td>\n");
  if (m_latAveOccMean>m_critical_strip_occ) sprintf(font_color, "red");
  else sprintf(font_color, "blue");
  fprintf(outfile, "<td><font color=%s> %.4e</font></td>\n", font_color, m_latAveOccMean);
  if (m_latAveOccMean>m_critical_strip_occ) sprintf(font_color, "red");
  else sprintf(font_color, "blue");
  fprintf(outfile, "<td><font color=%s>%.4e</font></td>\n", font_color, m_latAveOccMax);
  fprintf(outfile, "</tr>\n");
  fprintf(outfile, "</table>\n");
  fprintf(outfile, "<li><h4> Noise occupancy time average and transient maximum per Tower </h4>\n");
  fprintf(outfile, "<a href=\"%s_TowerAveStripOcc.gif\"><img src=\"%s_TowerAveStripOcc_s.gif\"></a>\n", m_prefix, m_prefix);
  //Layer Occupancy warning
  fprintf(outfile, "<li><h4> Layer Occupanacy time average > %.1e </h4>\n", m_critical_layer_occ);
  fprintf(outfile, "<a href=\"%s_hAveLayerOcc_sum.gif\"><img src=\"%s_hAveLayerOcc_sum_s.gif\"></a>\n", m_prefix, m_prefix);
  fprintf(outfile, "</ul>\n");
  
  ///// Properties per layer
  fprintf(outfile, "<li><h3> Layer profile of Strip Occupancy and Layer Occupancy </h3>\n");
    
  fprintf(outfile, "<ul>\n");
  fprintf(outfile, "  <li><h4> Strip Occupancy </h4>\n");
  fprintf(outfile, "  <table border=\"1\">\n");
  fprintf(outfile, "    <tr>\n");
  fprintf(outfile, "      <th>Layer-average Strip Occupancy <br> Average over the run</th>\n");
  fprintf(outfile, "      <th>Layer-average Strip Occupancy <br> Transient Maximum (&Delta;T= 2 s)</th>\n");
  fprintf(outfile, "    </tr>\n");
  fprintf(outfile, "    <tr>\n");
  fprintf(outfile, "      <td><a href=\"%s_hAveStripOcc.gif\"><img src=\"%s_hAveStripOcc_s.gif\"></a></td>\n", m_prefix, m_prefix);
  fprintf(outfile, "      <td><a href=\"%s_hMaxStripOcc.gif\"><img src=\"%s_hMaxStripOcc_s.gif\"></a></td>\n", m_prefix, m_prefix);
  fprintf(outfile, "    </tr>\n");
  fprintf(outfile, "  </table>\n");
  fprintf(outfile, "\n");
  fprintf(outfile, "  <li><h4> Layer Occupancy </h4>\n");
  fprintf(outfile, "  <table border=\"1\">\n");
  fprintf(outfile, "    <tr>\n");
  fprintf(outfile, "      <th>Layer Occupancy <br> Average over the run</th>\n");
  fprintf(outfile, "      <th>Layer Occupancy <br> Transient Maximum (&Delta;T= 2 s)</th>\n");
  fprintf(outfile, "    </tr>\n");
  fprintf(outfile, "    <tr>\n");
  fprintf(outfile, "      <td><a href=\"%s_hAveLayerOcc.gif\"><img src=\"%s_hAveLayerOcc_s.gif\"></a></td>\n", m_prefix, m_prefix);
  fprintf(outfile, "      <td><a href=\"%s_hMaxLayerOcc.gif\"><img src=\"%s_hMaxLayerOcc_s.gif\"></a></td>\n", m_prefix, m_prefix);
  fprintf(outfile, "    </tr>\n");
  fprintf(outfile, "  </table>\n");
  fprintf(outfile, "</ul>\n");
  
  fprintf(outfile, "<li><h3> Noise strip multiplicity</h3>\n");
  fprintf(outfile, "<ul>\n");
  fprintf(outfile, "<li><h4> Fraction of noise events with strip multiplicity > %d (suspicious as noise-flare)</h4>\n", m_critical_multi);
  fprintf(outfile, "<a href=\"%s_hLargeMultiRatio.gif\"><img src=\"%s_hLargeMultiRatio_s.gif\"></a></td>\n", m_prefix, m_prefix);
  fprintf(outfile, "</ul>\n");
    
  ///// Time variation per tower
  fprintf(outfile, "<li><h3> Time variation of noise occupancy per Tower</h3>\n");
  fprintf(outfile, "<table border=1>\n");
  fprintf(outfile, "<tr>\n");
  fprintf(outfile, "<th>Tower</th>\n");
  fprintf(outfile, "<th>Time Variation</th>\n");
  fprintf(outfile, "<th>Mean</th>\n");
  fprintf(outfile, "<th>Maximum</th>\n");
  fprintf(outfile, "</tr>\n");

  for(tower=0; tower<16; tower++) {
    fprintf(outfile, "<tr>\n");
    fprintf(outfile, "<th>Bay#%d</th>\n",tower);
    fprintf(outfile, "<td><a href=\"%s_TowerAveOccGr_Tower%d.gif\"><img src=\"%s_TowerAveOccGr_Tower%d_s.gif\"></a></td>\n", m_prefix, tower, m_prefix, tower);
    if (m_aveTwrOccList[tower]>m_critical_strip_occ) sprintf(font_color, "red");
    else sprintf(font_color, "blue");
    fprintf(outfile, "<td><font color=%s> %.4e </font></td>\n", font_color, m_aveTwrOccList[tower]);
    if (m_maxTwrOccList[tower]>m_critical_strip_occ) sprintf(font_color, "red");
    else sprintf(font_color, "blue");
    fprintf(outfile, "<td><font color=%s> %.4e </font></td>\n", font_color, m_maxTwrOccList[tower]);

    fprintf(outfile, "</tr>\n");
  } 
  fprintf(outfile, "</table>\n");

  fprintf(outfile, "</ul>\n");
    
  fprintf(outfile, "<hr>\n");
  fprintf(outfile, "<address></address>\n");
  fprintf(outfile, "</body> </html>\n");

  fclose(outfile);

}

void
TkrNoiseRep::generateNoisyLayerReport(){
  int tower, bilayer, view;
  int layer;

  for(tower=0; tower<16; tower++) {
    for(bilayer=0; bilayer<18; bilayer++) {
      for(view=0; view<2; view++) {
        if( (m_failStripOcc[tower][bilayer][view]==1) 
            || (m_failLayerOcc[tower][bilayer][view]==1) 
            || (m_failMultiRatio[tower][bilayer][view]==1) ){

          layer = bilayer*2+view;
          drawStripOccGr_perLayer(tower, layer);
          drawLayerOccGr_perLayer(tower, layer);
          drawStripHist(tower, layer);
          drawMultiHist(tower, layer);
          drawTotHist(tower, layer);
          generateLayerReport(tower, layer);
        }
      }
    }
  }
}

void
TkrNoiseRep::generateLayerReport(int tower, int layer){
 
  char tkrId[20], layerId[20];
  char repFilename[250];

  sprintf(tkrId, "LAT_Tower%d", tower);
  if (layer%2==0) sprintf(layerId, "X%d", layer/2);
  else sprintf(layerId, "Y%d", (layer-1)/2);

  sprintf(repFilename, "%s/%s_%s_%s_report.html", m_reportDirName, m_prefix, tkrId, layerId);
  FILE *outfile = fopen(repFilename, "w");

  fprintf(outfile,"<html> <head>\n");
  fprintf(outfile,"<title>%s-%s %s</title>\n",tkrId,layerId,m_prefix);
  fprintf(outfile,"</head>\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"<body>\n");
  fprintf(outfile,"<h1>%s-%s %s</h1>\n", tkrId,layerId,m_prefix);
  fprintf(outfile,"\n");
  fprintf(outfile,"<hr>\n");
  fprintf(outfile,"<table>\n");
    
  fprintf(outfile,"  <tr>\n");
  fprintf(outfile,"    <th>Strip Occupancy Time history</th>\n");
  fprintf(outfile,"    <th>Layer Occupancy Time history</th>\n");
  fprintf(outfile,"  </tr>\n");
  fprintf(outfile,"  <tr>\n");
  fprintf(outfile,"    <td><a href=\"%s_%s_%s_OccGr.gif\"><img src=\"%s_%s_%s_OccGr_s.gif\"></a></td>\n",       m_prefix,tkrId,layerId,m_prefix,tkrId,layerId);
  fprintf(outfile,"    <td><a href=\"%s_%s_%s_LyrOccGr.gif\"><img src=\"%s_%s_%s_LyrOccGr_s.gif\"></a></td>\n", m_prefix,tkrId,layerId,m_prefix,tkrId,layerId);
  fprintf(outfile,"  </tr>\n");
    
  fprintf(outfile,"  <tr>\n");
  fprintf(outfile,"    <th>Strip Profile</th>\n");
  fprintf(outfile,"    <th>Hit Multiplicity</th>\n");
  fprintf(outfile,"  </tr>\n");
  fprintf(outfile,"  <tr>\n");
  fprintf(outfile,"    <td><a href=\"%s_%s_%s_Strip.gif\"><img src=\"%s_%s_%s_Strip_s.gif\"></a></td>\n",m_prefix,tkrId,layerId,m_prefix,tkrId,layerId);
  fprintf(outfile,"    <td><a href=\"%s_%s_%s_Multi.gif\"><img src=\"%s_%s_%s_Multi_s.gif\"></a></td>\n",m_prefix,tkrId,layerId,m_prefix,tkrId,layerId);
  fprintf(outfile,"  </tr>\n");
  
  fprintf(outfile,"  <tr>\n");
  fprintf(outfile,"    <th>ToT Low side</th>\n");
  fprintf(outfile,"    <th>ToT High side</th>\n");
  fprintf(outfile,"  </tr>\n");
  fprintf(outfile,"  <tr>\n");
  fprintf(outfile,"    <td><a href=\"%s_%s_%s_ToT0.gif\"><img src=\"%s_%s_%s_ToT0_s.gif\"></a></td>\n",m_prefix,tkrId,layerId,m_prefix,tkrId,layerId);
  fprintf(outfile,"    <td><a href=\"%s_%s_%s_ToT1.gif\"><img src=\"%s_%s_%s_ToT1_s.gif\"></a></td>\n",m_prefix,tkrId,layerId,m_prefix,tkrId,layerId);
  fprintf(outfile,"  </tr>\n");
    
  fprintf(outfile,"</table>\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"\n");
  fprintf(outfile,"<hr>\n");
  fprintf(outfile,"<address></address>\n");
  fprintf(outfile,"</body> </html>\n");
    
  fclose(outfile);
}



