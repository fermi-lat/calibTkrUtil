import os, sys, array, time, string, math, datetime, glob

sys.path = [( os.path.join( os.environ.get( "ROOTSYS" ), "lib" ) )] + sys.path

import ROOT
import tkrUtils

ROOT.gSystem.Load("tkrPyRoot")
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch() # set ROOT into batch mode
canvas = ROOT.TCanvas( 'TKR', 'TKR monitor', 200, 10, 700, 500 )

#
# initialization of common variables, parameters
# 
nTower = tkrUtils.g_nTowers
nPlane = tkrUtils.g_nUniPlanes
nStrip = tkrUtils.g_nStrips

nsat = 15 # GTRC hit limit + 1

keys = [ ("TOT_Peak", "layer", "D"), ("TOT_LWidth", "layer", "D"), \
         ("TOT_GSigma", "layer", "D"), ("TOT_Entries", "layer", "D" ), \
         ("TOT_PeakError", "layer", "D"), ("TOT_LWidthError", "layer", "D"), \
         ("TOT_GSigmaError", "layer", "D"), ("TOT_FitProb", "layer", "D"), \
         ("TOT_FracLowTOT", "layer", "D" ),\
         ("layerEff", "layer", "D"), ("layerEff_err", "layer", "D"), \
         ("towerEff", "tower", "D"), ("towerEff_err", "tower", "D"), \
         ("layerdXY", "layer", "D"), ("layerdXY_err", "layer", "D"), \
         ("layerOcc", "layer", "D"), ("stripOcc", "layer", "D"), \
         ("fracSat", "layer", "D") ]
nelem = {}
nelem["layer"] = tkrUtils.g_nTowers * tkrUtils.g_nUniPlanes
nelem["tower"] = tkrUtils.g_nTowers
typecode = { "C":"c", "B":"b", "b":"B", "S":"i", "s":"I", "F":"f", "D":"d" }
alertTypes = [ "totPeak", "totFit", "layerEff", "layerOcc", "stripOcc", \
               "fracSat" ]
alertLevels = [ "ERR", "WARN", "NOTE", "INFO" ]
towerHistNames = ["layerEff", "layerOcc", "stripOcc", "fracSat", "layerdXY"]

#**************************************************************************
#**************************************************************************
#
# class TKR monitor
#
#**************************************************************************
#**************************************************************************
class TkrMonitor:
  def __init__(self, iname, oname, htmldir):
    #
    # set up directories and files.
    #
    self.inputRoot = ROOT.TFile( iname )
    print self.inputRoot.GetEndpointUrl()
    self.outputRoot = ROOT.TFile( oname )
    self.alertFile = aname
    self.logName = logname
    self.htmldir = htmldir
    if not os.path.exists( self.htmldir ):
      print "html directory, %s, does not exist. Create a new one." % htmldir
      os.mkdir( self.htmldir )
  
    self.totdir = os.path.join( self.htmldir, "totPlots" )
    if not os.path.exists( self.totdir ):
      print "No TOT plot directory, %s. Create a new one." % self.totdir
      os.mkdir( self.totdir )

    self.refRoot = []
    try:
      self.mondir = os.path.join( os.getenv("LATMonRoot"), "TKR", "TkrMonitor" )
    except:
      self.mondir = os.path.join( os.getenv( "CALIBTKRUTILROOT" ), "python" )

    if not os.path.exists( self.mondir ):
      print "TKR Monitor directory, %s, does not exist." % self.mondir
      sys.exit()
    path = os.path.join( self.mondir, "refRoot", "*ref*.root" )
    refs = glob.glob( path )
    for ref in refs:
      self.refRoot.append( ROOT.TFile( ref ) )
    
    #
    # initialize some parameters
    #
    self.values = {}
    for (key, type, vtype) in keys:
      self.values[key] = array.array( typecode[vtype], [-1]*nelem[type] )

    self.logs = {}
    for type in alertTypes:
      self.logs[type] = {}
      for level in alertLevels:
        self.logs[type][level] = []

    self.hists = {}
    self.alerts = []
    self.towerHists = []
    self.layerHists = []
    self.savedPictures = []
    for tower in range(tkrUtils.g_nTowers):
      self.towerHists.append( {} )
      for name in towerHistNames:
        hname = "%s_T%d" % (name, tower)
        title = "%s in tower %d" % (name, tower)
        self.towerHists[-1][name] = ROOT.TH1F(hname,title,nPlane,0,nPlane)

      self.alerts.append( [] )
      self.layerHists.append( [] )
      self.savedPictures.append( [] )
      for unp in range(tkrUtils.g_nUniPlanes):
        self.alerts[-1].append( [] )
        self.layerHists[-1].append( {} )
        self.savedPictures[-1].append( {} )

    self.hists["layerEff"] = ROOT.TH2F("layerEff_map", \
                                       "Layer efficiency map",\
                                       nTower,0,nTower, nPlane,0,nPlane)
    self.hists["stripOcc"] = ROOT.TH2F("stripOcc_map", \
                                       "Average strip occupancy map",\
                                       nTower,0,nTower, nPlane,0,nPlane)
    self.hists["layerOcc"] = ROOT.TH2F("layerOcc_map", \
                                        "Layer-OR occupancy map",\
                                        nTower,0,nTower, nPlane,0,nPlane)
    self.hists["fracSat"] = ROOT.TH2F("fracSat_map", \
                                      "frac. of saturated events map",\
                                      nTower,0,nTower, nPlane,0,nPlane)
    self.hists["towerEff"] = ROOT.TH1F("towerEff", "tower Efficiency", \
                                       nTower,0,nTower)
    self.hists["towerdX"] = ROOT.TH1F("towerdX", "tower delta X", \
                                      nTower,0,nTower)
    self.hists["towerdY"] = ROOT.TH1F("towerdY", "tower delta Y", \
                                      nTower,0,nTower)
    self.hists["totPeak"] = ROOT.TH1F("totPeakDist",\
                                      "TOT peak distributioin",100,0,10)
      
    self.readParamLimits()

  #
  #*********************
  # read parameter limits from xml file
  #*********************
  #
  def readParamLimitsXml(self):
    self.limits={}
    xname = os.path.join( self.mondir, "paramLimitsDefault.xml"  )
    dom = xml.dom.minidom.parse( xmlfile )
    topElm = dom.getElementsByTagName("ParamLimits")[0]
    tstuple = topElm.getAttribute("timestamp")
    limits = topElm.getElementsByTagName("limit")
    for limit in limits:
      type = str( limit.getAttribute("type") )
      level = str( limit.getAttribute("level") )
      value = float( limit.getAttribute("value") )
      if type not in alterType or name not in altertLevel:
        print "Invalid alert type:%s or level:%s." % (type,level)
      else: self.limits[name] = value
    print self.limits


  #
  #*********************
  # read parameter limits from file
  #*********************
  #
  def readParamLimits(self):
    self.limits={}
    tname = os.path.join( self.mondir, "limits.txt"  )
    input = open(tname, "r")
    for line in input:
      words = string.split(line)
      self.limits[words[0]] = float(words[1])
    print self.limits
    
  
  #
  #*********************
  # get time stamps and run IDs
  #*********************
  #
  def getTimeStamps(self):
    tree = self.inputRoot.FindObjectAny( "timeStamps" )
    startTime = array.array( 'd', [0] )
    endTime = array.array( 'd', [0] )
    firstRunId = array.array( 'L', [0] )
    lastRunId = array.array( 'L', [0] )
    
    tree.SetBranchAddress( "startTime", startTime )
    tree.SetBranchAddress( "endTime", endTime )
    tree.SetBranchAddress( "firstRunId", firstRunId )
    tree.SetBranchAddress( "lastRunId", lastRunId )
    
    self.startTime = array.array( 'd', [-1] )
    self.endTime = array.array( 'd', [1] )
    self.firstRunId = array.array( 'L', [0] )
    self.lastRunId = array.array( 'L', [0] )
    
    for i in range(tree.GetEntries()):
      tree.GetEntry( i )
      if self.startTime[0] < 0 or startTime[0] < self.startTime[0]:
        self.startTime[0] = startTime[0]
      if endTime[0] > self.endTime[0]: self.endTime[0] = endTime[0]
      if self.firstRunId[0] < 1 or firstRunId[0] < self.firstRunId[0]:
        self.firstRunId[0] = firstRunId[0]
      if lastRunId[0] > self.lastRunId[0]: self.lastRunId[0] = lastRunId[0]
      
    if firstRunId[0] == lastRunId[0]: self.runID = "%d" % firstRunId[0]
    else: self.runID = "%d-%d" % ( firstRunId[0], lastRunId[0] )

    print self.startTime, self.endTime, self.runID

  #
  #*********************
  # pefoerm analysis
  #*********************
  #
  def analyzeTKR(self):
    self.ffit = ROOT.defLangau( "langau", 0, 30 )
    self.ffit.SetParNames( "LWidth", "MP", "Area", "GSigma" )

    #
    # loop towers/layers ann fit TOT distributions
    #
    for tower in range(tkrUtils.g_nTowers):
      self.analyzeEfficiency( tower )
      for unp in range(tkrUtils.g_nUniPlanes):
      #for unp in range(1):
        self.fitTOT( tower, unp )
        self.analyzeOccupancy( tower, unp )
        
    self.analyzeTotPeak()
    self.analyzeDisplacement()

  #
  #*********************
  # fit TOT distributions 
  #*********************
  #
  def fitTOT(self, tower, unp):
    ielm = unp + tower*tkrUtils.g_nUniPlanes
    lname = tkrUtils.g_layerNames[unp]
    name = "chargeT%d%s" % (tower, lname)
    chargeHist = self.inputRoot.FindObjectAny( name )
    self.layerHists[tower][unp]["chargeHist"] = chargeHist;

    entries = chargeHist.Integral()
    mean = chargeHist.GetMean()
    rms = chargeHist.GetRMS()
    self.values["TOT_Entries"][ielm] = entries

    if entries<200 or mean==0.0 or rms==0.0:
      alert = "Entries %.0f, Mean: %.1f, RMS: %.1f skipped" \
              % (entries, mean, rms)                
      self.logAlerts( tower, unp, alert, "totFit" )
      return

    binWidth = chargeHist.GetBinWidth( 2 )
    bin = int(mean*0.5/binWidth) + 1
    fracBadTot = chargeHist.Integral(1,bin) / chargeHist.Integral()
    self.values["TOT_FracLowTOT"][ielm] = fracBadTot

    lowLim = mean - 1.4 * rms
    if fracBadTot > self.limits["fracBadTot"] and lowLim < mean*0.5:
      lowLim = mean*0.5
      alert = "large bad TOT fraction: %.3f > %.3f" \
              % (fracBadTot,self.limits["fracBadTot"])
      self.logAlerts( tower, unp, alert, "totFit" )

    self.ffit.SetParLimits( 0, 0.10, 0.5 )
    self.ffit.SetParLimits( 1, 0.0, mean*2 )
    self.ffit.SetParLimits( 2, 0.0, entries*1.5 )
    self.ffit.SetParLimits( 3, rms*0.2, rms*1.5 )
    self.ffit.SetParameters( rms*0.2, mean*0.75, entries*0.1, rms*0.4 )
    #self.ffit.SetRange( lowLim, mean+2.5*rms )
    #self.ffit.FixParameter( 4, m_RSigma )
    #self.ffit.FixParameter( 5, m_GFrac )
    chargeHist.Fit( "langau", "RBQ", "", lowLim, mean+2.5*rms )
    canvas.Update()

    #0:width 1:peak 2:total entries 3:width(sigma)
    par = self.ffit.GetParameters()
    error = self.ffit.GetParErrors()
    self.values["TOT_LWidth"][ielm] = par[0]
    self.values["TOT_Peak"][ielm] = par[1]
    self.values["TOT_GSigma"][ielm] = par[3]
    self.values["TOT_LWidthError"][ielm] = error[0]
    self.values["TOT_PeakError"][ielm] = error[1]
    self.values["TOT_GSigmaError"][ielm] = error[3]
    self.values["TOT_FitProb"][ielm] = self.ffit.GetProb()
    self.hists["totPeak"].Fill( par[1] )
    
  #
  #*********************
  # analyze TOT peak distribution and outliers
  #*********************
  #
  def analyzeTotPeak(self):
    totPeakDist = self.hists["totPeak"]
    fgaus = ROOT.TF1("fgaus", "gaus(0)")
    fgaus.SetParameters( totPeakDist.Integral(), \
                         totPeakDist.GetMean(), totPeakDist.GetRMS() )
    fgaus.SetLineColor(2)
    fgaus.SetLineWidth(2)
    totPeakDist.Fit("fgaus")

    mean = fgaus.GetParameter(1)
    rms = fgaus.GetParameter(2)
    limit = rms * self.limits["totPeakSigma"]

    self.totPeakAlerts = []
    for ielm in range(tkrUtils.g_nTowers*tkrUtils.g_nUniPlanes):
      value = self.values["TOT_Peak"][ielm]
      if value > 0 and abs(value-mean) > limit:
        # print ielm, values["TOT_Peak"][ielm]
        tower, unp = divmod( ielm, tkrUtils.g_nUniPlanes )
        alert = "TOT peak outlier: %.1f" % self.values["TOT_Peak"][ielm]
        self.logAlerts( tower, unp, alert, "totPeak", alertLevels[1] )
    

  #
  #*********************
  # analyze layer effciencies
  #*********************
  #
  def analyzeEfficiency(self, tower):
    hlhit = self.inputRoot.FindObjectAny("lhitT"+str(tower))
    hltrk = self.inputRoot.FindObjectAny("ltrkT"+str(tower))
    
    for unp in range(nPlane):
      ielm = unp + tower*nPlane
      lname = tkrUtils.g_layerNames[unp]
      
      #efficiency
      lhit = hlhit.GetBinContent(unp+1)
      ltrk = hltrk.GetBinContent(unp+1)
      eff = lhit / ltrk
      err = eff*(1-eff) / ltrk
      if err > 0.0: err = math.sqrt( err )
      else: err = 1.0 / ltrk
      if eff < self.limits["layerEff"]:
        alert = "layer efficiency, %.1f%s < %.1f%s" \
                % (eff*100, "%", self.limits["layerEff"]*100, "%" )
        self.logAlerts( tower, unp, alert, "layerEff" ,alertLevels[2] )

      key = "layerEff"
      self.values[key][ielm] = eff
      self.values[key+"_err"][ielm] = err
      self.towerHists[tower][key].SetBinContent(unp+1,eff)
      self.towerHists[tower][key].SetBinError(unp+1,err)
      self.hists[key].SetBinContent(tower+1,unp+1,1-eff)
      
    sumltrk = hltrk.Integral()
    if sumltrk > 0:
      teff = hlhit.Integral()/sumltrk
      terr = teff*(1-teff)/sumltrk
    else:
      teff = 0.0
      terr = 1.0E-5
    if terr >0.0:terr = math.sqrt(terr)
    else:terr = 1/sumltrk
    key = "towerEff"
    self.values[key][tower] = teff
    self.values[key+"_err"][tower] = terr
    self.hists[key].SetBinContent(tower+1,teff)
    self.hists[key].SetBinError(tower+1,terr)

    
  #
  #*********************
  # analyze layer displacements
  #*********************
  #
  def analyzeDisplacement(self):
    
    htresPX = self.inputRoot.FindObjectAny("tresProfX")
    htresPY = self.inputRoot.FindObjectAny("tresProfY")
    for tower in range(tkrUtils.g_nTowers):
      hresP = self.inputRoot.FindObjectAny("resProfT"+str(tower))
      for unp in range(nPlane):
        ielm = unp + tower*nPlane
        lname = tkrUtils.g_layerNames[unp]
    
        # displacement
        resP = hresP.GetBinContent(unp+1)
        resP_err = hresP.GetBinError(unp+1)
        self.values["layerdXY"][ielm] = resP
        self.values["layerdXY_err"][ielm] = resP_err
        self.towerHists[tower]["layerdXY"].SetBinContent(unp+1,resP)
        self.towerHists[tower]["layerdXY"].SetBinError(unp+1,resP_err)

      #
      # this parameter is probably useless.
      #
      tresPx = htresPX.GetBinContent(tower+1)
      tresPx_err = htresPX.GetBinError(tower+1)
      tresPy = htresPY.GetBinContent(tower+1)
      tresPy_err = htresPY.GetBinError(tower+1)
      self.hists["towerdX"].SetBinContent(tower+1, tresPx)
      self.hists["towerdX"].SetBinError(tower+1, tresPx_err)
      self.hists["towerdY"].SetBinContent(tower+1, tresPy)
      self.hists["towerdY"].SetBinError(tower+1, tresPy_err)
    
  #
  #*********************
  # analyze occupancies
  #*********************
  #
  def analyzeOccupancy(self, tower, unp):

    ielm = unp + tower*nPlane
    lname = tkrUtils.g_layerNames[unp]

    hmul = self.inputRoot.FindObjectAny( "hTkrNoiseMulT%d%s" %(tower,lname) )
    hmap = self.inputRoot.FindObjectAny( "hTkrHitMapT%d%s" %(tower,lname) )
    lexp = hmul.Integral()
    locc = hmul.Integral(2,128)
    lsat = hmul.Integral(nsat,128)
    socc = hmap.Integral()
    if lexp != None:
      locc /= lexp
      socc /= (lexp*nStrip)
      lsat /= lexp
    else:
      print "ERROR"
      locc = 0
      socc = 0
      lsat = 0
    if locc > self.limits["layerOcc"]:
      alert = "layer Occucpancy, %.1e > %.1e" \
              % (locc,self.limits["layerOcc"])
      self.logAlerts( tower, unp, alert, "layerOcc", alertLevels[0] )
    if socc > self.limits["stripOcc"]:
      alert = "strip Occucpancy, %.1e > %.1e" \
              % (socc,self.limits["stripOcc"])
      self.logAlerts( tower, unp, alert, "stripOcc", alertLevels[1] )
    if lsat > self.limits["fracSat"]:
      alert = "fraction of saturated events, %.1e > %.1e" \
              % (lsat,self.limits["fracSat"])
      self.logAlerts( tower, unp, alert, "fracSat", alertLevels[2] )
      
    self.values["layerOcc"][ielm] = locc
    self.values["stripOcc"][ielm] = socc
    self.values["fracSat"][ielm] = lsat
    self.towerHists[tower]["layerOcc"].SetBinContent(unp+1,locc)
    self.towerHists[tower]["stripOcc"].SetBinContent(unp+1,socc)
    self.towerHists[tower]["fracSat"].SetBinContent(unp+1,lsat)
    self.hists["layerOcc"].SetBinContent(tower+1,unp+1,locc)
    self.hists["stripOcc"].SetBinContent(tower+1,unp+1,socc)
    self.hists["fracSat"].SetBinContent(tower+1,unp+1,lsat)


  #
  #*********************
  # create html report for the given tower/layer
  #*********************
  #
  def htmlLayerReport(self, tower, unp):
    lname = tkrUtils.g_layerNames[unp]
    ldir = "T%d%s" % (tower, lname)
    for key in self.layerHists[tower][unp].keys():
      pname = "%s_T%d%s-%s.png" % (key, tower, lname, self.runID )
      path = os.path.join( htmldir, ldir )
      if not os.path.exists( path ): 
        print "html sub directory, %s, does not exist. Creat a new one." % path
        os.mkdir( path )
      path = os.path.join( path, pname )
      hist = self.layerHists[tower][unp][key]
      hist.Draw()
      canvas.SaveAs( path )
      self.savedPictures[tower][unp][key] = pname

    self.layerHists[tower][unp] = {} # reset to avoid duplicates
    return ldir
    
  #
  #*********************
  # create html report for the given tower
  #*********************
  #
  def htmlTowerReport(self, tower):
    tdir = "Tower%d" % tower
    for key in self.towerHists[tower].keys():
      pname = "%s_T%d-%s.png" % (key, tower, self.runID )
      path = os.path.join( htmldir, tdir )
      if not os.path.exists( path ):
        print "html sub directory, %s, does not exist. Creat a new one." % path
        os.mkdir( path )
      path = os.path.join( path, pname )
      hist = self.towerHists[tower][key]
      hist.Draw()
      if key[-3:] == "Occ" or key[:4] == "frac" : ROOT.gPad.SetLogy( 1 )
      if self.limits.has_key( key ):
        xmin = hist.GetBinLowEdge(1)
        xmax = hist.GetBinLowEdge( hist.GetNbinsX()+1 )
        limit = self.limits[key]
        tl = ROOT.TLine( xmin, limit, xmax, limit )
        tl.SetLineColor( 2 )
        tl.SetLineWidth( 2 )
        tl.Draw()
      canvas.SaveAs( path )
      ROOT.gPad.SetLogy( 0 )
      
    self.towerHists[tower] = {} # reset to avoid duplicates
    return tdir


  #
  #*********************
  # save histrograms and tress to root file
  #*********************
  #
  def saveRoot(self, fname):

    print "save ROOT file:", fname
    tf = ROOT.TFile( fname, "RECREATE")

    #
    # save root files
    #
    for key in self.hists.keys():
      self.hists[key].Write()

    # tower based histograms
    tf.mkdir( "Towers" )
    tf.cd( "Towers" )
    for tower in range(tkrUtils.g_nTowers):
      for key in self.towerHists[tower].keys():
        self.towerHists[tower][key].Write()
        
    tf.cd()

    # layer based histograms
    for key in ["chargeHist", "stripOcc", "stripEff"]:
      tf.mkdir( key )
      tf.cd( key )
      for tower in range(tkrUtils.g_nTowers):
        for unp in range(tkrUtils.g_nUniPlanes):
          if self.layerHists[tower][unp].has_key( key ):
            self.layerHists[tower][unp][key].Write()
      tf.cd()

    #
    # save TTree
    #
    tree = ROOT.TTree( "tkrMonitor", "tkrMonitor" )
    
    for (key, type, vtype) in keys:
      if type == "layer":
        leaflist = "%s[%d][%d]/%s" \
                   % (key, tkrUtils.g_nTowers, tkrUtils.g_nUniPlanes, vtype )
      elif type == "tower":
        leaflist = "%s[%d]/%s" % (key, tkrUtils.g_nTowers, vtype )
      tree.Branch(key, self.values[key], leaflist)

    tree.Branch( "startTime", self.startTime, "startTime/D" )
    tree.Branch( "endTime", self.endTime, "endTime/D" )
    tree.Branch( "firstRunId", self.firstRunId, "firstRunId/i" )
    tree.Branch( "lastRunId", self.lastRunId, "lastRunId/i" )

    tree.Fill()
    tree.Write()
    
    tf.Close()


  #
  #*********************
  # log alerts
  #*********************
  #
  def logAlerts(self, tower, unp, alert, type, level="INFO"):
    if type not in alertTypes:
      print "wrong alert type:", type
      sys.exit()
    if level not in alertLevels:
      print "wrong alert level:", level
      sys.exit()
      
    lname = tkrUtils.g_layerNames[unp]
    print "%s-%s: %s, T%d %s" % (level, type, alert, tower, lname)
    self.logs[type][level].append( (tower, unp, alert) )
    if level=="WARN" or level=="ERR":
      self.alerts[tower][unp].append( (alert, type, level) )


  #
  #*********************
  # save alert files and log files
  #*********************
  #
  def saveReports(self, aname, logname):
    print "log file:", logname
    lfile = open( logname, 'w')
    for level in alertLevels:
      for type in alertTypes:
        for tower, unp, alert in self.logs[type][level]:
          lname = tkrUtils.g_layerNames[unp]
          log = "%s %s: %s, T%d %s\n" \
                % ("%"+level+"%", type, alert, tower, lname)
          lfile.write( log )
          if level != "INFO":
            pass # use this for alert in the future
    lfile.close()

    self.saveHtmlReport()

    

  ###########################################
  #
  #####   html report   #####################
  #
  ###########################################

  def saveHtmlReport(self):

    #### main page
    path = os.path.join( htmldir, "index.html" )
    self.hout = open( path,"w")

    self.hout.write( "<html>\n" )
    self.hout.write( "<head><title>LAT-TKR Monitor Reprot</title><head>\n" )
    self.hout.write( "<body>\n" )
    self.htmlLine( "<B>LAT-TKR Monitor Reprot</B>", 3, "blue" )

    #
    # run info
    #
    self.hout.write( "<ul>" )
    self.htmlLine( "runID: %s" % self.runID )
    self.htmlLine( "intput file: %s" % self.inputRoot.GetName() )
    self.htmlLine( "Start/End Time: %.0f - %.0f" \
                    % (self.startTime[0], self.endTime[0]) )
    self.htmlLine( "Duration: %.0f second" \
                    % (self.endTime[0]-self.startTime[0]) )
    self.hout.write( "</ul>" )

    #
    # summary plots and alerts
    #
    items = { "layerEff":("layer Inefficinecy","COLZ"),\
              "layerOcc":("layer-OR occupancy","COLZ"),\
              "stripOcc":("Average Strip Occupancy","COLZ"),\
              "fracSat":("Fraction of Saturated Events","COLZ"),\
              }
    self.hout.write( '<table style="text-align: left; width: 100%;" ' \
                     + 'border="0" cellpadding="2" cellspacing="2">\n' )
    self.hout.write( '<tbody>\n' )
    for key in items:
      (title, dopt) = items[key]
      self.hout.write( '<tr><td>' )
      self.htmlLine( title )
      hist = self.hists[key]
      hist.Draw( dopt )
      pname = "%s-%s.png" % (hist.GetName(), self.runID )
      path = os.path.join( htmldir, pname )
      canvas.SaveAs( path )
      self.htmlImage( pname )
      self.hout.write( "</td>\n" )
      self.htmlAlerts( key )
      self.hout.write( "</tr>\n" )
    #
    # close out
    #
    self.hout.write( "</tbody>\n</table>\n" )
    self.hout.write( "</body>\n</html>" )
    self.hout.close()
    

  #*********************
  # misc. html functions
  #*********************
  def htmlImage( self, imageName, table=False ):
    txt = '<a href="%s"><img src="%s" width="300"></a>' \
          % (imageName, imageName)
    if table: txt = '<td>%s</td>' % txt
    self.hout.write( txt+"\n" )

  def htmlLine(self, line, fsize=2, color=None ):
    if color == None:
      self.hout.write( "<FONT SIZE=\"+%d\">%s</FONT><BR>\n" % (fsize, line) )
    else:
      self.hout.write( '<FONT SIZE="+%d" color="%s">%s</FONT><BR>\n' \
                       % (fsize, color, line) )

  def htmlLink(self, word, href, color=None):
    txt = '<a href="%s">%s</a>' % (href, word)
    if color == None: return txt
    else: return '<font color="%s">%s</font>' % (color, txt)
    

  #
  #*********************
  # close files
  #*********************
  #
  def htmlAlerts(self, key ):

    self.hout.write( "<td>\n" )
    for level in alertLevels:
      if len(self.logs[key][level])>0:
        if level == alertLevels[0]: color = "red"
        elif level == alertLevels[1]: color = "yellow"
        else: color = None
        self.htmlLine( level, 2, color )
        for tower, unp, alert in self.logs[key][level]:
          lname = tkrUtils.g_layerNames[unp]
          towerRef = self.htmlTowerReport( tower )
          towerLink = self.htmlLink( "Tower %d"%tower, towerRef )
          layerRef = self.htmlLayerReport( tower, unp )
          layerLink = self.htmlLink( "Layer %s"%lname, layerRef )
          self.htmlLine( "%s: %s %s" % (alert, towerLink, layerLink) )
    self.hout.write( "</td>\n" )
  
  #
  #*********************
  # close files
  #*********************
  #
  def close(self):
    self.inputRoot.Close()
    

#**************************************************************************
#**************************************************************************
#
# main
#
#**************************************************************************
#**************************************************************************
if __name__ == '__main__':

  #
  # decode command arguments
  # [input root file] [output root file] [html dir] [alert file] [log file]
  #
  timeStamp = time.strftime("%y%m%d-%H%M", time.gmtime() )
  timeStamp = time.strftime("%y%m%d", time.gmtime() ) # this temporary
  if len(sys.argv) > 1: iname = sys.argv[1]
  else: iname = "/nfs/slac/g/ki/ki16/htajima/GLAST/rootFiles/TKR.root"
  if len(sys.argv) > 2: oname = sys.argv[2]
  else: oname = "tkrMonitor-%s.root" % timeStamp
  if len(sys.argv) > 3: htmldir = sys.argv[3]
  else: htmldir = "./tkrReports-%s" % timeStamp
  if len(sys.argv) > 4: aname = sys.argv[4]
  else: aname = "tkrAlert-%s.xml" % timeStamp
  if len(sys.argv) > 5: logname = sys.argv[5]
  else: logname = "tkrMonitor-%s.log" % timeStamp

  tkrMonitor = TkrMonitor( iname, oname, htmldir )
  tkrMonitor.getTimeStamps()
  tkrMonitor.analyzeTKR()

  tkrMonitor.saveRoot( oname ) # save root before report.
  tkrMonitor.saveReports( aname, logname ) # this should be after saving ROOT

  tkrMonitor.close()
