# import python standard libraries
import os, sys, array, time, string, math, glob, xml.dom.minidom

sys.path = [( os.path.join( os.environ.get( "ROOTSYS" ), "lib" ) )] + sys.path

# import thrid party libraries
import ROOT

# inpot private libraies
import tkrUtils

# get tag and version numbers
__tag__  = "$Name:  $"
__version__  = "$Revision: 1.13 $"
tagv = "%s:%s" % (__tag__.split()[1], __version__.split()[1])

# ROOT initilization
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
nGTFE = tkrUtils.g_nFE
nChannel = tkrUtils.g_nChannels
#nTower = 4
#nPlane = 9

nsat = 15 # GTRC hit limit + 1

keys = [ ("TOT_Peak", "layer", "D"), ("TOT_LWidth", "layer", "D"), \
         ("TOT_GSigma", "layer", "D"), ("TOT_Entries", "layer", "D" ), \
         ("TOT_PeakError", "layer", "D"), ("TOT_LWidthError", "layer", "D"), \
         ("TOT_GSigmaError", "layer", "D"), ("TOT_FitProb", "layer", "D"), \
         ("TOT_FracLowTOT", "layer", "D" ),\
         ("layerEff", "layer", "D"), ("layerEff_err", "layer", "D"), \
         ("towerEff", "tower", "D"), ("towerEff_err", "tower", "D"), \
         ("trigEff", "tower", "D"), ("trigEff_err", "tower", "D"), \
         ("layerdXY", "layer", "D"), ("layerdXY_err", "layer", "D"), \
         ("layerOcc", "layer", "D"), ("stripOcc", "layer", "D"), \
         ("fracSat", "layer", "D") ]

nelem = {}
nelem["layer"] = nTower * nPlane
nelem["tower"] = nTower
typecode = { "C":"c", "B":"b", "b":"B", "S":"i", "s":"I", "F":"f", "D":"d" }
alertTypes = [ "totPeak", "totFit", "totPeakSigma", "fracBadTot", \
               "layerEff", "towerEff", "stripEff", "trigEff", "towerEff",\
               "LatOcc", "towerOcc", "layerOcc", "stripOcc", \
               "fracSat", "layerdXY", "layerdXYSigma", "fracNHstrip","deadGTFE" ]
alertLevels = [ "ERR", "WARN", "NOTE", "INFO" ]

latHistNames = { "layerEff":("Layer averaged hit inefficinecy ", "tower", "layer", "COLZ", "[%]"),\
                 "layerOcc":("Layer-OR occupancy", "tower", "layer", "COLZ", ""),\
                 "stripOcc":("Layer averaged Strip Occupancy", "tower", "layer", "COLZ", ""),\
                 "fracSat":("Fraction of Saturated Events", "tower", "layer", "COLZ", ""),\
                 "totPeak":("Layer averaged TOT peak distribution", "TOT peak [fC]", \
                            "Number of layers", "H", ""),\
                 "towerEff":("Tower averaged hit efficiency", "tower", \
                             "hit efficinency", " ", " [%]"),\
                 "trigEff":("Tower averaged trigger efficiency", "tower", \
                            "trigger efficiency", " ", " [%]"),\
                 "layerdXY":("Layer displacement from reference [micron]", \
                             "tower", "layer", "COLZ", ""),\
                 "fracNHstrip":("Fraction of no-hit strips", "tower", \
                                  "fraction of no-hit strips [%]", "", ""),\
                 "deadGTFE":("Number of Dead GTFE", "tower", \
                                  "number of dead GTFE", "", ""),\
                 }

towerHistNames = {"layerEff":("Layer Averaged Hit Efficiency", "layer", "hit efficiency", " [%]"),\
                  "layerOcc":("Layer-OR Occupancy","layer", "layer-OR occupancy", ""),\
                  "stripOcc":("Layer Averaged Strip Occupancy", "layer", "strip occupancy", ""),\
                  "fracSat":("Fraction of Saturated Events", "layer", \
                             "fraction of saturated events", ""),\
                  "layerdXY":("Layer Displacement from reference", "Layer", \
                              "layer displacement", " [micron]"),\
                  "totPeak":("TOT Peak distribution", "Layer", "TOT peak", " [fC]"),\
                  "stripOccT":("Tower averaged Strip Occupancy, Time-Variation","time[min]", \
                               "Tower averaged strip occupancy", ""),\
                  "fracNHstrip":("Fraction of no-hit strips", "layer", \
                                   "Fraction of no-hit strips", "[%]"),\
                  }

layerHistNames = {"chargeHist":("TOT Charge-deposit distribution", "charge [fC]", "Entries/bin"),\
                  "stripOcc":("Strip Occupancy", "strip", "strip Occupancy"),\
                  "stripEff":("GTFE averaged Strip Efficiency", "strip", \
                              "GTFE averaged strip efficiency [%]"),\
                  "layerOccT":("Layer-OR Occupancy Time-Variation","time[min]", \
                               "layer-OR occupancy", ""),\
                  "stripOccT":("Layer averaged Strip Occupancy, Time-Variation","time[min]", \
                               "layer averaged strip occupancy", ""),\
                  }
resScale = 1000
timebin = 120
red = "red"
yellow = "#e08000"



#**************************************************************************
#**************************************************************************
#
# class TKR monitor
#
#**************************************************************************
#**************************************************************************
class TkrMonitor:
  def __init__(self, iname, htmldir):
    #
    # set up directories and files.
    #
    print "open input root file: %s" % iname
    self.inputRoot = ROOT.TFile( iname )
    #print self.inputRoot.GetEndpointUrl()
    self.htmldir = htmldir
    if not os.path.exists( self.htmldir ):
      print "html directory, %s, does not exist. Create a new one." % htmldir
      os.mkdir( self.htmldir )
  
    self.xmldir = self.htmldir
    self.refRoot = []
    try:
      self.mondir = os.path.join( os.getenv("LATMonRoot"), "TKR", "TkrMonitor" )
    except:
      self.mondir = os.path.join( os.getenv( "CALIBTKRUTILROOT" ), "python" )
    if not os.path.exists( self.mondir ):
      print "TKR Monitor directory, %s, does not exist." % self.mondir
      sys.exit()
    path = os.path.join( self.mondir, "refRoot", "*ref*.root" )
    refs = glob.glob( path )[:2]
    for ref in refs:
      print "open ref foot file: %s" % ref
      self.refRoot.append( ROOT.TFile( ref ) )

    #
    # initialize some parameters
    #
    self.values = {}
    for (key, type, vtype) in keys:
      self.values[key] = array.array( typecode[vtype], [-1]*nelem[type] )

    self.latave = {}
    for key in ("trigEff", "towerEff", "fracNHstrip"):
      self.latave[key] = -1 

    self.towerave = {}
    for key in ("layerEff" ,"layerOcc", "stripOcc", "fracSat", "totPeak", "layerdXY", "fracNHstrip"):
      self.towerave[key] = {}
      for tower in range(nTower):
        self.towerave[key][tower] = -1
    
    self.logs = {}
    for type in alertTypes:
      self.logs[type] = {}
      for level in alertLevels:
        self.logs[type][level] = []

    self.striplogs = {}   
    self.badStrips = {}
    for tower in range(nTower):
      self.striplogs[tower] = {}
      self.badStrips[tower] = {}
      for unp in range(nPlane):
        self.striplogs[tower][unp] = {}
        self.badStrips[tower][unp] = []
        for type in ("stripOcc","stripEff"):
          self.striplogs[tower][unp][type] = []
    
        
    self.towtemp = []   
    self.laytemp = []   
    self.hists = {}
    self.alerts = []
    self.towerHistsRef =[]
    self.towerHists = []
    self.layerHists = []
    self.savedPictures = []
    for tower in range(nTower):
      self.towerHists.append( {} )
      self.towerHistsRef.append( {} )
      for name in towerHistNames:
        if name[-1] == "T":
          (bin, nbin) = (timebin, timebin*2)   
        else : (bin, nbin) = (nPlane, nPlane)
        hname = "%s_T%d" % (name, tower)
        title = "%s in tower %d -%s" % (name, tower, timeStamp)
        self.towerHists[-1][name] = ROOT.TH1F(hname,title,nbin,0,bin)
        hname = "%s_T%dRef" % (name, tower)
        title = "%s in tower %d -%s (ref)" % (name, tower, timeStamp)
        self.towerHistsRef[-1][name] = ROOT.TH1F(hname,title,nbin,0,bin)
        
                        
      self.alerts.append( [] )
      self.layerHists.append( [] )
      self.savedPictures.append( [] )
      for unp in range(nPlane):
        self.alerts[-1].append( [] )
        self.layerHists[-1].append( {} )
        self.savedPictures[-1].append( {} )
        name = "stripOcc"
        hname = "%s_T%d_L%d" % (name, tower, unp)
        title = "%s in tower %d Layer %d -%s" % (name, tower, unp, timeStamp)
        self.layerHists[-1][-1][name] = ROOT.TH1F(hname,title,nStrip,0,nStrip)

        name = "stripEff"
        hname = "%s_T%d_L%d" % (name, tower, unp)
        title = "%s in tower %d Layer %d -%s" % (name, tower, unp, timeStamp)
        self.layerHists[-1][-1][name] = ROOT.TH1F(hname,title,nGTFE,0,nStrip)

        for name in ("layerOccT","stripOccT"): 
          (bin, nbin) = (timebin, timebin*2)
          hname = "%s_T%d_L%d" % (name, tower, unp)
          title = "%s in tower %d Layer %d -%s" % (name, tower, unp, timeStamp)
          self.layerHists[-1][-1][name] = ROOT.TH1F(hname,title,nbin,0,bin)

        
                                
    self.hists["layerEff"] = ROOT.TH2F("layerEff_map", \
                                       "Layer inefficiency map",\
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
    self.hists["trigEff"] = ROOT.TH1F("trigEff", "tower Trigger Efficiency", \
                                      nTower,0,nTower)
    self.hists["towerdX"] = ROOT.TH1F("towerdX", "tower delta X", \
                                      nTower,0,nTower)
    self.hists["towerdY"] = ROOT.TH1F("towerdY", "tower delta Y", \
                                      nTower,0,nTower)
    self.hists["totPeak"] = ROOT.TH1F("totPeakDist",\
                                      "TOT peak distributioin",25,2,7)
    self.hists["layerdXY"] = ROOT.TH2F("layerdXY_map", \
                                       "Layer displacement map",\
                                       nTower,0,nTower, nPlane,0,nPlane)
    self.hists["fracNHstrip"] = ROOT.TH1F("fracNHstrip", "Fraction of Dead strips", \
                                       nTower,0,nTower)
    self.hists["deadGTFE"] = ROOT.TH1F("deadGTFE", "number of dead GTFE", \
                                       nTower,0,nTower)

    self.readParamLimitsXml()
    

  #
  #*********************
  # read parameter limits from xml file
  #*********************
  #
  def readParamLimitsXml(self):
    self.limits={}
    xname = os.path.join( self.mondir, "paramLimitsDefault.xml"  )
    ### xname = os.path.join( self.mondir, "paramLimitsTest.xml"  )
    print "read xml: %s" % xname
    dom = xml.dom.minidom.parse( xname )
    topElm = dom.getElementsByTagName("ParamLimits")[0]
    limits = topElm.getElementsByTagName("limit")
    for limit in limits:
      type = str( limit.getAttribute("type") )
      level = str( limit.getAttribute("level") )
      value = float( limit.getAttribute("value") )
      print type, level, value
      if type not in alertTypes or level not in alertLevels:
        print "Invalid alert type: %s or level: %s." % (type,level)
        self.limits[type] = value
      else: self.limits[type] = value
      print "mondir=",self.mondir
    

    
  #
  #*********************
  # get time stamps and run IDs
  #*********************
  #
  def getTimeStamps(self):
    tree = self.inputRoot.FindObjectAny( "timeStamps" )
    
    self.startTime = array.array( 'd', [-1] )
    self.endTime = array.array( 'd', [1] )
    self.firstRunId = array.array( 'L', [0] )
    self.lastRunId = array.array( 'L', [0] )
    
    for i in range(tree.GetEntries()):
      tree.GetEntry( i )
      if self.startTime[0] < 0 or tree.startTime < self.startTime[0]:
        self.startTime[0] = tree.startTime
      if tree.endTime > self.endTime[0]: self.endTime[0] = tree.endTime
      if self.firstRunId[0] < 1 or tree.firstRunId < self.firstRunId[0]:
        self.firstRunId[0] = tree.firstRunId
      if tree.lastRunId > self.lastRunId[0]: self.lastRunId[0] = tree.lastRunId
      
    if tree.firstRunId == tree.lastRunId: self.srunID = "%d" % tree.firstRunId
    else: self.srunID = "%d-%d" % ( tree.firstRunId, tree.lastRunId )

    tree = self.inputRoot.FindObjectAny( "tkrNoiseTree" )
    startTime = -1
    for i in range(tree.GetEntries()):
      tree.GetEntry(0)
      if startTime<0 or tree.startTime<startTime:
        startTime = tree.startTime
    self.runID = "%d" % startTime

    #Temporary
    if abs(self.firstRunId[0] - int(self.runID)) > 10000:
      self.runID = "%d" % self.firstRunId[0]

    print self.startTime, self.endTime, self.runID, self.srunID

  #
  #*********************
  # perform analysis
  #*********************
  #
  def analyzeTKR(self):
    self.ffit = ROOT.defLangau( "langau", 0, 30 )
    self.ffit.SetParNames( "LWidth", "MP", "Area", "GSigma" )
    
    self.analyzeTrigEff()
    #
    # loop towers/layers ann fit TOT distributions
    #
    for tower in range(nTower):
      self.analyzeLayerEfficiency( tower )
      self.analyzeStripEfficiency( tower )
      self.analyzeTimedep(tower)
      for unp in range(nPlane):
        self.fitTOT( tower, unp )
        self.analyzeOccupancy( tower, unp )

    self.analyzeTotPeak()
    self.analyzeDisplacement()

    saveHotStrips = False
    for tower in range(nTower):
      for unp in range(nPlane):
        if len(self.badStrips[tower][unp]) > 0: saveHotStrips = True

    if saveHotStrips:
      year, month, day, hour, minute = tkrUtils.getDateTime( str(int(self.startTime[0])) )
      tss = "%d/%d 20%02d, %d:%02d" % (month, day, year, hour, minute)
      year, month, day, hour, minute = tkrUtils.getDateTime( str(int(self.endTime[0])) )
      tse = "%d/%d 20%02d, %d:%02d" % (month, day, year, hour, minute)
      times = ( tss, tse )
      tkrUtils.writeTFEXml( self.badStrips, self.runID, self.xmldir, times,\
                            self.runID, tagv)
      tkrUtils.writeHotXml( self.badStrips, self.runID, self.xmldir, times,\
                            self.runID, tagv)
      # tkrUtils.writeHotXml(self.badStrips,self.runID,"/afs/slac.stanford.edu/u/ki/nishino")
        

  
  #
  #***********************
  # analyze reference data
  #***********************
  #
  def analyzeTKRref(self, key, tower, id):
    hlhit = self.refRoot[id].FindObjectAny("lhitT"+str(tower))
    hltrk = self.refRoot[id].FindObjectAny("ltrkT"+str(tower))
    hresP = self.refRoot[id].FindObjectAny("resProfT"+str(tower))
    for unp in range(nPlane):
      ielm = unp + tower*nPlane
      lname = tkrUtils.g_layerNames[unp]

      ## layer efficiency 
      if key == "layerEff":  
        lhit = hlhit.GetBinContent(unp+1)
        ltrk = hltrk.GetBinContent(unp+1)
        if ltrk == 0.0: ltrk = 1.0
        eff = lhit / ltrk
        err = eff*(1-eff) / ltrk
        if err > 0.0: err = math.sqrt( err )
        else: err = 1.0 / ltrk
        self.towerHistsRef[tower]["layerEff"].SetBinContent(unp+1,eff*100)
        self.towerHistsRef[tower]["layerEff"].SetBinError(unp+1,err*100)

      ## layer Occupancy
      if key == "layerOcc" or  "stripOcc" or  "fracSat":  
        hmul = self.refRoot[id].FindObjectAny( "hTkrNoiseMulT%d%s" %(tower,lname) )
        hmap = self.refRoot[id].FindObjectAny( "hTkrHitMapT%d%s" %(tower,lname) )
        lexp = hmul.Integral()
        locc = hmul.Integral(2,128)
        lsat = hmul.Integral(nsat,128)
        socc = hmap.Integral()
        if lexp != None:
          if lexp == 0.0: lexp = 1.0
          locc /= lexp
          socc /= (lexp*nStrip)
          lsat /= lexp
        else:
          locc = 0.0
          socc = 0.0
          lsat = 0.0
        self.towerHistsRef[tower]["layerOcc"].SetBinContent(unp+1,locc)
        self.towerHistsRef[tower]["stripOcc"].SetBinContent(unp+1,socc)
        self.towerHistsRef[tower]["fracSat"].SetBinContent(unp+1,lsat)

      ## displacement
      if key == "layerdXY":
        resP = hresP.GetBinContent(unp+1)
        resP_err = hresP.GetBinError(unp+1)
        resPdef = resScale * ( self.values["layerdXY"][ielm] - resP )
        resPdef_err = self.values["layerdXY_err"][ielm]**2 + resP_err**2
        if resPdef_err > 0:
          resPdef_err = resScale * math.sqrt( resPdef_err )
        self.towerHistsRef[tower]["layerdXY"].SetBinContent(unp+1,resPdef)
        self.towerHistsRef[tower]["layerdXY"].SetBinError(unp+1,resPdef_err)
        
    if key == "layerdXY":
      ymax = self.towerHistsRef[tower]["layerdXY"].GetMaximum()
      self.towerHists[tower]["layerdXY"].SetMaximum( ymax*1.5 )
      self.towerHists[tower]["layerdXY"].SetMinimum( -ymax*1.5 )

    if id == 0: self.towerHistsRef[tower][key].SetLineColor(4)
    else: self.towerHistsRef[tower][key].SetLineColor(5)
    self.towerHistsRef[tower][key].SetLineWidth(2)
    return self.towerHistsRef[tower][key] 

      
  #
  #*********************
  # analyze trigger efficiencies
  #*********************
  #
  def analyzeTrigEff(self):
    htrk = self.inputRoot.FindObjectAny( "sixInARowMIP" )
    htrg = self.inputRoot.FindObjectAny( "sixInARowWithTrigMIP" )
    for tower in range(nTower):
      ntrk = htrk.GetBinContent( tower+1 )
      ntrg = htrg.GetBinContent( tower+1 )
      if ntrk > 0:
        teff = ntrg/ntrk
        terr = teff*(1-teff)/ntrk
      else:
        teff = 0.0
        terr = 1.0E-5
      if terr >0.0: terr = math.sqrt(terr)
      else: terr = 1/ntrk
      key = "trigEff"
      self.values[key][tower] = teff
      self.values[key+"_err"][tower] = terr
      self.hists[key].SetBinContent(tower+1,teff*100)
      self.hists[key].SetBinError(tower+1,terr*100)

      if teff < self.limits["trigEff"]:
        alert = "trigger efficiency, %.1f%s < %.1f%s" % (teff*100, "%", self.limits["trigEff"]*100, "%" )
        self.logAlerts( tower, 0, alert, "trigEff" ,alertLevels[1] )   # 0 is dummy

      
    if htrk.Integral() > 0:
      self.latave["trigEff"] = (100*htrg.Integral()/htrk.Integral(), "trigger efficiency")
    else:
      self.latave["trigEff"] = (100.0, "trigger efficiency")      

  #
  #*********************
  #set Histogram label
  #*********************
  #
  def setHistLabel(self, hist, key, id):    
    if id == "tower": lst = towerHistNames[key]
    if id == "layer": lst = layerHistNames[key]
    if id == "lat": lst = latHistNames[key]
    hist.SetXTitle(lst[1])
    if id == "layer": hist.SetYTitle(lst[2])
    else: hist.SetYTitle(lst[2]+lst[-1])
    hist.GetXaxis().CenterTitle(1)
    hist.GetYaxis().CenterTitle(1)
    hist.SetLineWidth(3)

  #
  #**************************
  #draw limits into histogram
  #**************************
  #
  def setHistLimits(self, key, hist):
    xmin = hist.GetBinLowEdge(1)
    xmax = hist.GetBinLowEdge( hist.GetNbinsX()+1 )
    if key[-3:] == "Eff":
      limit = self.limits[key] * 100
    elif key == "fracNHstrip":
      limit = self.limits[key] * 100
    else: limit = self.limits[key]
    tl = ROOT.TLine( xmin, limit, xmax, limit )
    tl.SetLineColor( 2 )
    tl.SetLineWidth( 3 )
    return tl
  
  #
  #********************
  # set hist Max,Min
  #********************
  #                need review
  def setHistMaxMin( self, hist, limit ):
    ymax = hist.GetMaximum()
    if ymax < limit: ymax = limit * 1.2
    hist.SetMaximum( ymax )
    ymin = hist.GetMinimum()
    if ymin > limit: ymin = limit * 0.8
    hist.SetMinimum( ymin )

  #
  #*********************
  # fit TOT distributions 
  #*********************
  #
  def fitTOT(self, tower, unp):
    ielm = unp + tower*nPlane
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
    fracBadTot = chargeHist.Integral(1,bin) / entries
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
    self.ffit.SetLineColor(2)
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
    if par[1] > 0: self.hists["totPeak"].Fill( par[1] )
    self.towerHists[tower]["totPeak"].SetBinContent(unp+1,par[1])
    self.towerHists[tower]["totPeak"].SetBinError(unp+1,error[1])
        
    
  #
  #*********************
  # analyze TOT peak distribution and outliers
  #*********************
  #
  def analyzeTotPeak(self):
    totPeakDist = self.hists["totPeak"]
    fgaus = ROOT.TF1("fgaus", "gaus", 3.0,5.5)
    fgaus.SetParameters(0,totPeakDist.Integral())
    fgaus.SetParameters(1,totPeakDist.GetMean())
    fgaus.SetParameters(2,totPeakDist.GetRMS())
    #fgaus.SetParLimits( 2, 0.10, 0.3 )

    fgaus.SetLineColor(2)
    fgaus.SetLineWidth(2)
    totPeakDist.Fit("fgaus")

    mean = fgaus.GetParameter(1)
    rms = fgaus.GetParameter(2)
    chi2 = fgaus.GetChisquare()

    limit = rms * self.limits["totPeakSigma"]

    self.totPeakAlerts = []
    sum = 0
    for ielm in range(nTower*nPlane):
      value = self.values["TOT_Peak"][ielm]
      sum += value
      if (ielm+1)%nPlane == 0:
        self.towerave["totPeak"][(ielm+1)/nPlane-1] = sum/nPlane
        sum = 0
      if value > 0 and abs(value-mean) > limit:
        tower, unp = divmod( ielm, nPlane )
        alert = "TOT peak outlier: %.1f" % self.values["TOT_Peak"][ielm]
        self.logAlerts( tower, unp, alert, "totPeak", alertLevels[1] )
    

  #
  #*********************
  # analyze layer effciencies
  #*********************
  #
  def analyzeLayerEfficiency(self, tower):
    hlhit = self.inputRoot.FindObjectAny("lhitT"+str(tower))
    hltrk = self.inputRoot.FindObjectAny("ltrkT"+str(tower))
    
    for unp in range(nPlane):
      ielm = unp + tower*nPlane
      lname = tkrUtils.g_layerNames[unp]
      
      #efficiency
      lhit = hlhit.GetBinContent(unp+1)
      ltrk = hltrk.GetBinContent(unp+1)
      if ltrk==0.0: ltrk = 1.0
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
      self.towerHists[tower][key].SetBinContent(unp+1,eff*100)
      self.towerHists[tower][key].SetBinError(unp+1,err*100)
      self.hists[key].SetBinContent(tower+1,unp+1,(1-eff)*100)
      self.towerHists[tower][key].SetMaximum( 100.5 )
      self.hists[key].SetMaximum( 100.5 )


    if tower == 0:
      self.sumttrk = 0.0
      self.sumthit = 0.0
    self.sumttrk += hltrk.Integral()
    self.sumthit += hlhit.Integral()
    if tower == nTower-1 :
      if self.sumttrk > 0:
      self.latave["towerEff"] = (100*self.sumthit/self.sumttrk, "hit efficiency")
    else:
      self.latave["towerEff"] = (100.0, "hit efficiency")      

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
    self.hists[key].SetBinContent(tower+1,teff*100)
    self.hists[key].SetBinError(tower+1,terr*100)
    if teff < self.limits["towerEff"]:
        alert = "tower efficiency, %.1f%s < %.1f%s" % (teff*100, "%", self.limits["towerEff"]*100, "%" )
        self.logAlerts( tower, 0, alert, "towerEff" ,alertLevels[1] )   # 0 is dummy

    self.towerave["layerEff"][tower] = teff*100

  #
  #*******************************
  # analyze strip effciencies, no-hit strips, Dead GTFE 
  #*******************************
  #
  def analyzeStripEfficiency(self, tower):
    nt,nl = 0.0,0.0
    for unp in range(nPlane):
      ielm = unp + tower*nPlane
      lname = tkrUtils.g_layerNames[unp]
      
      heocc = self.inputRoot.FindObjectAny("eoccT"+str(tower)+lname)
      htocc = self.inputRoot.FindObjectAny("toccT"+str(tower)+lname)
      eocc,tocc = 0.0,0.0
      nl = 0.0
      for strip in range(nStrip):
        eocc += heocc.GetBinContent(strip+1)
        tocc += htocc.GetBinContent(strip+1)
        if heocc.GetBinContent(strip+1) == 0: # count no-hit strip
          nt += 1
          nl += 1 
        if (strip+1)%nChannel == 0:
          if eocc == 0.0:
            hmap = self.inputRoot.FindObjectAny( "hTkrHitMapT%d%s" %(tower,lname) )
            #print tower, lname, hmap.Integral(strip-nChannel-1, strip)
            if hmap.Integral(strip-nChannel-1, strip) == 0:            
              self.hists["deadGTFE"].Fill(tower)
              alert = "Dead GTFE(No.%d)" % int((strip+1)/nChannel - 1)
              self.logAlerts( tower, unp, alert, "deadGTFE" ,alertLevels[0] )
          if tocc == 0.0: tocc = 1.0
          seff = eocc/tocc
          serr = seff*(1-seff)/tocc
          serr = math.sqrt(serr)
          if seff < self.limits["stripEff"]:
            alert = "GTFE(No.%d) averaged strip efficiency, %.1f%s < %.1f%s" \
                    % ((strip+1)/nChannel - 1, seff*100, "%", self.limits["stripEff"]*100, "%" )
            self.striplogs[tower][unp]["stripEff"].append((alert, (strip+1)/nChannel - 1, \
                                                           "layerEff", alertLevels[1]))
            #self.logAlerts( tower, unp, alert, "layerEff" ,alertLevels[2] )
          self.layerHists[tower][unp]["stripEff"].SetBinContent((strip+1)/nChannel,seff*100)
          self.layerHists[tower][unp]["stripEff"].SetBinError((strip+1)/nChannel,serr*100)          
          self.layerHists[tower][unp]["stripEff"].SetMaximum( 100.5 )
          eocc,tocc = 0.0,0.0
      self.towerHists[tower]["fracNHstrip"].SetBinContent(unp+1,100*nl/nStrip)
      if nl/nStrip > self.limits["fracNHstrip"]:
        alert = "fraction of no-hit strip, %.1f%s > %.1f%s" \
                % (100*nl/nStrip, "%", self.limits["fracNHstrip"]*100, "%" )
        self.logAlerts( tower, unp, alert, "fracNHstrip" ,alertLevels[1] )
          
    self.hists["fracNHstrip"].SetBinContent(tower+1,100*nt/nPlane/nStrip)

    sumdead = self.towerHists[tower]["fracNHstrip"].Integral()
    self.towerave["fracNHstrip"][tower] = sumdead/nPlane 

    if tower == nTower-1 :
      self.latave["fracNHstrip"] = (self.hists["fracNHstrip"].Integral()/nTower, "fraction of no-hit strips")




    
    
  #
  #*********************
  # analyze layer displacements
  #*********************
  #
  def analyzeDisplacement(self):
    htresPX = self.inputRoot.FindObjectAny("tresProfX") # Now useless
    htresPY = self.inputRoot.FindObjectAny("tresProfY") # Now useless
    for tower in range(nTower):
      hresP = self.inputRoot.FindObjectAny("resProfT"+str(tower))
      for unp in range(nPlane):
        ielm = unp + tower*nPlane
        lname = tkrUtils.g_layerNames[unp]
        # displacement        
        resP = hresP.GetBinContent(unp+1)
        resP_err = hresP.GetBinError(unp+1)        
        self.values["layerdXY"][ielm] = resP
        self.values["layerdXY_err"][ielm] = resP_err
        self.towerHists[tower]["layerdXY"].SetBinContent(unp+1,0)
        self.towerHists[tower]["layerdXY"].SetBinError(unp+1,0)
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

      hist = self.analyzeTKRref("layerdXY", tower, 0)
      sumresP = 0
      for unp in range(nPlane):
        ielm = unp + tower*nPlane
        lname = tkrUtils.g_layerNames[unp]
        # displacement
        resP = hist.GetBinContent(unp+1)
        error = hist.GetBinError(unp+1)
        sumresP += abs(resP)
        #self.hists["layerdXY"].SetBinContent(tower+1, unp+1, abs(resP))
        self.hists["layerdXY"].SetBinContent(tower+1, unp+1, resP)
        if abs(resP) > error * self.limits["layerdXYSigma"]:
          alert = "layer displacement, | %.1f/%.1f | > %.1f" \
                  % (resP, error, self.limits["layerdXYSigma"])
          self.logAlerts( tower, unp, alert, "layerdXY" ,alertLevels[2] )
      self.towerave["layerdXY"][tower] = sumresP/nPlane                        
        
        



  #
  #*********************
  # analyze occupancies
  #*********************
  #
  def analyzeOccupancy(self, tower, unp):
    ielm = unp + tower*nPlane
    if ielm%nPlane == 0:
      self.sumlexp = 0
      self.sumlocc = 0
      self.sumlsat = 0
      self.sumsocc = 0
    lname = tkrUtils.g_layerNames[unp]
    hmul = self.inputRoot.FindObjectAny( "hTkrNoiseMulT%d%s" %(tower,lname) )
    hmap = self.inputRoot.FindObjectAny( "hTkrHitMapT%d%s" %(tower,lname) )
    #hexp = self.inputRoot.FindObjectAny( "hTkrExposT%d%s" %(tower,lname) )
    lexp = hmul.Integral()
    locc = hmul.Integral(2,128)
    lsat = hmul.Integral(nsat,128)
    socc = hmap.Integral()

    self.sumlexp += lexp
    self.sumlocc += locc
    self.sumlsat += lsat
    self.sumsocc += socc
    if (ielm+1)%nPlane == 0:
      if self.sumlexp==0.0: self.sumlexp=1.0
      self.sumlocc /= self.sumlexp
      self.sumsocc /= (self.sumlexp*nStrip)
      self.sumlsat /= self.sumlexp
      self.towerave["layerOcc"][(ielm+1)/nPlane-1] = self.sumlocc
      self.towerave["stripOcc"][(ielm+1)/nPlane-1] = self.sumsocc
      self.towerave["fracSat"][(ielm+1)/nPlane-1] = self.sumlsat      
    if lexp != None:
      if lexp==0.0: lexp==1.0
      locc /= lexp
      socc /= (lexp*nStrip)
      lsat /= lexp
    else:
      locc,socc,lsat = 0,0,0
    self.maskHotstrip(tower, unp)
    if locc > self.limits["layerOcc"]:
      alert = "layer Occucpancy, %.1e > %.1e" \
              % (locc,self.limits["layerOcc"])
      self.logAlerts( tower, unp, alert, "layerOcc", alertLevels[0] )
    if socc > self.limits["stripOcc"]:
      alert = "layer averaged strip Occucpancy, %.1e > %.1e" \
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

    #strip Occupancy
    for strip in range(nStrip):
      shit = hmap.GetBinContent(strip+1)
      if lexp > 0.0:
        stocc = shit/lexp
      else:
        stocc = 0
      self.layerHists[tower][unp]["stripOcc"].SetBinContent(strip+1,stocc)
      if stocc > self.limits["stripOcc"]: 
        alert = "strip:%d occupancy, %.1e > %.1e" \
                              % (strip, stocc,self.limits["stripOcc"])
        self.logAlerts( tower, unp, alert, "stripOcc", alertLevels[1] )
        self.striplogs[tower][unp]["stripOcc"].append((alert, strip, "stripOcc", alertLevels[1]))


  #
  #*********************
  # time-dependence of Noise Occupancy
  #*********************
  #
  def analyzeTimedep(self, tower):
    httexp = self.inputRoot.FindObjectAny( "hTkrTowerSumExpTwr%d-R%s" % (tower,self.srunID) )
    httocc = self.inputRoot.FindObjectAny( "hTkrTowerOccTwr%d-R%s" % (tower,self.srunID) )


    for k in range(httexp.GetNbinsX()+1):
        ttexp = httexp.GetBinContent(k+1)
        ttocc = httocc.GetBinContent(k+1)
        if ttexp != 0:
          self.towerHists[tower]["stripOccT"].SetBinContent(k+1, ttocc/ttexp)
        else:
          self.towerHists[tower]["stripOccT"].SetBinContent(k+1, 0)
          
    for unp in range(nPlane):
      ielm = unp + tower*nPlane
      lname = tkrUtils.g_layerNames[unp]
      htexp = self.inputRoot.FindObjectAny( "hTkrExposT%d%s-R%s" %(tower,lname,self.srunID) )
      htlocc = self.inputRoot.FindObjectAny( "hTkrLayerOccT%d%s-R%s" %(tower,lname,self.srunID) )
      htsocc = self.inputRoot.FindObjectAny( "hTkrStripOccT%d%s-R%s" %(tower,lname,self.srunID) )
      for k in range(htexp.GetNbinsX()+1):
        texp = htexp.GetBinContent(k+1)
        tlocc = htlocc.GetBinContent(k+1)
        tsocc = htsocc.GetBinContent(k+1)
        if texp != 0:
          self.layerHists[tower][unp]["layerOccT"].SetBinContent(k+1, tlocc/texp)
          self.layerHists[tower][unp]["stripOccT"].SetBinContent(k+1, tsocc/texp)
        else:
          self.layerHists[tower][unp]["layerOccT"].SetBinContent(k+1, 0)
          self.layerHists[tower][unp]["stripOccT"].SetBinContent(k+1, 0)

  #
  #***********************
  # create Hot strip
  #***********************
  #
  def maskHotstrip(self, tower, unp):
    lname = tkrUtils.g_layerNames[unp]
    hWmap = self.inputRoot.FindObjectAny( "hTkrWHitMapT%d%s" %(tower,lname) )
    hmul = self.inputRoot.FindObjectAny( "hTkrNoiseMulT%d%s" %(tower,lname) )
    lexp = hmul.Integral()
    if lexp == 0: lexp = 1
    locc = hWmap.Integral()/lexp
    lsat = hmul.Integral(nsat,128)/lexp
    if lsat < self.limits["fracSat"]:
      while locc > self.limits["layerOcc"]:
        vmax,smax = -1,-1
        for strip in range(nStrip):
          socc = hWmap.GetBinContent(strip+1)
          if socc > vmax:
            if strip not in self.badStrips[tower][unp]:
              vmax = socc
              smax = strip
        self.badStrips[tower][unp].append(smax)
        print "mask strip: T%d L%d %d, %.1e" % (tower, unp, smax, vmax/lexp)
        locc = (locc*lexp-vmax)/lexp
      for strip in range(nStrip):
        socc = hWmap.GetBinContent(strip+1) / lexp
        if socc > self.limits["stripOcc"] \
               and strip not in self.badStrips[tower][unp]:
          self.badStrips[tower][unp].append(strip)
          print "mask strip: T%d L%d %d, %.1e" % (tower, unp, strip, socc)
                  

  #
  #*********************
  # create html report for the given layer
  #*********************
  #
  def htmlLayerReport(self, tower, unp, alert):
    ielm = unp + tower*nPlane
    lname = tkrUtils.g_layerNames[unp]
    tdir = "Tower%d" % tower
    ldir = "T%d%s" % (tower, lname)
    path = os.path.join( self.htmldir, tdir, ldir )
    if not os.path.exists( path ):
      print "html sub directory, %s, does not exist. Creat a new one." % path
      os.mkdir( path )
    if ielm not in self.laytemp:
      path = os.path.join(self.htmldir, tdir, ldir, "T%d%s.html" % (tower, lname))
      self.hout1 = open(path, "w")
      self.hout1.write( '<html>\n' )
      self.hout1.write( '<head><title> Tower%d-%s summary </title></head><BR>\n' %(tower, lname))
      self.hout1.write( '<body>\n' )
      self.htmlLine( self.hout1,'<B> Tower%d-%s summary </B>\n' % (tower, lname))
      layerLink = self.htmlLink( "Back to Top Page", "../../index.html")
      self.htmlLine(self.hout1, "%s" % layerLink)
      layerLink = self.htmlLink( "Back to Tower%d" %tower, "../tower%d.html" %tower)
      self.htmlLine(self.hout1, "%s" % layerLink)      
      self.hout1.write( '<table style="text-align: left; width: 100%;" ' \
                        + 'border="0" cellpadding="2" cellspacing="2">\n' )
      self.hout1.write( '<tbody>\n' )
                                          
    for key in self.layerHists[tower][unp].keys():
      pname = "%s_T%d%s-%s.png" % (key, tower, lname, self.runID )
      path = os.path.join( self.htmldir, tdir, ldir )
      path = os.path.join( path, pname )
      hist = self.layerHists[tower][unp][key]
      self.setHistLabel(hist, key, "layer")
      if key == "stripOcc" or key[-4:] == "OccT": ROOT.gPad.SetLogy( 1 )
      hist.Draw()
      if key == "stripOcc" or key == "stripEff":
          tl = self.setHistLimits(key, self.layerHists[tower][unp][key])
          tl.Draw("SAME")
      canvas.SaveAs( path )
      ROOT.gPad.SetLogy( 0 )
      self.savedPictures[tower][unp][key] = pname
      self.hout1.write( '<tr><td>' )
      self.htmlLine(self.hout1, layerHistNames[key][0])
      self.htmlImage( self.hout1, pname)
      self.hout1.write( '</td>\n<td>\n' )

      if key == "chargeHist":
        if "TOT" in alert:
          color = yellow
          self.htmlLine(self.hout1, "WARN",2,color)
        else:
          color = "NONE"
          self.htmlLine(self.hout1, "INFO",2,color)
        self.htmlLine(self.hout1, "TOT_LWidth: %.2f " % self.values["TOT_LWidth"][ielm])
        self.htmlLine(self.hout1, "TOT_Peak: %.1f" % self.values["TOT_Peak"][ielm],2,color)
        self.htmlLine(self.hout1, "TOT_GSigma: %.1f " % self.values["TOT_GSigma"][ielm])
        self.htmlLine(self.hout1, "TOT_FitProb: %.1e " % self.values["TOT_FitProb"][ielm])


      elif key == "stripOcc":
        if len(self.striplogs[tower][unp][key]) > 0:
          if len(self.striplogs[tower][unp][key]) > 10:
            self.htmlLine(self.hout1, "WARN: Noise flare in following strips",2,yellow)
            first = -100
            last = -100
            txt = ""
            for alert, strip, type, level in self.striplogs[tower][unp][key]:
              if first < 0:
                first = strip
                last = strip
              elif strip > last+1:
                if len(txt) > 0: txt += ", "
                if first == last: txt += "%d" % first
                else: txt += "%d-%d" % (first,last)
                first = strip
                last = strip
              elif strip == last+1: last = strip
            self.htmlLine(self.hout1, txt)
            self.htmlLine(self.hout1, "Number of affected strips: %d" % \
                          len(self.striplogs[tower][unp][key]) ) 
          else:
            self.htmlLine(self.hout1, "WARN",2,yellow)
            for alert, strip, type, level in self.striplogs[tower][unp][key]:
              if alert in "0.0%":
                self.htmlLine(self.hout1, alert, 2, red)
              else:
                self.htmlLine(self.hout1,"strip%d:%s" %(strip,alert))
            if len(self.badStrips[tower][unp]) > 0:
              for s in self.badStrips[tower][unp]:
                self.htmlLine(self.hout1,"Hot strip%d: Need to be Masked" %(s),2,red)
                
      else :
        if self.striplogs[tower][unp].has_key( key ):
          if len(self.striplogs[tower][unp][key]) > 0:
            self.htmlLine(self.hout1, "WARN",2,yellow)
            for alert, strip, type, level in self.striplogs[tower][unp][key]:
              if alert in "0.0%":
                self.htmlLine(self.hout1, alert, 2, red)
              else:
                self.htmlLine(self.hout1,"%s" %(alert))
                
      self.hout1.write( '</td></tr>\n' )
      
    if ielm not in self.laytemp:
      self.hout1.write( '</tbody>\n</table>\n' )
      self.hout1.write( '</body>\n</html>\n' )
                  
    self.laytemp.append(ielm)
    self.layerHists[tower][unp] = {} # reset to avoid duplicates
    return "%s/%s/T%d%s.html" % (tdir,ldir,tower,lname)
    
  #
  #*********************
  # create html report for the given tower
  #*********************
  #
  
  def htmlTowerReport(self, tower, unp):
    tdir = "Tower%d" % tower
    path = os.path.join( self.htmldir, tdir )
    if not os.path.exists( path ):
      print "html sub directory, %s, does not exist. Creat a new one." % path
      os.mkdir( path )     
    if tower not in self.towtemp:
      path = os.path.join(self.htmldir,tdir, "tower%d.html" % tower)
      self.hout1 = open(path, "w")
      self.hout1.write( '<html>\n' )
      self.hout1.write( '<head><title> Tower%d summary </title></head><BR>\n' % tower )
      self.hout1.write( '<body>\n' )
      self.htmlLine( self.hout1, '<B> Tower%d summary </B>\n' % tower)
      layerLink = self.htmlLink( "Back to Top Page", "../index.html")
      self.htmlLine( self.hout1, "%s" % layerLink)
      self.hout1.write( '<table style="text-align: left; width: 100%;" ' \
                        + 'border="0" cellpadding="2" cellspacing="2">\n' )
      self.hout1.write( '<tbody>\n' )
      
    for key in self.towerHists[tower].keys():
      pname = "%s_T%d-%s.png" % (key, tower, self.runID )
      path = os.path.join( self.htmldir, tdir )
      path = os.path.join( path, pname )     
      hist = self.towerHists[tower][key]
      self.setHistLabel(hist, key, "tower")
      hist.Draw()
      self.analyzeTKRref(key, tower, 0).Draw("SAME")
      #self.analyzeTKRref(key, tower, 1).Draw("SAME")
      if self.limits.has_key( key ):
        tl = self.setHistLimits(key, hist)
        tl.Draw("SAME")
      if key[-3:] == "Occ" or key[:5] == "fracS" or key[-4:] == "OccT": ROOT.gPad.SetLogy( 1 )  
      canvas.SaveAs( path )
      ROOT.gPad.SetLogy( 0 )
      self.hout1.write( '<tr><td>' )
      self.htmlLine(self.hout1, towerHistNames[key][0])
      self.htmlImage( self.hout1, pname )
      self.hout1.write( '</td>\n<td>\n' )
      
      label = towerHistNames[key][2]
      if key == "fracNHstrip" and self.latave["fracNHstrip"] > 50.0:
        self.htmlLine( self.hout1, "Poor statistics, No displayed", 2, "blue")
      if self.logs.has_key(key):        
        for level in alertLevels:
          if len(self.logs[key][level])>0:
            if level == alertLevels[0]: color = red
            elif level == alertLevels[1]: color = yellow
            else: color = None
            if self.logs[key][level][0][0] == tower:
              self.htmlLine( self.hout1, level, 2, color )           
            for t, u, a in self.logs[key][level]:
              if t == tower:
                lname = tkrUtils.g_layerNames[u]
                if key == "fracNHstrip":
                  if self.latave["fracNHstrip"] < 50.0:
                    layerLink = self.htmlLink( "Layer T%d%s" %(t,lname),\
                                               "T%d%s/T%d%s.html" %(t, lname, t, lname))
                    self.htmlLine( self.hout1, "T%d%s %s %s" %(t, lname, a, layerLink), 2, "None" )  
                elif key == "layerdXY":
                  self.htmlLine( self.hout1, "T%d%s %s" %(t, lname, a), 2, "None" )
                else:
                  layerLink = self.htmlLink( "Layer T%d%s" %(t,lname),\
                                             "T%d%s/T%d%s.html" %(t, lname, t, lname))
                  self.htmlLine( self.hout1, "T%d%s %s %s" %(t, lname, a, layerLink), 2, "None" )
      if self.towerave.has_key(key): 
         if key == "layerEff" or key == "fracNHstrip":
           txt = "Tower averaged %s: %.1f%s" % (label, self.towerave[key][tower], "%")
         elif key == "totPeak": txt = "Tower averaged %s: %.2f [fC]" % (label, self.towerave[key][tower]) 
         elif key == "layerdXY": txt = "Tower averaged %s: %.1f [micron]" \
              % (label, self.towerave[key][tower]) 
         else: txt = "Tower averaged %s: %.1e" % (label, self.towerave[key][tower]) 
         self.htmlLine( self.hout1, txt, 2, "green")
      self.hout1.write( '</td></tr>\n' )        
    if tower not in self.towtemp:
      self.hout1.write( '</tbody>\n</table>\n' )
      self.hout1.write( '</body>\n</html>\n' )
    self.towerHists[tower] = {} # reset to avoid duplicates
    self.towtemp.append(tower)  

    return "Tower%d/tower%d.html" % (tower,tower)

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
    for tower in range(nTower):
      for key in self.towerHists[tower].keys():
        self.towerHists[tower][key].Write()
    tf.cd()
    # layer based histograms
    for key in ["chargeHist", "stripOcc", "stripEff"]: 
      tf.mkdir( key )
      tf.cd( key )
      for tower in range(nTower):
        for unp in range(nPlane):
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
                   % (key, nTower, nPlane, vtype )
      elif type == "tower":
        leaflist = "%s[%d]/%s" % (key, nTower, vtype )
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
    if unp == -1:print "%s-%s: %s, T%d" % (level, type, alert, tower)
    else:print "%s-%s: %s, T%d %s" % (level, type, alert, tower, lname)
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
    lfile.write( "runID: %s\n" % self.runID )
    lfile.write( "intput file: %s\n" % self.inputRoot.GetName() )
    lfile.write( "reference files:\n" )
    for rfile in self.refRoot:
      lfile.write( rfile.GetName() + "\n" )
    gmt = tkrUtils.getGMT( self.startTime[0] )
    lfile.write( "Start Time: %s\n" \
                    % time.strftime( "%Y/%m/%d, %H:%M:%S", gmt ) )
    lfile.write( "Duration: %.0f seconds\n" \
                    % (self.endTime[0]-self.startTime[0]) )
    lfile.write( "Script version: %s\n" \
                    % __version__.split()[1] )
    lfile.write( "Script tag: %s\n" \
                    % __tag__.split()[1] )
    lfile.write( "****** alerts *****\n" )
    for level in alertLevels:
      for type in alertTypes:
        for tower, unp, alert in self.logs[type][level]:
          lname = tkrUtils.g_layerNames[unp]
          if "trigger efficiency" in alert or "tower efficiency" in alert:
            log = "%s %s: %s, T%d \n" \
                  % ("%"+level+"%", type, alert, tower)
          else:
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
    path = os.path.join( self.htmldir, "index.html" )
    self.hout = open( path,"w")

    self.hout.write( "<html>\n" )
    self.hout.write( "<head><title>LAT-TKR Monitor Report</title><head>\n" )
    self.hout.write( "<body>\n" )
    self.htmlLine( self.hout, "<B>LAT-TKR Monitor Report</B>", 3, "blue" )

    #
    # run info
    #
    self.hout.write( "<ul>" )
    self.htmlLine( self.hout, "<B>runID</B>: %s" % self.runID )
    self.htmlLine( self.hout, "<B>intput file</B>: %s" % self.inputRoot.GetName() )
    self.htmlLine( self.hout, "<B>reference files</B>:" )
    for rfile in self.refRoot:
      self.htmlLine( self.hout, rfile.GetName() )
    gmt = tkrUtils.getGMT( self.startTime[0] )
    self.htmlLine( self.hout, "<B>Start Time</B>: %s" \
                    % time.strftime( "%Y/%m/%d, %H:%M:%S", gmt ) )
    self.htmlLine( self.hout, "<B>Duration</B>: %.0f seconds" \
                    % (self.endTime[0]-self.startTime[0]) )
    self.htmlLine( self.hout, "<B>Script version</B>: %s" \
                    % __version__.split()[1] )
    self.htmlLine( self.hout, "<B>Script tag</B>: %s" \
                    % __tag__.split()[1] )
    self.hout.write( "</ul>" )

    #
    # summary plots and alerts
    #
    self.hout.write( '<table style="text-align: left; width: 100%;" ' \
                     + 'border="0" cellpadding="2" cellspacing="2">\n' )
    self.hout.write( '<tbody>\n' )
    for key in latHistNames:
      title = latHistNames[key][0]
      dopt = latHistNames[key][3]
      self.hout.write( '<tr><td>' )
      self.htmlLine( self.hout, title )
      hist = self.hists[key]
      self.setHistLabel(hist, key, "lat")
      hist.Draw( dopt )
      if dopt == " ":
        tl = self.setHistLimits(key, hist)
        tl.Draw("SAME")
      pname = "%s-%s.png" % (hist.GetName(), self.runID )
      path = os.path.join( self.htmldir, pname )
      canvas.SaveAs( path )
      self.htmlImage( self.hout, pname )
      self.hout.write( "</td>\n" )
      #if key != "towerEff":
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

  def htmlImage( self, hout, imageName, table=False ):
    txt = '<a href="%s"><img src="%s" width="300"></a>' \
          % (imageName, imageName)
    if table: txt = '<td>%s</td>' % txt
    hout.write( txt+"\n" )
                                                                                
  def htmlLine(self, hout, line, fsize=2, color=None):
    if color == None:
      hout.write( "<FONT SIZE=\"+%d\">%s</FONT><BR>\n" % (fsize, line) )
    else:
      hout.write( '<FONT SIZE="+%d" color="%s">%s</FONT><BR>\n' \
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
      if self.logs.has_key(key):
        if len(self.logs[key][level])>0:
          if level == alertLevels[0]: color = red
          elif level == alertLevels[1]: color = yellow
          else: color = None
          self.htmlLine( self.hout, level, 2, color)
          for tower, unp, alert in self.logs[key][level]:         
            if key == "fracNHstrip" and self.latave["fracNHstrip"] > 50.0:
              self.htmlLine( self.hout, "Poor statistics, No displayed", 2, "blue")
              break
            lname = tkrUtils.g_layerNames[unp]
            towerRef = self.htmlTowerReport( tower, unp )
            towerLink = self.htmlLink( "Tower %d"%tower, towerRef )
            layerRef = self.htmlLayerReport( tower, unp, alert )
            layerLink = self.htmlLink( "Layer %s"%lname, layerRef )
            if key == "towerEff" or key == "trigEff":
              self.htmlLine( self.hout, "%s: %s" % (alert, towerLink) )
            else:self.htmlLine( self.hout, "%s: %s %s" % (alert, towerLink, layerLink) )
    if self.latave.has_key(key): 
         txt = "Lat averaged %s: %.1f%s" % (self.latave[key][1], self.latave[key][0], "%") 
         self.htmlLine( self.hout, txt, 2, "green")
    self.hout.write( "</td>\n" )
  
  #
  #*********************
  # close files
  #*********************
  #
  def close(self):
    self.inputRoot.Close()
    for rfile in self.refRoot: rfile.Close()


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

  tkrMonitor = TkrMonitor( iname, htmldir )
  tkrMonitor.getTimeStamps()
  tkrMonitor.analyzeTKR()
  
  tkrMonitor.saveRoot( oname ) # save root before report.
  tkrMonitor.saveReports( aname, logname ) # this should be after saving ROOT

  tkrMonitor.close()
