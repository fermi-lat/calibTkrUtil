import sys, math, array, os, xml.dom.minidom, pickle, gzip
sys.path = [( os.path.join( os.environ.get( "ROOTSYS" ), "lib" ) )] + sys.path
import ROOT
import tkrUtils

ROOT.gROOT.SetBatch() # set ROOT into batch mode

nsatS = 15 # GTRC hit limit + 1
nsatL = 15

def calcEfficiency( nhit, entry ):
  if nhit > 0 and entry > 0:
    eff = nhit / entry
    err = eff*(1-eff) / entry
    if err >0.0: err = math.sqrt(err)
    else: err = 0.9/entry
    return eff, err
  else:
    return 1E-10, 0.9E-10


def getEfficiency( rfile ):
  #
  # trigger efficiency
  #
  hist = rfile.FindObjectAny( "sixInARowWithTrigMIP" )
  nhit = hist.Integral()
  hist = rfile.FindObjectAny( "sixInARowMIP" )
  entry = hist.Integral()
  teff, terr = calcEfficiency( nhit, entry )

  #
  # hit efficiency
  #
  nhit = 0
  entry = 0
  for tower in range(tkrUtils.g_nTowers):
    hist = rfile.FindObjectAny( "lhitT%d" % tower )
    nhit += hist.Integral()
    hist = rfile.FindObjectAny( "ltrkT%d" % tower )
    entry += hist.Integral()
  heff, herr = calcEfficiency( nhit, entry )

  return float(teff*100), float(terr*100), float(heff*100), float(herr*100)


def getOccupancy( rfile ):
  ( lexp, locc, socc, fnsat, fhsat ) = ( 0.0, 0.0, 0.0, 0.0, 0.0 )
  nexp = rfile.FindObjectAny( "nTrack" ).Integral() \
         * tkrUtils.g_nTowers * tkrUtils.g_nUniPlanes
  for tower in range(tkrUtils.g_nTowers):
    hist = rfile.FindObjectAny( "numHitGTRCT%d" % tower )
    fhsat += hist.GetBinContent( nsatL )
    if nsatS != nsatL: fhsat += hist.GetBinContent( nsatS )
    for unp in range(tkrUtils.g_nUniPlanes):
      lname = tkrUtils.g_layerNames[unp]
      hmul = rfile.FindObjectAny( "hTkrNoiseMulT%d%s" %(tower,lname) )
      hmap = rfile.FindObjectAny( "hTkrHitMapT%d%s" %(tower,lname) )
      lexp += hmul.Integral()
      locc += hmul.Integral(2,128)
      if unp<11: fnsat += hmul.Integral(nsatL,128)
      else: fnsat += hmul.Integral(nsatS,128)
      socc += hmap.Integral()
  locc, lerr = calcEfficiency( locc, lexp )
  socc, serr = calcEfficiency( socc, lexp*tkrUtils.g_nStrips)
  fnsat, fnerr = calcEfficiency( fnsat, lexp )
  fhsat, fherr = calcEfficiency( fhsat, nexp )
  return locc, lerr, socc, serr, fnsat, fnerr, fhsat, fherr
  

def analyzeOccupancy(rfile,startTime):
  errlog[startTime]=[]  
  for tower in range(tkrUtils.g_nTowers):
    for unp in range(tkrUtils.g_nUniPlanes):
      lname = tkrUtils.g_layerNames[unp]
      hmul = rfile.FindObjectAny( "hTkrNoiseMulT%d%s" %(tower,lname) )
      hmap = rfile.FindObjectAny( "hTkrHitMapT%d%s" %(tower,lname) )
      lexp = hmul.Integral()
      locc = hmul.Integral(2,128)
      if lexp != None:
        if lexp==0.0:lexp = 1.0
        locc /= lexp
      else:
        locc = 0.0
      if locc > limits["LayerOcc"]:
         errlog[startTime].append((tower, unp, locc))
  return errlog

def getRunInfo( rfile ):
   tree = rfile.FindObjectAny( "timeStamps" )

   startTime = -1
   endTime = -1
   firstRunId = -1
   lastRunId = -1

   for i in range(tree.GetEntries()):
     tree.GetEntry( i )
     #if startTime < 0 or tree.startTime < startTime:
     #  startTime = tree.startTime
     if tree.endTime > endTime: endTime = tree.endTime
     if firstRunId < 1 or tree.firstRunId < firstRunId:
       firstRunId = tree.firstRunId
     if tree.lastRunId > lastRunId: lastRunId = tree.lastRunId

   if firstRunId == lastRunId: runID = "%d" % firstRunId
   else: runID = "%d-%d" % ( firstRunId, lastRunId )

   tree = rfile.FindObjectAny( "tkrNoiseTree" )
   for i in range(tree.GetEntries()):
     tree.GetEntry( i )
     if startTime < 0 or tree.startTime < startTime:
       startTime = tree.startTime

   if abs(startTime-firstRunId) >100000:  # Temporary
     runID = int(startTime)               #

   return startTime, endTime, runID 


def htmlImage(hout, imageName, table=False ):
  txt = '<a href="%s"><img src="%s" width="180"></a>' \
        % (imageName, imageName)
  if table: txt = '<td>%s</td>' % txt
  hout.write( txt + "\n" )
  
def htmlLine(hout, line, fsize=2, color=None):
  if color == None:
    hout.write( "<FONT SIZE=\"+%d\">%s</FONT><BR>\n" % (fsize, line) )
  else:
    hout.write( '<FONT SIZE="+%d" color="%s">%s</FONT><BR>\n' \
                % (fsize, color, line) )
    
def htmlLink(word, href, color=None):
  txt = '<a href="%s">%s</a>' % (href, word)
  if color == None: return txt
  else: return '<font color="%s">%s</font>' % (color, txt)


def readParamLimitsXml():
  limits={}
  xname = os.path.join( "/afs/slac.stanford.edu/u/ki/nishino/slac/calibTkrUtil/calibTkrUtil-02-02-03/python/paramLimitsDefault.xml"  )
  print "read xml: %s" % xname
  dom = xml.dom.minidom.parse( xname )
  topElm = dom.getElementsByTagName("ParamLimits")[0]
  Limits = topElm.getElementsByTagName("limit")
  for limit in Limits:
    type = str( limit.getAttribute("type") )
    value = float( limit.getAttribute("value") )
    if type == "layerOcc":limits["LayerOcc"]=value
    if type == "towerEff":limits["HitEff"]=value*100
    if type == "trigEff":limits["TrigEff"]=value*100
    if type == "stripOcc":limits["StripOcc"]=value
  return limits    


#**************************************************************************
#
# main
#
#**************************************************************************
#**************************************************************************
if __name__ == '__main__':
  if len(sys.argv)>1:
    suf = sys.argv[1]
  else:
    print "please specify data set name."
    sys.exit()
    
  fname = "summaryRuns-%s.txt" % suf
  #rdir = os.path.join( os.getenv("LATMonRoot"), "TKR", "trend", suf )
  htmldir = "/afs/slac.stanford.edu/u/ki/nishino/trendMonitor-%s" %suf  
  #htmldir = "/afs/slac.stanford.edu/u/ki/nishino/public_html"  
  #if not os.path.exists( rdir ): os.mkdir( rdir )
  if not os.path.exists( htmldir ): os.mkdir( htmldir )
    
  ts = 236000000
  (timestamps, durations, runid, starttime) = ( [], [], [], [] )
  (teffs, terrs, heffs, herrs ) = ( [], [], [], [] )
  (loccs, lerrs, soccs, serrs ) = ( [], [], [], [] )
  (fnsats, fnerrs, fhsats, fherrs ) = ( [], [], [], [] )
  (totpars, totperrs) = ([],[])

  file = open( fname )
  runs = []
  for line in file:
    runs.append( line.split()[0] )
  
  #
  # read results from previous analyses
  #
  
  pname = "./runSummary-" + suf + ".pkl.gz"
  #pname = "pickles/runSummary-" + suf + ".pkl.gz"
  try:
    pfile = gzip.open( pname, "rb" )
    print "read pickle file:", pname
    summary = pickle.load( pfile )
    pfile.close()
  except:
    print "no pickle file or error in", pname
    summary = {}
    summary["files"]=[]

  update = False  
  limits = readParamLimitsXml()
  errlog = {}
  for rname in runs:
    if rname in summary["files"]:
      print "skip:",rname
      continue
    update = True   
    summary["files"].append( rname )
    rfile = ROOT.TFile.Open( rname )
    if rfile.IsZombie():
      print "rname is not proper ROOT path."
    else:
      print "open:", rfile.GetName()
      startTime, endTime, runID = getRunInfo( rfile )
      runinfo = getRunInfo( rfile )
      if abs(startTime-236224264)<100: # TkrBuf_taperTrcBuf
        nsatS = 12
        nsatL = 21
      elif abs(startTime-236223904)<100: # TkrBuf_alternatingSplit
        nsatS = 28
        nsatL = 28
      else:
        nsatS = 15
        nsatL = 15      
      errlog = analyzeOccupancy(rfile,startTime)
      eff = getEfficiency( rfile )
      #print "Trigger Efficiency: %.1f+-%.1f%s" % (teff,terr, "%")
      #print "Hit Efficiency: %.1f+-%.1f%s" % (heff,herr, "%")
      occ = getOccupancy( rfile )
      #print "occupancies: %.1e, %.1e" % (locc, socc)
      #print "saturation fractions: %.1e, %.1e" % (fnsat, fhsat)
      
      hist = rfile.FindObjectAny( "chargeAll" )
      func = hist.GetFunction( "langau2" )
      for i in range(6):
        totpars.append( func.GetParameter(i) )
        totperrs.append( func.GetParError(i) )
      if summary.has_key( startTime ):
        print "replace information for %d." % startTime
      
      summary[startTime] = ( rname, runinfo, eff, occ, totpars, totperrs, errlog)
    rfile.Close()

  #
  # save results
  #
  if update:
    pfile = gzip.open( pname, "wb" )
    pickle.dump( summary, pfile )
    pfile.close()
    print "updated pickle file:", pname
  
  keys = [("Trigger Efficiency (%)", "TrigEff", teffs, terrs),\
          ("Hit Efficiency (%)", "HitEff", heffs, herrs), \
          ("layer-OR occupancy", "LayerOcc", loccs, lerrs),\
          ("strip occupancy", "StripOcc", soccs, serrs),\
          ("TRC noise saturation fraction", "FracNoiseSat", fnsats, fnerrs),\
          ("TRC hit saturation fraction", "FracHitSat", fhsats, fherrs),\
          ]

  skey = summary.keys()
  skey.sort()
  for key in skey:
    if key == "files": continue
    (rname, runinfo, eff, occ, totpars, totperrs, errlog) = summary[key]
    startTime, endTime, runID = runinfo    
    timestamps.append( (startTime+endTime)*0.5-ts )
    durations.append( (endTime-startTime)*0.5 )    
    runid.append( runID )
    teff, terr, heff, herr = eff
    teffs.append( teff )
    terrs.append( terr )
    heffs.append( heff )
    herrs.append( herr )
    locc,lerr, socc,serr, fnsat,fnerr, fhsat,fherr = occ
    loccs.append( locc )
    lerrs.append( lerr )
    soccs.append( socc )
    serrs.append( serr )
    fnsats.append( fnsat )
    fnerrs.append( fnerr )
    fhsats.append( fhsat )
    fherrs.append( fherr )

  rfile = ROOT.TFile( "flightTrending.root", "RECREATE" )
  xp = array.array( "d", timestamps )
  ex = array.array( "d", durations )
  
  for label, name, vlist, elist in keys:
    yp = array.array( "d", vlist )
    if elist==None:
      gr = ROOT.TGraph( len(vlist), xp, yp )
    else:
      ey = array.array( "d", elist )
      gr = ROOT.TGraphErrors( len(heffs), xp, yp, ex, ey )
    gr.SetLineWidth(3)
    gr.SetLineColor(2)
    gr.SetTitle( label )
    gr.GetXaxis().SetTitle( "MET - 236E6" );
    gr.GetYaxis().SetTitle( label )
    if gr.GetYaxis().GetXmax() > 100: gr.SetMaximum( 100 )    
    gr.Draw( "AP" )
    if limits.has_key( name ):
      tl = ROOT.TLine( gr.GetXaxis().GetXmin(), limits[name], \
                       gr.GetXaxis().GetXmax(), limits[name] )
      tl.SetLineColor(4)
      tl.SetLineWidth(2)    
      tl.Draw("SAME")
    if gr.GetYaxis().GetXmax() < 1:
      if gr.GetYaxis().GetXmax()/gr.GetYaxis().GetXmin() < 10:
        gr.SetMaximum( gr.GetYaxis().GetXmax()*3 )
        gr.SetMinimum( gr.GetYaxis().GetXmin()/3 )
      ROOT.gPad.SetLogy( 1 )
    else: ROOT.gPad.SetLogy( 0 )
    ROOT.gPad.Update()
    ROOT.gPad.SaveAs( "%s/%sTrend-%s.gif" % (htmldir,name,suf) )
    ROOT.gPad.SaveAs( "%s/%sTrend-%s.pdf" % (htmldir,name,suf) )
    if gr.GetYaxis().GetXmax() > 92 and gr.GetYaxis().GetXmin() < 90:
      gr.SetMinimum( 90 )
      gr.Draw( "AP" )
      ROOT.gPad.SaveAs( "%s/%sMagTrend-%s.gif" % (htmldir,name,suf) )
      ROOT.gPad.SaveAs( "%s/%sMagTrend-%s.pdf" % (htmldir,name,suf) )

  ROOT.gDirectory.Write()
  print "Saved root file: %s" % rfile.GetName()
  rfile.Close()

      

  #
  # html Report
  #
  path = os.path.join( htmldir, "index.html" )
  hout = open( path, "w")

  hout.write( "<html>\n" )
  hout.write( "<head><title>LAT-TKR Performance Trending</title><head>\n" )
  hout.write( "<body><BR>\n" )

  htmlLine( hout, "LAT-TKR Performance Trending", 3, "blue")
  htmlLine( hout, "Input file:%s" %fname)
  htmlLine( hout, "run ID: %s - %s" %(runid[0], runid[len(summary["files"])-1]))
  hout.write( '<BR><table style="text-align: left; width: 100%;" ' \
              + 'border="1" cellpadding="2" cellspacing="2">\n' )
  hout.write( '<tbody>\n' )
  hout.write("<tr>")

  hout.write("<td></td>")  
  for label, name, vlist, elist in keys:
    hout.write("<td>")
    htmlLine(hout, label)
    hout.write("</td>")
  hout.write("</tr>")
  
  hout.write("<tr>")
  hout.write("<td></td>")
  for label, name, vlist, elist in keys:
    hout.write("<td>")
    htmlImage(hout, "%sTrend-%s.gif" % (name,suf))
    txt = htmlLink("download PDF", "%sTrend-%s.pdf" % (name,suf) )
    htmlLine(hout, txt, 1)
    hout.write("</td>")    
  hout.write("</tr>")

  hout.write("<tr><td>")
  htmlLine(hout, " runID ")
  hout.write("</td></tr>")

  val,err = {},{}
  #for key in summary.keys(): 
  skey = summary.keys()
  skey.sort()
  for key in skey:
    if key == "files":continue
    (rname, runinfo, eff, occ, totpars, totperrs, errlog) = summary[key]
    if totperrs[1] > 0.3*totpars[1]: continue

    startTime, endTime, runID = runinfo
    val["TrigEff"], err["TrigEff"], val["HitEff"], err["HitEff"] = eff
    val["LayerOcc"],err["LayerOcc"],val["StripOcc"],err["StripOcc"],\
        val["FracNoiseSat"],err["FracNoiseSat"],val["FracHitSat"],err["FracHitSat"] = occ
    hout.write("<tr>")
    hout.write("<td>")
    htmlLine(hout, runID, 1.8)
    hout.write("</td>")
    for label, name, vlist, elist in keys:
      hout.write("<td>")
      if limits.has_key( name ):
        if "Eff" in name:
          if val[name]<limits[name]:color = "red"
          else:color = "black"
          htmlLine(hout, "%.1f(+-%.1f)%s" % (val[name],err[name], "%"), 1.8, color)
        else:
          if val[name]>limits[name]:color = "red"
          else:color = "black"
          htmlLine(hout, "%.1e(+-%.1e)" % (val[name],err[name]), 1.8, color)      
      else:
        htmlLine(hout, "%.1e(+-%.1e)" % (val[name],err[name]), 1.8)   
      if name == "LayerOcc":
        if len(errlog[startTime]) > 0:
          for tower, unp, locc in errlog[startTime]:
            lname = tkrUtils.g_layerNames[unp]
            htmlLine(hout, "ERR T%d%s: %.1e" % (tower,lname,locc), 1.8, "red")
      hout.write("</td>")
    hout.write("</tr>")

  
  hout.write( "</body>\n</html>" )
  hout.close()
