import os, glob, xml.dom.minidom, math, time, datetime, sys, cStringIO
from xml.dom import xmlbuilder

import numarray, ROOT, pyfits


g_nTowers = 16
g_nTCC = 8
g_nTRC = 9
g_nFE = 24
g_nUniPlanes = 36
g_nChannels = 64
g_nStrips = g_nChannels * g_nFE
g_Serials = ["TkrFMA", "TkrFM2", "TkrFM14", "TkrFM15", \
             "TkrFMB", "TkrFM1", "TkrFM12", "TkrFM13", \
             "TkrFM5", "TkrFM3", "TkrFM7", "TkrFM9", \
             "TkrFM6", "TkrFM4", "TkrFM10", "TkrFM11"]
g_layerNames = [ "Y0", "X0", "X1", "Y1", "Y2", "X2", "X3", "Y3", "Y4", "X4", \
                 "X5", "Y5", "Y6", "X6", "X7", "Y7", "Y8", "X8", "X9", "Y9", \
                 "Y10", "X10", "X11", "Y11", "Y12", "X12", "X13", "Y13", \
                 "Y14", "X14", "X15", "Y15", "Y16", "X16", "X17", "Y17" ]
g_TrayNames = [ "0top", "1bot", "1top", "2bot", "2top", "3bot", "3top", \
               "4bot", "4top", "5bot", "5top", "6bot", "6top", \
               "7bot", "7top", "8bot", "8top", "9bot", "9top", \
               "10bot", "10top", "11bot", "11top", "12bot", "12top", \
               "13bot", "13top", "14bot", "14top", "15bot", "15top", \
               "16bot", "16top", "17bot", "17top", "18bot" ]


nfsDir = "/nfs/farm/g/glast/"
rawDataDirs  = ["/EM2/rawData/", "/Integration/rawData/"]
volumes = [ "u12", "u14", "u16", "u19", "u21", "u25", "u29", "u34", "u38", "u40", "u41", "u52" ]

def getLayerNameFromCCRC( tcc, trc ):
  view = ["Y", "X", "Y", "X"]
  layer = 2*trc + tcc/4
  return view[tcc/2]+str(layer)

def getCCRCFromLayerName( lname ):
  layer = int( lname[1:] )
  if lname[0] == 'X':
    if layer%2==0: tcc = (3,2)
    else: tcc = (6,7)
  else:
    if layer%2==0: tcc = (0,1)
    else: tcc = (5,4)
  trc = layer/2
  return (tcc, trc)

def getLayerNameFromSPT( spt ):
  if spt[1] == 'y': lname = "X"
  else: lname = "Y"
  rc = int( spt[2] )
  if spt[0] == "+": lname += str( 2*rc )
  else: lname += str( 2*rc+1 )
  return lname

def getSPTFromLayerName( lname ):
  if lname[0] == "X": side = "y"
  else: side = "x"
  layer = int( lname[1:] )
  rc = str( int( layer/2 ) )
  if layer%2 == 0:
    spt = "+" + side + rc
  else:
    spt = "-" + side + rc
  return spt

def getSPTFromISPT( ispt ): # this index is only used Eric's config root file.
  rc = ispt % g_nTRC
  cc = ispt / g_nTRC
  if cc == 1 or cc == 2: side = "y"
  else: side = "x"
  if cc < 2: spt = "+" + side + str(rc)
  else: spt = "-" + side + str(rc)
  return spt

def getMaskFromStripList( stripList, feId=0 ):
  vmask = 0L
  offset = int(feId)*g_nChannels
  for ch in range(g_nChannels):
    strip = ch + offset
    if strip in stripList: continue
    vmask = vmask|(1L<<ch)
  return "0x%x" % vmask  

  
def getStripListFromMask( mask, feId="0" ):
  if mask[:2] == "0x": vmask = int( mask[2:], 16 )
  else: vmask = int( mask )
  fe = int( feId )
  
  stripList = []

  for bit in range(g_nChannels):
    if vmask&(1L<<bit) == 0:
      stripList.append( bit + fe*g_nChannels )

  return stripList


def getGMT( timeStamp ):
  gmt0 = time.gmtime( 0 ) # GMT base time
  baseTimeDelta = datetime.datetime( 2001, 1, 1 ) \
                  - datetime.datetime( gmt0[0], gmt0[1], gmt0[2] )
  t0 = baseTimeDelta.days * 24 * 3600.0
  return time.gmtime(timeStamp+t0)


def getTimeStamps( xmlfile ):
  dom = xml.dom.minidom.parse( xmlfile )
  topElm = dom.getElementsByTagName("testReport")[0]
  tstuple = topElm.getAttribute("timestamp")
  tsarray = map(int, tstuple.replace(",", " ").strip("()").split( ))
  tsarray[0] -= 2000
  timestamp = int( "%2d%02d%02d%02d%02d" % tuple(tsarray[:5]) )
  timestamp = int( "%2d%02d%02d" % tuple(tsarray[:3]) )
  tstuple = dom.getElementsByTagName("StartTime")[0].childNodes[0].data
  tsarray = map(int, tstuple.replace(",", " ").strip("()").split( ))
  tsarray[0] -= 2000
  startTime = int( "%2d%02d%02d%02d%02d" % tuple(tsarray[:5]) )
  tstuple = dom.getElementsByTagName("EndTime")[0].childNodes[0].data
  tsarray = map(int, tstuple.replace(",", " ").strip("()").split( ))
  tsarray[0] -= 2000
  endTime = int( "%2d%02d%02d%02d%02d" % tuple(tsarray[:5]) )

  #print timestamp, startTime, endTime
  
  return (timestamp, startTime, endTime)

  
def findRunDirs( runs, keyword=None ):
  rundirlist = []
  if keyword == None:
    xname = "/*.xml"
    rname = "/*.root"
    fits = "/*.fits"
  else:
    xname = "/%s*.xml" % keyword
    rname = "/%s*.root" % keyword
    fits = "/%s*.fits" % keyword
    

  for run in runs:
    found = False
    for volume in volumes: # check variuous volumes
      for rawData in rawDataDirs:
        dir = nfsDir + volume + rawData
        runpath = dir + run
        if os.access(runpath, os.R_OK): # check readability of the pass
          # check if LICOS
          if os.access(runpath+"/LICOS", os.R_OK):
            xmls = glob.glob( runpath + "/LICOS/Analysis" + xname ) \
                   + glob.glob( runpath + "/LICOS/Analysis" + fits ) \
                   + glob.glob( runpath + "/LICOS/Analysis" + rname )
            #print runpath, xmls
            if len(xmls) > 0:
              runpath = glob.glob( runpath + "/LICOS/Analysis" )
            else:
              xmls = glob.glob( runpath + "/LICOS/analysis_TEM*" + xname )
              if len(xmls) > 0:
                runpath = glob.glob( runpath + "/LICOS/analysis_TEM*" )
            if len(xmls) > 0:
              # get time stamp from rcReport
              report = dir+run+"/LICOS/rcReport.out"
              if os.access(report, os.R_OK):
                # get time stamp from rcReport
                timestamp = getTimeStamps( report )
                rundirlist.append( [timestamp, run, xmls, runpath] )
                found = True
              else:
                print "rcReport does not exist in", report
            else: # no xml files
              print "xml files cannot be found under", runpath+"/LICOS"
          else: # pre-LICOS
            report = runpath+"/rcReport.out"
            if os.access(report, os.R_OK):
              # get time stamp from rcReport
              timestamps = getTimeStamps( report )
              xmls = glob.glob( runpath + xname )
              if len(xmls)==0:
                print "xml file does not exist in", runpath
              else:
                rundirlist.append( [timestamps, run, xmls, [runpath] ] )
                found = True
            else:
              print "rcReport does not exist in", runpath
          if found: print run + " is found in " + dir
        if found: break
      if found: break
    if not found: print "%s does not exist in volumes:" % run, volumes

  return rundirlist


def getBadStripsList( xmlfile, badStripsList=None, hotType=False ):

  if badStripsList == None:
    badStripsList = []
    for tw in range(g_nTowers):
      badStripsList.append( [] )
      for unp in range(g_nUniPlanes):
        badStripsList[-1].append( {} )
    
  # read dead strip xml file
  dom = xml.dom.minidom.parse( xmlfile )
  topElm = dom.getElementsByTagName("badStrips")[0]
  type = str( topElm.getAttribute("badType") )
  if type == "hot": hot = hotType
  else: hot = False
  towers = topElm.getElementsByTagName("tower")
  for tower in towers:
    #counts={}
    hwserial = str( tower.getAttribute("hwserial") )
    if hwserial not in g_Serials:
      print "Invalid hardware: ", hwserial
      continue
    row = int( str( tower.getAttribute("row") ) )
    col = int( str( tower.getAttribute("col") ) )
    tem = col + row*4
    if hwserial != g_Serials[tem]:
        print "wrong serial name %s for bay %d, %s" \
              % ( hwserial, tem, g_Serials[tem] )
    badStrips = badStripsList[tem]
    if hot:
      for unp in range(g_nUniPlanes):
        for nBad in range(3):
          if not badStrips[unp].has_key( nBad ):
            badStrips[unp][nBad]=[]
          badStrips[unp][nBad].append( [] ) # add new test results

    # loop planes
    planes = tower.getElementsByTagName("uniplane")
    for plane in planes:
      tray = str( plane.getAttribute("tray") )
      which = str( plane.getAttribute("which") )
      unp = g_TrayNames.index( tray+which )

      if hot:
        nOnbdTrig = str( plane.getAttribute("nOnbdTrig") )
        nOnbdData = str( plane.getAttribute("nOnbdData") )
        nBad = -1
        if nOnbdTrig == "true" and nOnbdData == "false": nBad = 0
        elif nOnbdTrig == "false" and nOnbdData == "true": nBad = 1
        elif nOnbdTrig == "true" and nOnbdData == "true": nBad = 2
        if nBad < 0: continue
      else:
        nBad = int( plane.getAttribute("howBad") )
        #if not counts.has_key( nBad ): counts[nBad] = 0
        if not badStrips[unp].has_key( nBad ):
          badStrips[unp][nBad]=[]
        badStrips[unp][nBad].append( [] ) # add new test results
      try:
        stripElm = plane.getElementsByTagName("stripList")[0]
        stripList = str( stripElm.getAttribute("strips") )
        strips = stripList.split( " " )
        #counts[nBad] += len(strips)
        for strip in strips:
          badStrips[unp][nBad][-1].append( int( strip ) )
      except: pass

  return badStripsList

def parseTFEFile( xmlFile, badStripsList=None ):
  
  if badStripsList == None:
    badStripsList = []
    for tw in range(g_nTowers):
      badStripsList.append( [] )
      for unp in range(g_nUniPlanes):
        badStripsList[-1].append( {} )
        for nBad in range(2):
          badStripsList[-1][-1][nBad] = []
          

  for tw in range(g_nTowers):
    for unp in range(g_nUniPlanes):
      for nBad in range(2):
        badStripsList[tw][unp][nBad].append( [] )

  # read TFE xml file
  dom = xml.dom.minidom.parse( xmlFile )
  topElm = dom.getElementsByTagName("LATC_XML")[0]
  towers = topElm.getElementsByTagName("TEM")
  for tower in towers:
    tem = int( tower.getAttribute("ID") )
    #print tem, g_Serials[tem]
    badStrips = badStripsList[tem]

    # loop SPT
    planes = tower.getElementsByTagName("SPT")
    for plane in planes:
      spt = str( plane.getAttribute("ID") )
      lname = getLayerNameFromSPT( spt )
      unp = g_layerNames.index( lname )
      fes = plane.getElementsByTagName("TFE")
      for fe in fes:
        feId = str( fe.getAttribute("ID") )
        masks = [ "trig_enable", "data_mask" ]
        for maskName in masks:
          maskElm = fe.getElementsByTagName(maskName)[0]
          stripList = getStripListFromMask( maskElm.childNodes[0].data, feId )
          nBad = masks.index( maskName )
          #print tem, unp, feId, nBad, stripList
          for strip in stripList:
            badStrips[unp][nBad][-1].append( int( strip ) )
      #print tem, g_layerNames[unp], badStrips[unp][0][-1], badStrips[unp][1][-1]

  return badStripsList


def convDac( dac ):
  # reverse the bit order
  revbits = 0
  for ibit in range(14):
    digit = 1<<ibit
    if digit & dac:
      revbits |= 1<<(13-ibit)
  # calibDAC is bit 0 to 5
  calibDac = revbits & 0x3F
  # calibRange is bit 6
  calibRange = ( revbits >> 6 ) & 0x1
  # threshold DAC is bit 7 to 12
  thrDac = (revbits >> 7 ) & 0x3F
  # threshold range is bit 13
  thrRange = (revbits >> 13 ) & 0x1
  return (thrRange, thrDac, calibRange, calibDac)


def encodeDac(thrRange, thrDac, calibRange, calibDac):

  revbits = (thrRange & 0x1) << 13
  revbits += (thrDac & 0x3F) << 7
  revbits += (calibRange & 0x1) << 6
  revbits += (calibDac & 0x3F)

  dac = 0
  for ibit in range(14):
    digit = 1<<ibit
    if digit & revbits:
      dac |= 1<<(13-ibit)

  return dac


def getDataThr( run ):
  runid = run[1]
  rundirs = run[3]
  datalist = []
  for rundir in rundirs:
    xmlfiles = glob.glob( rundir + "/TkrThrDispersion_data_thr_*" + runid + ".xml" )
    for fname in xmlfiles:
      # get serial ID and DAC values from xml file
      datalist += readThrXml( fname )
  return datalist


def readThrXml( xmlfile ):

  params = []
  print "open ", xmlfile

  # read charge scale xml file
  dom = xml.dom.minidom.parse(xmlfile)
  topElm = dom.getElementsByTagName("threshold")[0]
  towers = topElm.getElementsByTagName("tower")
  for tower in towers:
    hwserial = str( tower.getAttribute("hwserial") )
    if hwserial not in g_Serials:
      print "Invalid hardware: ", hwserial
      continue
    
    # loop planes
    planes = tower.getElementsByTagName("uniplane")
    if len(planes) > g_nUniPlanes:
      print "Too many uniplanes:", len(planes), xmlfile
      params.append( ( hwserial, None ) )
      continue
    
    # initialize arrays for tot parameters
    thr = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
    for plane in planes:
      tray = str( plane.getAttribute("tray") )
      which = str( plane.getAttribute("which") )
      unp = g_TrayNames.index( tray+which )
      strips = plane.getElementsByTagName("stripThresh")
      if len(strips) == 0:
        strips = plane.getElementsByTagName("stripDataThresh")        
      for strip in strips:
        id = int( strip.getAttribute("id") )
        thr[unp][id] = float( strip.getAttribute("data_thr") )
        
    params.append( ( hwserial, thr ) )

  return params


def getThrDac( run ):
  runid = run[1]
  rundirs = run[3]
  daclist = []
  for rundir in rundirs:
    # get serial ID and DAC values from xml file
    daclist += readThrDacFromXml( rundir + "/TkrThresholdCal_" + runid + ".xml" )
    continue
  
    # get DAC values from tnt file
    dactnt = readThrDacFromTnt( runid, rundir )
    # check tnt and xml DAC values are consistent
    diff = dacxml - dactnt
    dacmax = dacxml.max() # get truncation points for xml file
    dacmin = dacxml.min()
    for unp in range( g_nUniPlanes ):
      for fe in range( g_nFE ):
        if dactnt[unp][fe] != dacxml[unp][fe]: # xml DAC is truncated.
          if (dacxml[unp][fe]==dacmax and dactnt[unp][fe]>dacmax) \
             or (dacxml[unp][fe]==dacmin and dactnt[unp][fe]<dacmin):
            pass
          else:
            print "inconsistent thr DAC in xml:",  \
                  dacxml[unp][fe], dactnt[unp][fe], dacmin, dacmax
  return  daclist

def readThrDacFromTnt( runid, rundir ):
  prefix = "/TkrThresholdCal_FeThr_Layer"
  dacarray = numarray.zeros( (g_nUniPlanes,g_nFE), numarray.Float32 )
  for unp in range(g_nUniPlanes):
    lname = g_layerNames[unp]
    fname = rundir + prefix  + lname + "_" + runid + ".tnt"
    file = open( fname )
    lines = file.readlines()
    file.close()
    for line in lines:
      array = line.split( '\t' )
      try:
        fe = int( array[1] )
        temp = array[2].split( '\r' )
        dac = int( temp[0] )
        dacarray[unp][fe] = dac
      except: pass
        
  return dacarray


def readThrDacFromXml( xmlfile ):
  daclist = []

  # read dead strip xml file
  print "open ", xmlfile
  dom = xml.dom.minidom.parse( xmlfile )
  topElm = dom.getElementsByTagName("configuration")[0]
  latElm = topElm.getElementsByTagName("GLAT")[0]
  temlist = latElm.getElementsByTagName("GTEM")
  for tem in temlist:
    serialnos = tem.getElementsByTagName("serialnos")[0]
    hwserial = str( serialnos.getAttribute("tkr") )
    if hwserial not in g_Serials:
      print "Invalid hardware: ", hwserial
      continue
    gtcclist = tem.getElementsByTagName("GTCC")
    dacarray = numarray.zeros( (g_nUniPlanes,g_nFE), numarray.Float32 )

    # loop GTCC
    for gtcc in gtcclist:
      tcc = int( str( gtcc.getAttribute("ID") ) )
      gtrclist = gtcc.getElementsByTagName("GTRC")
      # loop GTRC
      for gtrc in gtrclist:
        trc = int( str( gtrc.getAttribute("ID") ) )
        lname = getLayerNameFromCCRC( tcc, trc )
        unp = g_layerNames.index( lname )
        gtfelist = gtrc.getElementsByTagName("GTFE")
        # loop GTFE
        for gtfe in gtfelist:
          tfe = int( str( gtfe.getAttribute("ID") ) )
          dac = str( gtfe.getElementsByTagName("dac")[0].firstChild.data )
          thrdac = ( convDac( int(dac,0) ) )[1]
          dacarray[unp][tfe] = thrdac
          #print hwserial, tcc, trc, unp, tfe, dac, thrdac
    daclist.append( (hwserial, dacarray) )

  return daclist


def modThrDacFromXml( infile, outfile, daclist ):

  dac = {}
  for (hwserial, dacarray) in daclist: dac[hwserial] = dacarray

  # read dead strip xml file
  print "open input file:", infile
  dom = xml.dom.minidom.parse( infile )
  topElm = dom.getElementsByTagName("configuration")[0]
  latElm = topElm.getElementsByTagName("GLAT")[0]
  temlist = latElm.getElementsByTagName("GTEM")
  for tem in temlist:
    serialnos = tem.getElementsByTagName("serialnos")[0]
    hwserial = str( serialnos.getAttribute("tkr") )
    if hwserial not in g_Serials:
      print "Invalid hardware: ", hwserial
      continue
    temId = g_Serials.index( hwserial )
    tem.setAttribute( "ID", str(temId) )
    gtcclist = tem.getElementsByTagName("GTCC")
    dacarray = dac[hwserial]

    # loop GTCC
    for gtcc in gtcclist:
      tcc = int( str( gtcc.getAttribute("ID") ) )
      gtrclist = gtcc.getElementsByTagName("GTRC")
      # loop GTRC
      for gtrc in gtrclist:
        trc = int( str( gtrc.getAttribute("ID") ) )
        lname = getLayerNameFromCCRC( tcc, trc )
        unp = g_layerNames.index( lname )
        gtfelist = gtrc.getElementsByTagName("GTFE")
        # loop GTFE
        for gtfe in gtfelist:
          tfe = int( str( gtfe.getAttribute("ID") ) )
          thrdac = hex( encodeDac( 0, int(dacarray[unp][tfe]), 0, 0 ) )
          gtfe.getElementsByTagName("dac")[0].firstChild.data = thrdac
          #print hwserial, tcc, trc, unp, tfe, dac, thrdac

  print "open output:",outfile
  outf = open( outfile, "w" )
  outf.write( dom.toxml( ) )
  outf.close()
  
  return


def getDateTime( time ):
  times = time.split( " " )
  if len(times) > 1: # convert from string
    month = times[0].split( "/" )[0]
    day = times[0].split( "/" )[1]
    year = times[1].split( "," )[0][-2:]
    hour = times[2].split(":")[0]
    minute = times[2].split(":")[1]
    return (year, month, day, hour, minute)
  else: # convert from MET
    MET = int( time )
    delta = datetime.timedelta( seconds=MET )
    date = datetime.datetime( 2001, 1, 1 ) + delta
    return (date.year-2000, date.month, date.day, date.hour, date.minute)
  

def readChargeScaleFromXml( xmlfile ):

  # read charge scale xml file
  dom = xml.dom.minidom.parse(xmlfile)
  topElm = dom.getElementsByTagName("chargeScale")[0]
  generic = topElm.getElementsByTagName("generic")[0]
  inputSample = generic.getElementsByTagName("inputSample")[0]
  startTime = str( inputSample.getAttribute("startTime") )
  times = startTime.split( " " )
  month = times[0].split( "/" )[0]
  day = times[0].split( "/" )[1]
  year = times[1].split( "," )[0][-2:]
  hour = times[2].split(":")[0]
  minute = times[2].split(":")[1]
  timeStamp = int( "%02d%02d%02d" % ( int(year), int(month), int(day) ) )
  startTime = int( "%02d%02d%02d%02d%02d" \
               % ( int(year), int(month), int(day), int(hour), int(minute) ) )
  stopTime = str( inputSample.getAttribute("stopTime") )
  times = stopTime.split( " " )
  month = times[0].split( "/" )[0]
  day = times[0].split( "/" )[1]
  year = times[1].split( "," )[0][-2:]
  hour = times[2].split(":")[0]
  minute = times[2].split(":")[1]
  stopTime = int( "%02d%02d%02d%02d%02d" \
               % ( int(year), int(month), int(day), int(hour), int(minute) ) )
  params = []
  towers = topElm.getElementsByTagName("tower")
  for tower in towers:
    hwserial = str( tower.getAttribute("hwserial") )
    if hwserial not in g_Serials:
      print "Invalid hardware: ", hwserial
      continue

    # loop planes
    planes = tower.getElementsByTagName("uniplane")
    if len(planes) > g_nUniPlanes:
      print "Too many uniplanes:", len(planes), xmlfile
      return [ ( hwserial, None ) ]
    
    # initialize arrays for tot parameters
    csarray = numarray.zeros( (g_nUniPlanes,g_nFE), numarray.Float32 )
    for plane in planes:
      tray = str( plane.getAttribute("tray") )
      which = str( plane.getAttribute("which") )
      unp = g_TrayNames.index( tray+which )
      fes = plane.getElementsByTagName("gtfeScale")
      for fe in fes:
        id = int( fe.getAttribute("id") )
        csarray[unp][id] = float( fe.getAttribute("chargeScale") )
        
    params.append( ( hwserial, csarray ) )

  #print startTime, stopTime
  return (timeStamp, startTime, stopTime, params)


def compareResults( rawValues, refValues, limits=None ):

  # create a mask for zero elements
  mask = (rawValues!=0) * (refValues!=0) * 1.0
  if limits != None:
    (minValue,maxValue) = limits
    mask *= (rawValues<maxValue) * (refValues<maxValue)
    mask *= (rawValues>minValue) * (refValues>minValue)
  rawValues *= mask
  refValues *= mask

  # calculate mean and rms of difference
  diff = rawValues - refValues
  mean = diff.mean()
  rms = diff.stddev()
  #print "mean %.1f, rms %.2f" % (mean,rms)
  maxdif = mean + rms * 4.0 + 1E-10
  mindif = mean - rms * 4.0 - 1E-10
  # check 3.5 sigma outliers 
  mask *= (diff<maxdif)
  mask *= (diff>mindif)
  rawValues *= mask
  refValues *= mask
  
  # calculate number of valid elements
  num = mask.sum()
  foutliers = 1 - float(num) / (rawValues.getshape()[0] * rawValues.getshape()[1])
  if num == 0:
    return (0.0, 0.0, -1.0, -1.0, -1.0, foutliers )
    
  nums = numarray.sum( numarray.transpose( mask ) )
  masknums = (nums!=0)
  nums += (1-masknums)*1E-10

  # calculate simple mean and rms
  # we cannot use mean() and stdev() method to exclude outliers
  mean = rawValues.sum() / num
  rmssq = (rawValues*rawValues).sum() / num - mean*mean
  if rmssq > 0.0: rms = math.sqrt( rmssq )
  else: rms = 0.0

  # take ratio
  refmin = (1-mask) * 1.0E-10
  ratios = mask*rawValues/(refValues+refmin)

  # calculate mean and rms of ratios
  meanR = ratios.sum() / num
  sq = ratios*ratios
  rmssq = sq.sum() / num - meanR*meanR
  if rmssq > 0.0: rmsR = math.sqrt( rmssq )
  else: rmsR = 0.0

  # calculate layer rms of ratios
  meanL = numarray.sum( numarray.transpose( ratios ) ) / nums
  rmssq = numarray.sum( numarray.transpose( sq ) )/nums - meanL*meanL
  masknums *= (rmssq>0)
  rmssq += (1-masknums)*1E-10 # avoid sqrt of zero
  rmsL = numarray.sqrt( rmssq )
  if masknums.sum() > 0:
    rmsLayerR = (rmsL*masknums).sum() / masknums.sum()
  else:
    rmsLayerR = 0.0
  
  return (mean, meanR, rms, rmsR, rmsLayerR, foutliers)


def readTotParams( fname ):
  if fname[-4:] == ".xml":
    params = readTotParamsFromXml( fname )
  elif fname[-5:] == ".root":
    params = readTotParamsFromRoot( fname )
  elif fname[-5:] == ".fits":
    params = readTotParamsFromFits( fname )
  else:
    print "Invalid file extention:", fname
    params = []

  return params
  

def readTotParamsFromFits( fname ):
  print "Reading", fname
  fits = pyfits.open( fname )
  pmap = {}
  next = len(fits)

  for iext in range(1,next):
    ext = fits[ iext ]
    #print ext.header["EXTNAME"], fname
    tower = int( ext.header["TOWERID"] )
    tray = ext.header["TRAY"]
    which = ext.header["WHICH"]
    
    hwserial = g_Serials[tower]
    if not pmap.has_key( hwserial ):
      # initialize arrays for tot parameters
      p0array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
      p1array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
      p2array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
      chisqarray = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
      pmap[hwserial] = [p0array, p1array, p2array, chisqarray]
    unp = g_TrayNames.index( tray+which )
    #print unp, ext.data
    idlist = ext.data.field( "STRIP" )
    p0list = ext.data.field( "OFFSET" )
    p1list = ext.data.field( "SLOPE" )
    p2list = ext.data.field( "QUAD_TERM" )
    chisqlist = ext.data.field( "CHI2" )
    for (id,p0,p1,p2,chisq) in zip(idlist,p0list,p1list,p2list,chisqlist):
      strip = int( id )
      pmap[hwserial][0][unp][strip] = float(p0)
      pmap[hwserial][1][unp][strip] = float(p1)
      pmap[hwserial][2][unp][strip] = float(p2)
      pmap[hwserial][3][unp][strip] = float(chisq)

  params = []
  for hwserial in pmap.keys():
    params.append( ( hwserial, pmap[hwserial] ) )

  return params


def readTotParamsFromXml( xmlfile ):

  params = []

  # read charge scale xml file
  dom = xml.dom.minidom.parse(xmlfile)
  topElm = dom.getElementsByTagName("tot")[0]
  towers = topElm.getElementsByTagName("tower")
  for tower in towers:    
    hwserial = str( tower.getAttribute("hwserial") )
    if hwserial not in g_Serials:
      print "Invalid hardware: ", hwserial
      continue

    # loop planes
    planes = tower.getElementsByTagName("uniplane")
    if len(planes) > g_nUniPlanes:
      print "Too many uniplanes:", len(planes), xmlfile
      return [ ( hwserial, None ) ]
    
    # initialize arrays for tot parameters
    p0array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
    p1array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
    p2array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
    chisqarray = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
    for plane in planes:
      tray = str( plane.getAttribute("tray") )
      which = str( plane.getAttribute("which") )
      unp = g_TrayNames.index( tray+which )
      strips = plane.getElementsByTagName("stripTot")
      for strip in strips:
        id = int( strip.getAttribute("id") )
        p0array[unp][id] = float( strip.getAttribute("intercept") )
        p1array[unp][id] = float( strip.getAttribute("slope") )
        p2array[unp][id] = float( strip.getAttribute("quad") )
        chisqarray[unp][id] = float( strip.getAttribute("chi2") )
        
    params.append( ( hwserial, [ p0array, p1array, p2array, chisqarray] ) )

  return params


def readTotParamsFromTnt( runid, rundir ):
  prefix = "/TkrTotGainNt_Layer"
  # initialize arrays for tot parameters
  p0array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
  p1array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
  p2array = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )
  chisqarray = numarray.zeros( (g_nUniPlanes,g_nStrips), numarray.Float32 )

  # loop through uniplanes
  for unp in range(g_nUniPlanes):
    lname = g_layerNames[unp]
    fname = rundir + prefix  + lname + "_" + runid + ".tnt"
    file = open( fname )
    lines = file.readlines()
    file.close()
    for line in lines:
      array = line.split( )
      try:
        strip = int( array[0] )
        p0array[unp][strip] = float( array[2] )
        p1array[unp][strip] = float( array[3] )
        p2array[unp][strip] = float( array[4] )
        chisqarray[unp][strip] = float( array[5] )
      except: pass
        
  return [ p0array, p1array, p2array, chisqarray ]


def getTotParams( run ):
  runid = run[1]
  rundirs = run[3]
  paramslist = []
  for rundir in rundirs:
    # get serial ID and tot parameter values from xml file
    xmlfile = rundir + "/TkrTotGain_" + runid + ".xml"
    if not os.access(xmlfile, os.R_OK): # check readability of the path
      xmlfile = rundir + "/TkrTotGain_" + str(int(runid)) + ".xml"
    params = readTotParamsFromXml( xmlfile )
    if params[0] == None: # bad xml, fall back to tnt.
      hwserial = params[1]
      totparams = readTotParamsFromTnt( runid, rundir )
      paramslist.append( ( hwserial, totparams ) )
      print runid, hwserial
    else:
      paramslist += params
  return  paramslist


def saveParams( params, paramNames ):
  cv = ROOT.TCanvas()
  varNames = ["mean", "mean-ratio-tower", "rms", "rms-ratio-tower", "mean-rms-ratio-layer", "fraction-outliers"]
  tslist = params.keys()
  tslist.sort()
  for pindex in range(len(paramNames)):
    logname = "Tkr%sTrend.txt" % paramNames[pindex]
    log = open( logname, "w" )
    for vindex in range(len(varNames)):
      line = "parameter name, %s, %s" \
             % ( paramNames[pindex], varNames[vindex] )
      print line
      log.write( line+"\n" )
      line = "serial ID, " + str(g_Serials).strip("[]").replace("'","").replace( "Tkr", "" )  
      print line
      log.write( line+"\n" )
      graphlist = [None]*g_nTowers
      gmax = -1E10
      gmin = 1E10
      tmax = 0.
      for timeStamp in tslist:
        line = "%06d" % timeStamp
        year = (timeStamp/10000) + 2000
        month = (timeStamp%10000)/100
        day = (timeStamp%100)
        deltaDays = datetime.datetime( year, month, day ) \
               - datetime.datetime( 2005, 1, 1 )
        days = deltaDays.days
        if days > tmax: tmax = days
        for towerId in range(g_nTowers):
          if params[timeStamp][towerId] is None:
            line += ", NA"
          else:
            param = params[timeStamp][towerId][pindex][vindex]
            if param == 0: line += ", 0"
            elif param == 1: line += ", 1"
            else:
              nfp = int( 3.0-math.log10( param ) )
              if nfp < 0: nfp = 0
              if nfp == 4: nfp = 3
              if nfp > 4: line += ", %.1e" % param
              else: line += ", %.*f" % (nfp, param)
            if graphlist[towerId] == None:
              graphlist[towerId] = ROOT.TGraph()
            ip = graphlist[towerId].GetN()
            graphlist[towerId].SetPoint( ip, days, param )
            if param > gmax: gmax = param
            if param < gmin: gmin = param
        print line
        log.write( line+"\n" )
      ROOT.gPad.Clear()
      gmax += abs( gmax * 0.1 )
      gmin -= abs( gmin * 0.1 )
      hr = cv.DrawFrame( 0, gmin, tmax, gmax)
      hr.SetXTitle( "time from 2005/Jan/1 (days)" )
      hr.SetYTitle( "%s (%s)" % (paramNames[pindex], varNames[vindex]) )
      for towerId in range(g_nTowers):
        if towerId==g_nTowers/2:
          ROOT.gPad.Clear()
          hr = cv.DrawFrame( 0, gmin, tmax, gmax)
          hr.SetXTitle( "time from 2005/Jan/1 (days)" )
          hr.SetYTitle( "%s (%s)" % (paramNames[pindex], varNames[vindex]) )
        hist = graphlist[towerId]
        #hist.SetMaximum( gmax )
        #hist.SetMinimum( gmin )
        hist.SetLineWidth( 2 )
        hist.SetLineColor( towerId%8+1 )
        hist.Draw( "LP" )
        if towerId%(g_nTowers/2) == g_nTowers/2-1:
          tnum = "T%d-%d" % (towerId-(g_nTowers/2)+1,towerId)
          ROOT.gPad.SaveAs( "pdf/%s-%s-%s.pdf" \
                            %(paramNames[pindex], varNames[vindex], tnum) )
      ROOT.gPad.Update()

    log.close()
    print "Log saved to", logname


def saveParamsToXml( params, paramNames, times, topDir=None, calibType=None ):

  timeStamp = time.strftime("%y%m%d-%H%M%S", time.gmtime() )
  if calibType == None: calibType = paramNames[0]
  if topDir == None:
    topDir = os.path.join( os.environ.get('TKRDATA'), 'trend', calibType )
    if not os.path.exists(topDir): os.mkdir( topDir )
    topDir = os.path.join( topDir, timeStamp )
    if not os.path.exists(topDir): os.mkdir( topDir )

  varNames = ["mean", "mean-ratio-tower", "rms", "rms-ratio-tower", \
              "mean-rms-ratio-layer", "fraction-outliers"]
  tslist = params.keys()
  tslist.sort()
  for tStamp in tslist:
    fileName = "%s-%06d-%s.xml" % (calibType, tStamp, timeStamp)
    file = os.path.join( topDir, fileName )
    print "saving", file
    handler = XmlHandler(file, "trend")
    doc_elem = handler.getDoc().documentElement
    
    doc_elem.setAttribute( "CalibStatus", "OK" )
    doc_elem.setAttribute( "CalibProcLebel", "DEV" )
    doc_elem.setAttribute( "creationDate", timeStamp )
    doc_elem.setAttribute( "Creator", os.getenv('USER') )
    doc_elem.setAttribute( "Software", sys.argv[0] )
    doc_elem.setAttribute( "Description", "" )
    doc_elem.setAttribute( "Version", timeStamp )
    doc_elem.setAttribute( "timestamp", "%06d"%tStamp )
    doc_elem.setAttribute( "startTime", "%010d"%times[tStamp][0] )
    doc_elem.setAttribute( "stopTime", "%010d"%times[tStamp][1] )
    doc_elem.setAttribute( "calibType", calibType )
    
    for towerId in range(g_nTowers):
      tower_elem = handler.goToTower( (towerId, "Tower") )
      if params[tStamp][towerId] != None:
        for pindex in range(len(paramNames)):
          paramName = paramNames[pindex]
          #vparams = ["%06d"%tStamp]
          #vparams = []
          #for param in params[tStamp][towerId][pindex]:
          #  vparams.append( param )
          #handler.addElement( tower_elem, paramName, ["timeStamp"]+varNames, \
          #                    vparams )
          handler.addElement( tower_elem, paramName, varNames, \
                              params[tStamp][towerId][pindex] )

    handler.saveAndClose()
    print "XML saved to", file

  
class XmlHandler:
  """
  """

  def __init__(self, fileName, type):
    if type is 'tot' or type is 'threshold':
      dtd = 'tkrCalibration.dtd'
    else:
      dtd = 'badStrips.dtd'
    try:
      dtdfile = os.path.join(os.getenv('CALIBTKRUTILROOT'), 'src', dtd)
      if not os.path.exists(dtdfile):
        dtdfile = os.path.join(os.getenv("LATMonRoot"),"TKR","TkrMonitor",dtd)
    except:
      dtdfile = os.path.join(os.getenv("LATMonRoot"), "TKR", "TkrMonitor", dtd)
    if not os.path.exists(dtdfile):
      print 'FATAL: Did not find the DTD file ' + dtdfile
      sys.exit()
    
    print "DTD:", dtdfile
    self.fileName = fileName
    
    #check if it is a new doc creation that is requested
    try:
      self.createNewDOM(dtdfile,type)
    except:
      print 'DOM creation failure'
      raise

  def createNewDOM(self,dtdfile,type):
    f = file(dtdfile,'r')
    string = f.read()
    if type is 'tot':
      template = ( \
          '<!DOCTYPE tot [\n' + string.replace('\n', '\r\n') +''']> \
          <tot > \
          </tot> \
      ''')
    elif type is 'threshold':
      template = ( \
          '<!DOCTYPE threshold [\n' + string.replace('\n', '\r\n') +''']> \
          <threshold > \
          </threshold> \
      ''')
    elif type == "TFE":
      template = ( \
          '<?xml version="1.0" ?> \
          <LATC_XML > \
          </LATC_XML> \
      ''')
    elif type == "trend":
      template = ( \
          '<TKR_trend > \
          </TKR_trend> \
      ''')
    else:
      template = ( \
          '<!DOCTYPE badStrips [\n' + string.replace('\n', '\r\n') +''']> \
          <badStrips > \
          </badStrips> \
      ''')

    builder = xmlbuilder.DOMBuilder()
    source = xmlbuilder.DOMInputSource()

    source.byteStream = cStringIO.StringIO(template)
    self.doc = builder.parse(source)

    doc_elem = self.doc.documentElement
    if type is 'dead' or type is 'hot':
      doc_elem.setAttribute('badType',type)


  def loadExistingFile(self, fileName, the_type):
    #I need to validate it
    self.doc = xml.dom.minidom.parse(open(fileName,'rw'))
    checktype = the_type
    if (the_type is 'dead') or (the_type is 'hot'):
      checktype = str(self.doc.documentElement.getAttribute('badType'))

    if checktype != the_type:
      print fileName + ' exists but is inconsistent with given type: ' + type
      raise

  def createGenericElement(self, instrument, timestamp, runid, calType, creator, version):
    generic = self.doc.createElement('generic')
    if instrument != "LAT" and instrument != "TWR": instrument = "TWR"
    generic.setAttribute('instrument', instrument)
    generic.setAttribute('timestamp',  timestamp)
    generic.setAttribute('runId',  runid)
    generic.setAttribute('calType',    calType)
    generic.setAttribute('creatorName',  creator)
    generic.setAttribute('creatorVersion',  version)
    generic.setAttribute('fmtVersion', "NA")
    return generic

  def createInSampleElement(self, start, stop, trig, mode, src):
    iS = self.doc.createElement('inputSample')
    iS.setAttribute('startTime', start)
    iS.setAttribute('stopTime',  stop)
    iS.setAttribute('triggers',  trig)
    iS.setAttribute('mode',      mode)
    iS.setAttribute('source',    src)
    new_text = self.doc.createTextNode('\n Output from LSR bad strips alg, running on MC data')
    iS.appendChild(new_text)
    return iS

  def addElement(self, element, name, attNames, attValues, content=None ):
    """
    add dom elements with some attributes and content
    """
    newElem = self.doc.createElement(name)
    for (attName,attValue) in zip(attNames,attValues):
      newElem.setAttribute(attName,str(attValue))
      self.doc.documentElement.appendChild(newElem)
      
    if content != None:
      textNode = self.doc.createTextNode( str(content) )
      newElem.appendChild( textNode )
    
    element.appendChild(newElem)
    return newElem


  def goToTower(self, tem, nOnbdCalib='', nOnbdTrig='', nOnbdData='' ):
    """
    Check current DOM for given tower element, and create it if
    not found
    """
    (tower, twserial) = tem
    if twserial == "TEM": # TEM mode
      towerElem = self.doc.createElement('TEM')      
      towerElem.setAttribute('ID',str(tower))
      self.doc.documentElement.appendChild(towerElem)
      return towerElem
    elif twserial == "Tower": # Tower mode
      towerElem = self.doc.createElement('Tower')      
      towerElem.setAttribute('ID',str(tower))
      self.doc.documentElement.appendChild(towerElem)
      towerElem.setAttribute('hwserial',g_Serials[tower])
      self.doc.documentElement.appendChild(towerElem)
      return towerElem
    
    row = str(tower/4)
    col = str(tower%4)
    for subnode in self.doc.childNodes:
      for sub2 in subnode.childNodes:
        if sub2.nodeName == 'tower':
          if (str(sub2.getAttribute('row'))==row) & (str(sub2.getAttribute('col'))==col):
            return sub2

    towerElem = self.doc.createElement('tower')
    towerElem.setAttribute('row',row)
    towerElem.setAttribute('col',col)
    if nOnbdCalib is not '': towerElem.setAttribute( "nOnbdCalib", nOnbdCalib )
    if nOnbdTrig is not '': towerElem.setAttribute( "nOnbdTrig", nOnbdTrig )
    if nOnbdData is not '': towerElem.setAttribute( "nOnbdData", nOnbdData )
    
    towerElem.setAttribute('hwserial',twserial)
    self.doc.documentElement.appendChild(towerElem)
    return towerElem


  def goToSPT(self, towerElem, unp, trigList, dataList ):
    """
    Check current DOM for given uniplane element, and create it if
    not found.
    Add bad strip list
    """

    # find if this SPT already exist
    lname = g_layerNames[unp]
    spt = getSPTFromLayerName( lname )
    sptElm = None
    for node in towerElem.childNodes:
      if node.nodeName == 'SPT':
        if node.getAttribute('ID') == str(spt):
          sptElm = node
          break

    # create a new SPT element if it does not exist.
    if sptElm == None:
      (tcc,trc) = getCCRCFromLayerName( lname )
      comment = 'Layer %s, unp %d, CC pair (%d,%d), RC %d' \
                % (lname,unp,tcc[0],tcc[1],trc)
      header = self.doc.createComment(comment)
      towerElem.appendChild(header)
      sptElm = self.doc.createElement( 'SPT' )
      sptElm.setAttribute( 'ID', str(spt) )
      towerElem.appendChild(sptElm)

    # loop over strips
    for type,stripList in zip(["trig_enable","data_mask"],[trigList,dataList]):
      chList = []
      for fe in range(g_nFE): chList.append( [] )
      for strip in stripList:
        fe = int( strip / g_nChannels )
        ch = int( strip % g_nChannels )
        chList[fe].append( ch )
      for fe in range(g_nFE):
        if len(chList[fe]):
          self.goToTFE( sptElm, fe, chList[fe], type )
    #print sptElm.toxml()
    return sptElm


  def goToTFE(self, sptElm, fe, chList, maskType ):
    """
    Check current DOM for given TFE element, and create it if
    not found.
    Add mask.
    """
    # find if this TFE already exist
    tfeElm = None
    for node in sptElm.childNodes:
      if node.nodeName == 'TFE':
        if node.getAttribute('ID') == str(fe):
          tfeElm = node
          break

    # create a new TFE element if it does not exist.
    if tfeElm == None:
      tfeElm = self.doc.createElement( 'TFE' )
      tfeElm.setAttribute( 'ID', str(fe) )
      textNode = self.doc.createTextNode( "0xffffffffffffffff" )
      tmask = self.doc.createElement( "trig_enable" )
      tmask.appendChild( textNode )
      textNode = self.doc.createTextNode( "0xffffffffffffffff" )
      dmask = self.doc.createElement( "data_mask" )
      dmask.appendChild( textNode )
      textNode = self.doc.createTextNode( "0x0" )
      cmask = self.doc.createElement( "calib_mask" )
      cmask.appendChild( textNode )
      tfeElm.appendChild(tmask)
      tfeElm.appendChild(dmask)
      tfeElm.appendChild(cmask)
      sptElm.appendChild(tfeElm)
    
    for node in tfeElm.childNodes:
      if node.nodeName == maskType:
        mask = node.childNodes[0].data
        #stripList = getStripListFromMask( maskElm.childNodes[0].data, feId )
        stList = getStripListFromMask( mask )
        for ch in chList:
          if ch not in stList: stList.append( ch )
        textNode = self.doc.createTextNode( getMaskFromStripList( stList ) )
        old = node.replaceChild( textNode, node.childNodes[0] )

    
  def goToUniplane(self, towerElem, unp, howBad='', \
                   nOnbdCalib='false', nOnbdTrig='false', nOnbdData='false' ):
    """
    Check current DOM for given uniplane element, and create it if
    not found
    """
    tray      = g_TrayNames[unp][:-3]
    whichside = g_TrayNames[unp][-3:]
    for node in towerElem.childNodes:
      if node.nodeName == 'uniplane':
        if (node.getAttribute('tray')==str(tray)) \
           and (node.getAttribute('which')==whichside) \
           and (node.getAttribute('howBad')==howBad):
          return node
        
    lname = g_layerNames[unp]
    (tcc,trc) = getCCRCFromLayerName( lname )
    comment = 'Layer %s, unp %d, CC pair (%d,%d), RC %d' \
              % (lname,unp,tcc[0],tcc[1],trc)
    header = self.doc.createComment(comment)
    towerElem.appendChild(header)
    uniplane = self.doc.createElement( 'uniplane' )
    uniplane.setAttribute( 'tray', str(tray) )
    uniplane.setAttribute( 'which', whichside )
    if howBad is not '':
      uniplane.setAttribute( 'howBad', howBad )
      uniplane.setAttribute( "nOnbdCalib", nOnbdCalib )
      uniplane.setAttribute( "nOnbdTrig", nOnbdTrig )
      uniplane.setAttribute( "nOnbdData", nOnbdData )
    towerElem.appendChild(uniplane)
    return uniplane


  def addStripList(self, uniplane, stripList):
    if stripList == []:
      return
    for node in uniplane.childNodes:
      if node.nodeName == 'stripList':
        existingList = node.getAttribute('strips')
        existingList = str(existingList).split(' ')
        stringList = []
        for l in stripList:
          stringList.append(str(l))
        allList = existingList+ stringList
        result = []
        for a in allList:
          if a not in result:
            result.append(a)
        node.setAttribute('strips',reformatStripList(result))
        return
    stripsElem = self.doc.createElement('stripList')
    if stripList.__class__() == []:
      stripList = reformatStripList(stripList)
    stripsElem.setAttribute('strips',stripList)
    uniplane.appendChild(stripsElem)

  def addStripSpan(self, uniplane, start, stop):
    """
    """
    for node in uniplane.childNodes:
      if node.nodeName == 'stripSpan':
        oldstart = atoi(str(node.getAttribute('first')))
        oldstop  = atoi(str(node.getAttribute('last')))
        newstart = oldstart
        newstop = stop
        if start<oldstart: newstart = start
        if stop>oldstop: newstop = stop
        node.setAttribute('first',str(newstart))
        node.setAttribute('last',str(newstop))
        return

    strips = self.doc.createElement('stripSpan')
    strips.setAttribute('first',str(start))
    strips.setAttribute('last',str(stop))
    uniplane.appendChild(strips)


  def addAttributes(self,motherElem,childName,attribute_dict):
    """
    """
    child_elem = self.doc.createElement(childName)
    for key in attribute_dict.keys():
      child_elem.setAttribute(key,str(attribute_dict[key]))
    motherElem.appendChild(child_elem)
    
##   def addStripAttributeDict(self,uniplane, dict):
##     """
##     """
##     strip_elem = self.doc.createElement('strip')
##     for key in dict.keys():
##       strip_elem.setAttribute(key,str(dict[key]))
##     uniplane.appendChild(strip_elem)

  def saveAs(self,filepath):
    self.fileName = filepath
    saveAndClose()

  def save(self):
    if os.path.exists(self.fileName): os.remove(self.fileName)
    #else:
    out = self.doc.toprettyxml( "  " )
    f = file(self.fileName,'w')
    f.write(out)
    return f

  def saveAndClose(self):
    f = self.save()
    f.close()

  def close(self):
    self.saveAndClose()
    
  def getDoc(self):
    return self.doc


def initHotXml( file, badType, runId, times, inst="LAT", tagv="1.1", \
                source="LPA", trigger="Nominal" ):
  print "initializing", file
  handler = XmlHandler(file, badType)
  doc_elem = handler.getDoc().documentElement
  elem = handler.createGenericElement(inst,\
                                      time.strftime("%y%m%d%H%M%S", time.gmtime()),\
                                      runId, "HotStrips", \
                                      sys.argv[0], tagv )
  elem.appendChild( handler.createInSampleElement( \
      times[0], times[1], \
      trigger, "NA", source ) )
  doc_elem.appendChild(elem)

  return handler
  
def addHotStrips( handler, tower, badStripMap ):
  tower_elem = handler.goToTower( (tower,g_Serials[tower]) )
  for unp in range(g_nUniPlanes-1,-1,-1):
    list = badStripMap[tower][unp]
    if len(list) > 0:
      uniplane_elem = handler.goToUniplane( tower_elem, unp, '1',\
                                            'false', 'true', 'true' )
      handler.addStripList(uniplane_elem,list)
    else:
      uniplane_elem = handler.goToUniplane( tower_elem, unp, '1',\
                                            'false', 'false', 'false' )
  


def writeHotXml( badStripMap, runId, topDir=None, times=None, timeStamp=None, \
                 tagv="1.1", source="LPA", trigger="Nominal" ):
  badType = 'hot'
  if times == None:
    times = (time.strftime("%a %m/%d %Y, %H:%M GMT", time.gmtime()), \
             time.strftime("%a %m/%d %Y, %H:%M GMT", time.gmtime()) )

  if timeStamp == None:
    timeStamp = time.strftime("%y%m%d", time.gmtime() )

  if topDir == None:
    topDir = os.path.join( os.environ.get('TKRDATA'), 'TE403' )
  topDir = os.path.join(topDir, timeStamp)
  if not os.path.exists( topDir ): os.mkdir( topDir )

  fileName = timeStamp +'_HotStrips_TKR.xml'
  file = os.path.join( topDir, fileName )
  handler = initHotXml( file, badType, runId, times, "LAT", tagv, source, trigger )
  
  for tower in range(g_nTowers):
    hwserial = g_Serials[tower]
    fileName = timeStamp +'_HotStrips_' + hwserial + '.xml'
    fpath = os.path.join(topDir, hwserial)
    if not os.path.exists(fpath): os.mkdir( fpath )
    file = os.path.join( fpath, fileName )
    thandler = initHotXml( file, badType, runId, times, "TWR", tagv, source, trigger )
    addHotStrips( thandler, tower, badStripMap )
    addHotStrips( handler, tower, badStripMap )
    thandler.saveAndClose()
    
  handler.saveAndClose()


#***********************************
# write to LATC TFE file
#***********************************
def writeTFEXml( badStripMap, runId, topDir=None, times=None, timeStamp=None, tagv="1.1" ):
  badType = 'hot'
  if times == None:
    times = (time.strftime("%a %m/%d %Y, %H:%M GMT", time.gmtime()), \
             time.strftime("%a %m/%d %Y, %H:%M GMT", time.gmtime()) )

  if timeStamp == None:
    timeStamp = time.strftime("%y%m%d", time.gmtime() )

  if topDir == None:
    topDir = os.path.join( os.environ.get('TKRDATA'), 'TE403' )
  topDir = os.path.join(topDir, timeStamp)
  if not os.path.exists( topDir ): os.mkdir( topDir )
  fileName = timeStamp +'_intDefaultLatc_TFE.xml'
  file = os.path.join( topDir, fileName )
  print "saving", file
  handler = XmlHandler(file, "TFE")
  doc_elem = handler.getDoc().documentElement
  
  doc_elem.setAttribute( "latteInputFile", runId )
  doc_elem.setAttribute( "convertDate", \
                         time.strftime("%a %b %d, %Y %H:%M:%S GMT", \
                                       time.gmtime()))
  doc_elem.setAttribute( "convertUser", os.getenv('USER') )
  doc_elem.setAttribute( "convertPC", os.getenv('HOSTNAME') )
  doc_elem.setAttribute( "latteInputType", "Hot Strip xml File" )
  doc_elem.setAttribute( "latteCfgVersion", tagv )
    
  for tower in range(g_nTowers):
    tower_elem = handler.goToTower( (tower, "TEM") )
    for unp in range(g_nUniPlanes-1,-1,-1):
      stList = badStripMap[tower][unp]
      if len(stList) > 0:
        sptElm = handler.goToSPT( tower_elem, unp, stList, stList )
  handler.saveAndClose()


def reformatStripList(stripList):
  """\brief Converts a list of channels into a string for being included in the .xml file.

  The string contains all the strip numbers, separated by a space.
  """
  String = ''
  for channel in stripList:
    String = String + str(channel) + ' '
  # Remove the last space.
  String = String[:-1]
  return String
