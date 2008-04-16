import os, sys, array
sys.path = [( os.path.join( os.environ.get( "ROOTSYS" ), "lib" ) )] + sys.path

import ROOT
import tkrUtils

ROOT.gSystem.Load("tkrPyRoot")

#
# initialization
# 
keys = [ ("TOT_Peak", "layer", "D"), ("TOT_LWidth", "layer", "D"), \
         ("TOT_GSigma", "layer", "D"), ("TOT_Entries", "layer", "D" ), \
         ("TOT_PeakError", "layer", "D"), ("TOT_LWidthError", "layer", "D"), \
         ("TOT_GSigmaError", "layer", "D"), ("TOT_FitProb", "layer", "D"), \
         ("TOT_FracLowTOT", "layer", "D" ) ]
nelem = {}
nelem["layer"] = tkrUtils.g_nTowers * tkrUtils.g_nUniPlanes
typecode = { "C":"c", "B":"b", "b":"B", "S":"i", "s":"I", "F":"f", "D":"d" }

rname = "/nfs/slac/work/htajima/outputs/TE403_080416-072720.root"
tname = "tkrMonitorTrees.root"

#
# setup stuuf
#
values = {}
for (key, type, vtype) in keys:
  values[key] = array.array( typecode[vtype], [-1]*nelem[type] )

ffit = ROOT.defLangau( "langau", 0, 30 )
ffit.SetParNames( "LWidth", "MP", "Area", "GSigma" )

rf = ROOT.TFile( rname )
tree = rf.FindObjectAny( "timeStamps" )
tStartTime = array.array( 'd', [0] )
tEndTime = array.array( 'd', [0] )
tree.SetBranchAddress( "startTime", tStartTime )
tree.SetBranchAddress( "endTime", tEndTime )
startTime = array.array( 'd', [-1] )
endTime = array.array( 'd', [1] )
for i in range(tree.GetEntries()):
  tree.GetEntry( i )
  if startTime[0] < 0 or tStartTime[0] < startTime[0]:
    startTime[0] = tStartTime[0]
  if tEndTime[0] > endTime[0]: endTime[0] = tEndTime[0]
print startTime, endTime
#
# loop towers/layers ann fit TOT distributions
#
for tower in range(tkrUtils.g_nTowers):
  for unp in range(tkrUtils.g_nUniPlanes):
    ielm = unp + tower*tkrUtils.g_nUniPlanes
    lname = tkrUtils.g_layerNames[unp]
    name = "chargeT%d%s" % (tower, lname)
    chargeHist = rf.FindObjectAny( name )

    entries = chargeHist.Integral()
    mean = chargeHist.GetMean()
    rms = chargeHist.GetRMS()
    values["TOT_Entries"][ielm] = entries

    if entries<200 or mean==0.0 or rms==0.0:
      print "T%d %s, Entries %.0f, Mean: %.1f, RMS: %.1f skipped" \
            % (tower, lname, entries, mean, rms)
      continue

    binWidth = chargeHist.GetBinWidth( 2 )
    bin = int(mean*0.5/binWidth) + 1
    fracBadTot = chargeHist.Integral(1,bin) / chargeHist.Integral()
    values["TOT_FracLowTOT"][ielm] = fracBadTot
    #m_fracBatTot.Fill( fracBadTot )

    lowLim = mean - 1.4 * rms
    if fracBadTot > 0.05 and lowLim < mean*0.5:
      lowLim = mean*0.5
      print "WARNIG, large bad TOT fraction: %.3f, T%d %s." \
            % ( fracBadTot, tower, lname )

    ffit.SetParLimits( 0, 0.10, 0.5 )
    ffit.SetParLimits( 1, 0.0, mean*2 )
    ffit.SetParLimits( 2, 0.0, entries*1.5 )
    ffit.SetParLimits( 3, rms*0.2, rms*1.5 )
    ffit.SetParameters( rms*0.2, mean*0.75, entries*0.1, rms*0.4 )
    #ffit.SetRange( lowLim, mean+2.5*rms )
    #ffit.FixParameter( 4, m_RSigma )
    #ffit.FixParameter( 5, m_GFrac )
    chargeHist.Fit( "langau", "RBQ", "", lowLim, mean+2.5*rms )
    ROOT.gPad.Update()

    #0:width 1:peak 2:total entries 3:width(sigma)
    par = ffit.GetParameters()
    error = ffit.GetParErrors()
    values["TOT_LWidth"][ielm] = par[0]
    values["TOT_Peak"][ielm] = par[1]
    values["TOT_GSigma"][ielm] = par[3]
    values["TOT_LWidthError"][ielm] = error[0]
    values["TOT_PeakError"][ielm] = error[1]
    values["TOT_GSigmaError"][ielm] = error[3]
    values["TOT_FitProb"][ielm] = ffit.GetProb()

rf.Close()

#
# save TTree
#
tf = ROOT.TFile( tname, "RECREATE")
tree = ROOT.TTree( "tkrMonitor", "tkrMonitor" )

for (key, type, vtype) in keys:
  if type == "layer":
    leaflist = "%s[%d][%d]/%s" \
               % (key, tkrUtils.g_nTowers, tkrUtils.g_nUniPlanes, vtype )
  elif type == "tower":
    leaflist = "%s[%d]/%s" % (key, tkrUtils.g_nTowers, vtype )
  tree.Branch(key, values[key], leaflist)

tree.Branch( "startTime", startTime, "startTime/D" )
tree.Branch( "endTime", endTime, "endTime/D" )

tree.Fill()
tree.Write()
tf.Close()
