# import python standard libraries
import sys, os

sys.path = [( os.path.join( os.environ.get( "ROOTSYS" ), "lib" ) )] + sys.path

# import thrid party libraries
import ROOT

classNames = ["TH1F", "TH2F", "TProfile", "TTree"]


#**********************************
# find obejects in the directory
#**********************************
def findObjects( dir, (dirMap, objMap) ):
  tlist = dir.GetListOfKeys()
  obj = tlist.First()
  while obj:
    name = obj.GetName()
    cname = obj.GetClassName()
    #print name, cname
    if cname[:4] == "TDir":
      if not dirMap.has_key( name ):
        dirMap[name] = ( {}, {} )
      ROOT.gDirectory.cd( name )
      findObjects( ROOT.gDirectory, dirMap[name] )
      ROOT.gDirectory.cd( ".." )
    elif cname in classNames:
      if not objMap.has_key(name): objMap[name] = cname
    else: print "unknown object: %s, %s" % (name, cname)
    obj = tlist.After( obj )
  print dir.GetName(), len(objMap)


#**********************************
# merge obejects in the directory
#**********************************
def mergeObjects( objMap, dir, inFiles ):
  #print "+++ %s +++" % dir.GetName(), objMap.keys()
  print "# of objects in %s: %d" % ( dir.GetName(), len(objMap) )
  for name in objMap.keys():
    cname = objMap[name]
    #print dir.GetName(), name, cname
    if cname == "TTree":
      tlist = ROOT.TList()
      for inFile in inFiles:
        tlist.Add( inFile.FindObjectAny( name ) )
      dir.cd()
      tree = ROOT.TTree.MergeTrees( tlist )
      tree.Write()
    else:
      dir.cd()
      hist = None
      for inFile in inFiles:
        try: thist = inFile.FindObjectAny( name )
        except: continue
        if thist != None:
          if hist == None: hist = thist.Clone( name )
          else: hist.Add( thist )
            
      if hist != None: hist.Write()
      else: print "%s missing." % name


#**********************************
# find obejects in the directory
#**********************************
def mergeDirectories( dirMap, dir, inFiles ):
  #print "--- %s ---" % dir.GetName(), dirMap.keys()
  print "# of dirs in %s: %d" % ( dir.GetName(), len(dirMap) )
  for key in dirMap.keys():
    print "***** %s *****" % key
    tdir = dir.mkdir( key )
    ( tdirMap, objMap ) = dirMap[key]
    mergeObjects( objMap, tdir, inFiles )
    mergeDirectories( tdirMap, tdir, inFiles )
  dir.cd()
  print "==== %s completed. ====" % dir.GetName()



#**************************************************************************
#
# main
#
#**************************************************************************
#**************************************************************************
if __name__ == '__main__':

  #
  # decode command arguments
  # [out root file] [input root files]
  #
  if len(sys.argv) > 2:
    oname = sys.argv[1]
    outRoot = ROOT.TFile( oname, "RECREATE" )
    inFiles = []
    (dirMap, objMap) = ( {}, {} )
    for iname in sys.argv[2:]:
      if not os.path.exists( iname ):
        print "file %s does not exist." % iname
        sys.exit()
      inFiles.append( ROOT.TFile( iname ) )
      findObjects( inFiles[-1], (dirMap, objMap) )
  else:
    print "usage: mergeTkrRootFiles.py [output] [input1] [input2]..."
    sys.exit()

  mergeObjects( objMap, outRoot, inFiles )
  mergeDirectories( dirMap, outRoot, inFiles )

  outRoot.Close()
  for inFile in inFiles: inFile.Close()
  
