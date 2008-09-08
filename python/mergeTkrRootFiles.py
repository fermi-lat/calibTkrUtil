# import python standard libraries
import sys, os, math, time, array

sys.path = [( os.path.join( os.environ.get( "ROOTSYS" ), "lib" ) )] + sys.path

# import thrid party libraries
import ROOT

classNames = ["TH1F", "TH2F", "TProfile", "TTree"]


#**********************************
# find obejects in the directory
#**********************************
def findObjects( dir, (dirMap, objMap), treeList ):
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
      findObjects( ROOT.gDirectory, dirMap[name], treeList )
      ROOT.gDirectory.cd( ".." )
    elif cname in classNames:
      if not objMap.has_key(name): objMap[name] = cname
      if cname == "TTree": treeList[name] = None
    else: print "unknown object: %s, %s" % (name, cname)
    obj = tlist.After( obj )
  tlist.Delete()
  #print dir.GetName(), len(objMap)


#**********************************
# merge obejects in the directory
#**********************************
def mergeObjects( objMap, dir, inFiles ):
  #print "+++ %s +++" % dir.GetName(), objMap.keys()
  #print "# of objects in %s: %d" % ( dir.GetName(), len(objMap) )
  for name in objMap.keys():
    cname = objMap[name]
    #print dir.GetName(), name, cname
    if name == "totInfo": continue
    elif cname == "TTree":
      continue
      print "merge", name
      dir.cd()
      tlist = ROOT.TList()
      for inFile in inFiles:
        tree = inFile.FindObjectAny( name )
        if tree: tlist.Add( tree )
      dir.cd()
      try:
        tree = ROOT.TTree.MergeTrees( tlist )
        tree.Write()
        tlist.Delete()
        print name, "saved. size:", tree.GetEntries()
      except:
        print "merge erorr:", name
    else:
      dir.cd()
      hist = None
      for inFile in inFiles:
        try: thist = inFile.FindObjectAny( name )
        except: continue
        if thist != None:
          if hist == None: hist = thist.Clone( name )
          else: hist.Add( thist )
          thist.Delete()
            
      if hist != None: hist.Write()
      else: print "%s missing." % name


#**********************************
# find obejects in the directory
#**********************************
def mergeDirectories( dirMap, dir, inFiles, path="" ):
  #print "--- %s ---" % path, dirMap.keys()
  #print "# of dirs in %s: %d" % ( dir.GetName(), len(dirMap) )
  for key in dirMap.keys():
    tdir = dir.mkdir( key )
    ( tdirMap, objMap ) = dirMap[key]
    print "***** %s: %d %d *****" % (path+":"+key, len(objMap), len(tdirMap) )
    mergeObjects( objMap, tdir, inFiles )
    mergeDirectories( tdirMap, tdir, inFiles, path+":"+key )
  dir.cd()
  print "==== %s completed. ====" % path



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
  t0 = time.time()
  maxFiles = 500
  tdir = "."
  oname = None
  fnames = []
  if len(sys.argv) > 2:
    iarg = 1
    while iarg < len(sys.argv):
      if sys.argv[iarg] == "-tdir":
        iarg += 1
        tdir = sys.argv[iarg]
        print "temporary directory:", tdir
      elif sys.argv[iarg] == "-max":
        iarg += 1
        maxFiles = int(sys.argv[iarg])
        print "max files:", maxFiles
      elif oname == None:
        oname = sys.argv[iarg]
        print "output file:", oname
      else:
        iname = sys.argv[iarg]
        if not os.path.exists( iname ):
          print "file %s does not exist." % iname
          sys.exit()
        fnames.append( iname )
      iarg += 1
        
  if oname == None or len(fnames)==0:
    print "usage: mergeTkrRootFiles.py [options] [output] [input1] [input2]..."
    print " [options] include -tdir for temporary directory"
    print " and -max for maximum files per process"
    sys.exit()

  #
  # analyze directory structure and find all objects in all files
  #
  treeList = {}
  (dirMap, objMap) = ( {}, {} )
  for iname in fnames:
    inFile = ROOT.TFile( iname )
    findObjects( inFile, (dirMap, objMap), treeList )
    inFile.Close()

  #
  # if total # of input files are larger than maximum
  # input files are divided.
  #
  inFiles = []
  tnames = []
  if len(fnames) > maxFiles:
    ts = time.strftime("%y%m%d%H%M%S", time.gmtime() )
    ndiv = int(len(fnames)/maxFiles) + 1
    for idiv in range(ndiv):
      kf = idiv*len(fnames)/ndiv
      kl = (idiv+1)*len(fnames)/ndiv
      print "merge", fnames[kf:kl]
      print "# of files", len(fnames[kf:kl])
      tFiles = []
      for fname in fnames[kf:kl]:
        tFiles.append( ROOT.TFile( fname ) )
      tname = os.path.join( tdir, "temp-%s-%d.root"%(ts,idiv) )
      tnames.append( tname )
      tRoot = ROOT.TFile( tname, "RECREATE" )
      mergeObjects( objMap, tRoot, tFiles )
      mergeDirectories( dirMap, tRoot, tFiles, "TEMP%d"%idiv )
      tRoot.Close()
      for tFile in tFiles: tFile.Close()
      print "Elapsed time: %.1f seconds" % (time.time()-t0) 
    #fnames = tnames

  #
  # produce final output
  #
  for key in treeList.keys():
    treeList[key] = ROOT.TList()
  rFiles = []
  for fname in fnames:
    rFiles.append( ROOT.TFile( fname ) )
    # trees need to be picked up here to work properly
    for key in treeList.keys():
      tree = rFiles[-1].FindObjectAny( key )
      if tree:
        treeList[key].Add( tree )
        #print key, tree.GetEntries()
        
  print "# of trees to merge", len(rFiles)
  outRoot = ROOT.TFile( oname, "RECREATE" )
  print "open:", oname
  for key in treeList.keys():
    try:
      tree = ROOT.TTree.MergeTrees( treeList[key] )
      tree.Write()
      treeList[key].Delete()
      print key, "saved. size:", tree.GetEntries()
    except:
      print "merge erorr:", key
      sys.exit()
  if len(inFiles) == 0:
    print "# of file to merge", len(rFiles)
    mergeObjects( objMap, outRoot, rFiles )
    mergeDirectories( dirMap, outRoot, rFiles, "OUT" )
  else:
    print "# of file to merge", len(inFiles)
    sys.exit()
    mergeObjects( objMap, outRoot, inFiles )
    mergeDirectories( dirMap, outRoot, inFiles, "OUT" )

  print "close output file:", outRoot.GetName()
  outRoot.Close()
  print "close input files"
  for inFile in inFiles: inFile.Close()
  for rFile in rFiles: rFile.Close()

  if len(tnames)>0:
    print "clean up temporary files"
    for tname in tnames:
      os.remove( tname )
  
  print "Elapsed time: %.1f seconds" % (time.time()-t0) 
