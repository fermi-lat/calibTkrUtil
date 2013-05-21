# -*- python -*-
# $Header: 
# Authors: Johann Cohen-Tanugi <cohen@slac.stanford.edu>
# Version: calibTkrUtil-02-10-01
Import('baseEnv')
Import('listFiles')
Import('packages')
if baseEnv['PLATFORM'] != "win32":
  progEnv = baseEnv.Clone()
  libEnv = baseEnv.Clone()
  swigEnv = baseEnv.Clone()

  locIncs = listFiles(['src/tkrPyRoot/*.h'])
  libEnv.Tool('addLinkDeps', package='calibTkrUtil', toBuild='rootlib')
  libEnv.AppendUnique(CPPPATH = ['src/test'])
  calibTkrUtilRootcint = libEnv.Rootcint('src/tkrPyRoot/tkrPyRoot_rootcint',
                                         ['src/tkrPyRoot/tkrPyRoot.h',
                                          'src/tkrPyRoot/LinkDef.h'],
                                         includes = [''],
                                         localIncludes = locIncs,
                                         packageName = 'calibTkrUtil')
  libEnv['rootcint_node'] = calibTkrUtilRootcint

  tkrPyRootLib = libEnv.RootDynamicLibrary('tkrPyRoot',
                                           listFiles(['src/*.cxx'])
                                           + ['src/tkrPyRoot/tkrPyRoot_rootcint.cxx']
                                           + ['src/tkrPyRoot/tkrPyRoot.cxx'])

  tkrPyRootSwigLib = swigEnv.SwigLibrary('_tkrPyRoot', 'src/tkrPyRoot.i')

  progEnv.Tool('calibTkrUtilLib')
  progEnv.AppendUnique(CPPPATH = ['src/test'])
  tkrRootAnalysis = progEnv.Program('tkrRootAnalysis',
                                    listFiles(['src/test/*.cxx']))
  #listFiles(['src/test/*.cxx', 'src/*.cxx',
  progEnv.Tool('registerTargets', package = 'calibTkrUtil',
               rootcintSharedCxts = [[tkrPyRootLib, libEnv]],
               swigLibraryCxts = [[tkrPyRootSwigLib, swigEnv]], 
               binaryCxts = [[tkrRootAnalysis, progEnv]],
               includes = listFiles(['calibTkrUtil/*.h']))




