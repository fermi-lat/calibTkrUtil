# -*- python -*-
# $Header: 
# Authors: Johann Cohen-Tanugi <cohen@slac.stanford.edu>
# Version: calibTkrUtil-02-09-01
import os
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
swigEnv = baseEnv.Clone()

libEnv.Tool('calibTkrUtilLib', depsOnly = 1)
libEnv.AppendUnique(CPPPATH = ['src/test'])
calibTkrUtilRootcint = libEnv.Rootcint('src/tkrPyRoot/tkrPyRoot_rootcint',
                                       ['src/tkrPyRoot/tkrPyRoot.h',
                                        'src/tkrPyRoot/LinkDef.h'],
                                       includes = [''])
libEnv['rootcnt_node'] = calibTkrUtilRootcint
 
calibTkrUtil = libEnv.SharedLibrary('calibTkrUtil', listFiles(['src/*.cxx']) + ['src/tkrPyRoot/tkrPyRoot_rootcint.cxx']
				+ ['src/tkrPyRoot/tkrPyRoot.cxx'])

lib_tkrPyRoot = swigEnv.SwigLibrary('lib_tkrPyRoot', 'src/tkrPyRoot.i')

progEnv.Tool('calibTkrUtilLib')
progEnv.AppendUnique(CPPPATH = ['src/test'])
tkrRootAnalysis = progEnv.Program('tkrRootAnalysis',
                                  listFiles(['src/test/*.cxx', 'src/*.cxx',
                                             'src/tkrPyRoot/*.cxx']))
progEnv.Tool('registerTargets', package = 'calibTkrUtil',
             rootcintSharedCxts = [[calibTkrUtil, libEnv]],
             swigLibraryCxts = [[lib_tkrPyRoot, swigEnv]], 
	binaryCxts = [[tkrRootAnalysis, progEnv]],
	includes = listFiles(['calibTkrUtil/*.h']))




