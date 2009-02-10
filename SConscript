# -*- python -*-
# $Header: 
# Authors: Johann Cohen-Tanugi <cohen@slac.stanford.edu>
# Version: calibTkrUtil-02-09-00
import os
Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
swigEnv = baseEnv.Clone()

libEnv.Tool('calibTkrUtilLib', depsOnly = 1)
libEnv.AppendUnique(CPPPATH = ['src/test'])
calibTkrUtilRootcint = libEnv.Rootcint('src/tkrPyRoot/tkrPyRoot_rootcint', ['src/tkrPyRoot/tkrPyRoot.h',
			'src/tkrPyRoot/LinkDef.h'],
			includes = [''])
 
calibTkrUtil = libEnv.SharedLibrary('calibTkrUtil', listFiles(['src/*.cxx']) + ['src/tkrPyRoot/tkrPyRoot_rootcint.cxx']
				+ ['src/tkrPyRoot/tkrPyRoot.cxx'])

lib_tkrPyRoot = swigEnv.SharedLibrary('lib_tkrPyRoot', 'src/tkrPyRoot.i')

progEnv.Tool('calibTkrUtilLib')
progEnv.AppendUnique(CPPPATH = ['src/test'])
tkrRootAnalysis = progEnv.Program('tkrRootAnalysis', listFiles(['src/test/*.cxx', 'src/*.cxx', 'src/tkrPyRoot/*.cxx']))
progEnv.Tool('registerObjects', package = 'calibTkrUtil', libraries = [calibTkrUtil, lib_tkrPyRoot], 
	binaries = [tkrRootAnalysis],
	includes = listFiles(['calibTkrUtil/*.h']))




