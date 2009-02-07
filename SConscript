# -*-python-*-
# $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/SConscript,v 1.15 2009/02/06 21:30:30 glastrm Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: AnalysisNtuple-02-54-00

import platform

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
tempEnv = baseEnv.Clone()

libEnv.Tool('AnalysisNtupleLib', depsOnly=1)
libEnv.AppendUnique(CPPDEFINES = ['GLEAM'])
libEnv.AppendUnique(CPPDEFINES = ['__i386'])
libEnv.AppendUnique(CPPDEFINES = ['EFC_FILTER'])
libEnv.AppendUnique(CPPDEFINES= ['f2cFortran'])
libEnv.AppendUnique(FORTRANFLAGS = ['-Wno-globals'])

if baseEnv['PLATFORM'] != 'win32':
	libEnv.AppendUnique(CPPDEFINES = ['MSG_ST'])
	libEnv.AppendUnique(CPPDEFINES = ['__unix__'])
	libEnv.AppendUnique(CPPDEFINES = ['WUT_IMPLEMENTATION_PX'])

if baseEnv['PLATFORM'] == 'win32':
	libEnv.AppendUnique(CPPDEFINES = ['WIN32'])
AnalysisNtuple = libEnv.SharedLibrary('AnalysisNtuple', listFiles(['src/*.cxx','src/Dll/*.cxx']))


progEnv.Tool('AnalysisNtupleLib')
if baseEnv['PLATFORM'] != 'win32':
	progEnv.Tool('addLibrary', library = ['dl'])
	progEnv.Tool('addLibrary', library = ['pthread'])

test_AnalysisNtuple = progEnv.GaudiProgram('test_AnalysisNtuple', ['src/AnalysisNtupleAlg.cxx'],test=1)

progEnv.Tool('registerObjects', package = 'AnalysisNtuple', libraries = [AnalysisNtuple],
	     testApps = [test_AnalysisNtuple], includes = listFiles(['AnalysisNtuple/*.h']))




