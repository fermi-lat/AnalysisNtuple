# -*-python-*-
# $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/AnalysisNtuple/SConscript,v 1.3 2008/08/19 14:08:07 ecephas Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
import platform

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

libEnv.Tool('AnalysisNtupleLib', depsOnly=1)
libEnv.AppendUnique(CPPDEFINES = ['GLEAM'])
libEnv.AppendUnique(CPPDEFINES = ['__i386'])
libEnv.AppendUnique(CPPDEFINES = ['EFC_FILTER'])

if baseEnv['PLATFORM'] == 'Linux':
	libEnv.AppendUnique(CPPDEFINES = ['MSG_ST'])
	libEnv.AppendUnique(CPPDEFINES = ['__unix__'])
	libEnv.AppendUnique(CPPDEFINES = ['WUT_IMPLEMENTATION_PX'])

if baseEnv['PLATFORM'] == 'win32':
	libEnv.AppendUnique(CPPDEFINES = ['WIN32'])

AnalysisNtuple = libEnv.SharedLibrary('AnalysisNtuple', listFiles(['src/*.cxx','src/Dll/*.cxx']))


progEnv.Tool('AnalysisNtupleLib')
if baseEnv['PLATFORM'] == 'Linux':
	progEnv.Tool('addLibrary', library = ['dl'])
	progEnv.Tool('addLibrary', library = ['pthread'])

test_AnalysisNtuple = progEnv.Program('test_AnalysisNtuple', ['src/test/AnalysisNtupleAlg.cxx'])

progEnv.Tool('registerObjects', package = 'AnalysisNtuple', libraries = [AnalysisNtuple],
	     testApps = [test_AnalysisNtuple], includes = listFiles(['AnalysisNtuple/*.h']))
