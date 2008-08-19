# -*-python-*-
# $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/AnalysisNtuple/SConscript,v 1.2 2008/08/18 17:18:54 ecephas Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
import platform

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()

progEnv.Replace(CPPDEFINES = ['GLEAM'])
progEnv.AppendUnique(CPPDEFINES = ['__i386'])
progEnv.AppendUnique(CPPDEFINES = 'EFC_FILTER')

if baseEnv['PLATFORM'] != 'win32':
	progEnv.AppendUnique(CPPDEFINES = 'MSG_ST')
	progEnv.AppendUnique(CPPDEFINES = '__unix__')
	progEnv.AppendUnique(CPPDEFINES = 'WUT_IMPLEMENTATION_PX')
	libEnv.Tool('addLibrary', library = ['dl'])
	libEnv.Tool('addLibrary', library = ['pthread'])

if baseEnv['PLATFORM'] == 'win32':
	progEnv.AppendUnique(CPPDEFINES = 'WIN32')


libEnv.Tool('AnalysisNtupleLib', depsOnly=1)
AnalysisNtuple = libEnv.SharedLibrary('AnalysisNtuple', listFiles(['src/*.cxx','src/Dll/*.cxx']))


progEnv.Tool('AnalysisNtupleLib')
test_AnalysisNtuple = progEnv.Program('test_AnalysisNtuple', ['src/test/AnalysisNtupleAlg.cxx'])

progEnv.Tool('registerObjects', package = 'AnalysisNtuple', libraries = [AnalysisNtuple],
	     testApps = [test_AnalysisNtuple], includes = listFiles(['AnalysisNtuple/*.h']))
