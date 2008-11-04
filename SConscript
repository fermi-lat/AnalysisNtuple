# -*-python-*-
# $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/AnalysisNtuple/SConscript,v 1.5 2008/10/31 15:42:28 ecephas Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
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
ctype_fix=tempEnv.Object(target = 'ctype_fix.o', source = ['src/flux_c_type_fix.cxx'])

if baseEnv['PLATFORM'] != 'win32':
	libEnv.AppendUnique(CPPDEFINES = ['MSG_ST'])
	libEnv.AppendUnique(CPPDEFINES = ['__unix__'])
	libEnv.AppendUnique(CPPDEFINES = ['WUT_IMPLEMENTATION_PX'])

if baseEnv['PLATFORM'] == 'win32':
	libEnv.AppendUnique(CPPDEFINES = ['WIN32'])
libEnv.AppendUnique(CPPFLAGS = ['-Wl,--wrap,__ctype_b'])
AnalysisNtuple = libEnv.SharedLibrary('AnalysisNtuple', listFiles(['src/*.cxx','src/Dll/*.cxx']))


progEnv.Tool('AnalysisNtupleLib')
progEnv.AppendUnique(CPPDEFINES= ['f2cFortran'])
progEnv.AppendUnique(FORTRANFLAGS = ['-Wno-globals'])

if baseEnv['PLATFORM'] != 'win32':
	progEnv.Tool('addLibrary', library = ['dl'])
	progEnv.Tool('addLibrary', library = ['pthread'])

test_AnalysisNtuple = progEnv.Program('test_AnalysisNtuple', ['src/AnalysisNtupleAlg.cxx'] + [ctype_fix])

progEnv.Tool('registerObjects', package = 'AnalysisNtuple', libraries = [AnalysisNtuple],
	     testApps = [test_AnalysisNtuple], includes = listFiles(['AnalysisNtuple/*.h']))
