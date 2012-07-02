# -*-python-*-
# $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/SConscript,v 1.93 2012/05/29 22:59:36 lsrea Exp $
# Authors: Leon Rochester <lsrea@slac.stanford.edu>
# Version: AnalysisNtuple-02-63-01

import platform

Import('baseEnv')
Import('listFiles')
Import('packages')
progEnv = baseEnv.Clone()
libEnv = baseEnv.Clone()
tempEnv = baseEnv.Clone()

libEnv.Tool('addLinkDeps', package='AnalysisNtuple', toBuild='component')
libEnv.AppendUnique(CPPDEFINES = ['GLEAM'])
libEnv.AppendUnique(CPPDEFINES = ['__i386'])
libEnv.AppendUnique(CPPDEFINES = ['EFC_FILTER'])


if baseEnv['PLATFORM'] != 'win32':
	libEnv.AppendUnique(CPPDEFINES = ['MSG_ST'])
	libEnv.AppendUnique(CPPDEFINES = ['__unix__'])
	libEnv.AppendUnique(CPPDEFINES = ['WUT_IMPLEMENTATION_PX'])
else:
	libEnv.AppendUnique(CPPDEFINES = ['WIN32'])

AnalysisNtuple = libEnv.ComponentLibrary('AnalysisNtuple',
					 listFiles(['src/*.cxx']))


progEnv.Tool('AnalysisNtupleLib')
if baseEnv['PLATFORM'] != 'win32':
	progEnv.Tool('addLibrary', library = ['dl'])
	progEnv.Tool('addLibrary', library = ['pthread'])

test_AnalysisNtuple = progEnv.GaudiProgram('test_AnalysisNtuple', [],
					   test=1, package='AnalysisNtuple')
#                                           ['src/AnalysisNtupleAlg.cxx'],

progEnv.Tool('registerTargets', package = 'AnalysisNtuple',
             libraryCxts = [[AnalysisNtuple,libEnv]],
	     testAppCxts = [[test_AnalysisNtuple, progEnv]],
             includes = listFiles(['AnalysisNtuple/*.h']),
	     data = listFiles(['data/*.txt']),
	     jo = listFiles(['src/test/*.txt']))
