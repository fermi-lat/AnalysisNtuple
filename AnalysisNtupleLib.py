# $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/AnalysisNtupleLib.py,v 1.6 2012/10/03 02:58:29 lsrea Exp $
def generate(env, **kw):
     if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AnalysisNtuple'])
        if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
	    env.Tool('findPkgPath', package = 'AnalysisNtuple') 
     env.Tool('OverlayEventLib')
     env.Tool('EventLib')
     env.Tool('LdfEventLib')
     env.Tool('ntupleWriterSvcLib')
     env.Tool('OnboardFilterTdsLib')
     env.Tool('CalibDataLib')
     env.Tool('CalUtilLib')
     #env.Tool('RootIoLib')
     env.Tool('GlastSvcLib')
     env.Tool('AcdUtilLib')
     env.Tool('astroLib')
     env.Tool('facilitiesLib')
     env.Tool('addLibrary', library=env['gaudiLibs'])
     env.Tool('addLibrary', library=env['rootLibs'])
     env.Tool('addLibrary', library=env['clhepLibs'])

# For packages which only need to provide access to header files not libs
# This is in support of non-installed headers on Windows
     if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'TkrUtil')
        env.Tool('findPkgPath', package = 'CalUtil')
def exists(env):
    return 1;
