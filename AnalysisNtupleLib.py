# $Header: /nfs/slac/g/glast/ground/cvs/GlastRelease-scons/AnalysisNtuple/AnalysisNtupleLib.py,v 1.3 2008/11/06 22:09:30 ecephas Exp $
def generate(env, **kw):
     if not kw.get('depsOnly', 0):
          env.Tool('addLibrary', library = ['AnalysisNtuple'])
     env.Tool('OverlayEventLib')
     env.Tool('EventLib')
     env.Tool('LdfEventLib')
     #env.Tool('TkrUtilLib')
     env.Tool('ntupleWriterSvcLib')
     #env.Tool('OnboardFilterLib')
     env.Tool('OnboardFilterTdsLib')
     #env.Tool('FluxSvcLib')
     env.Tool('CalibDataLib')
     #env.Tool('G4PropagatorLib')
     #####env.Tool('addLibrary', library=env['obfLibs'])
     #env.Tool('RootIoLib')
     env.Tool('GlastSvcLib')
     env.Tool('AcdUtilLib')
     env.Tool('astroLib')
     env.Tool('facilitiesLib')
     env.Tool('addLibrary', library=env['gaudiLibs'])
     env.Tool('addLibrary', library=env['rootLibs'])
     env.Tool('addLibrary', library=env['clhepLibs'])
def exists(env):
    return 1;
