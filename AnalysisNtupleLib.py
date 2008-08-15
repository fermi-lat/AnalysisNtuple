# $Header$
def generate(env, **kw):
     if not kw.get('depsOnly', 0):
          env.Tool('addLibrary', library = ['AnalysisNtuple'])
     env.Tool('addLibrary', library = ['dl'])
     env.Tool('addLibrary', library = ['pthread'])

def exists(env):
    return 1;
