#$Header: 
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['tkrPyRoot'])
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('commonRootDataLib')
    env.Tool('reconRootDataLib')
    env.Tool('digiRootDataLib')
def exists(env):
    return 1;
