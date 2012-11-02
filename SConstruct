import distutils.sysconfig
env = Environment(SWIGFLAGS = ['-Wall', '-c++', '-python'],
                  CPPPATH = [distutils.sysconfig.get_python_inc()],
                  SHLIBPREFIX = "")
env.Append(CPPFLAGS = "-DSELDON_DEBUG_LEVEL_4")
conf = Configure(env)
# Link to the appropriate version of Python.
for python_version in ["2.7", "2.6", "2.5", ""]:
    if conf.CheckLib("python" + python_version):
        break
env.SharedLibrary('_seldon.so', ['Seldon.cpp', 'seldon.i'])
