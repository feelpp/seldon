import distutils.sysconfig
env = Environment(SWIGFLAGS = ['-Wall', '-c++', '-python'],
                  CPPPATH = [distutils.sysconfig.get_python_inc()],
                  SHLIBPREFIX = "")
env.SharedLibrary('_seldon.so', ['Seldon.cpp', 'seldon.i'])
