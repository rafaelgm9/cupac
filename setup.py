import os
import subprocess
import shutil
import sys
import glob
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

CMAKE_EXE = os.environ.get('CMAKE_EXE', shutil.which('cmake'))


def check_for_cmake():
    if not CMAKE_EXE:
        print('cmake executable not found. '
              'Set CMAKE_EXE environment or update your path')
        sys.exit(1)


class CMakeExtension(Extension):
    """
    setuptools.Extension for cmake
    """

    def __init__(self, name, sourcedir='', src=[], hdr=[]):
        check_for_cmake()
        Extension.__init__(self, name, sources=src, depends=hdr)
        self.sourcedir = os.path.abspath(sourcedir)
        self.target = name.split('.')[-1]


class CMakeBuildExt(build_ext):
    """
    setuptools build_exit which builds using cmake & make
    You can add cmake args with the CMAKE_COMMON_VARIABLES environment variable
    """

    def build_extension(self, ext):
        check_for_cmake()
        if isinstance(ext, CMakeExtension):
            output_dir = os.path.abspath(
                os.path.dirname(self.get_ext_fullpath(ext.name)))

            build_type = 'Debug' if self.debug else 'Release'
            cmake_args = [CMAKE_EXE,
                          ext.sourcedir,
                          '-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + output_dir,
                          '-DCMAKE_BUILD_TYPE=' + build_type]
            cmake_args.extend(
                [x for x in
                 os.environ.get('CMAKE_COMMON_VARIABLES', '').split(' ')
                 if x])

            env = os.environ.copy()
            if not os.path.exists(self.build_temp):
                os.makedirs(self.build_temp)
            subprocess.check_call(cmake_args,
                                  cwd=self.build_temp,
                                  env=env)
            subprocess.check_call(['make', '-j', ext.target],
                                  cwd=self.build_temp,
                                  env=env)
            print()
        else:
            super().build_extension(ext)


sources = glob.glob(os.path.join('src', '*.cu'))
headers = glob.glob(os.path.join('src', '*.cuh'))

# module = Extension('cupac',
#                    define_macros=[('MAJOR_VERSION', '1'),
#                                   ('MINOR_VERSION', '0')],
#                    sources=sources,
#                    depends=headers)

module = CMakeExtension('cupac', src=sources, hdr=headers)

setup(name='cupac',
      version='1.0',
      description='Module for counting pairs of particles in a cubic volume.',
      author='Rafael Garcia',
      author_email='rgarciamar@email.arizona.edu',
      ext_modules=[module],
      cmdclass={'build_ext': CMakeBuildExt},
      packages=['cupac'],
      install_requires=['numpy'],
      zip_safe=False)
