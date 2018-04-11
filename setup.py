import sys
import imp
import os
import subprocess
from setuptools import setup, Extension, find_packages
from setuptools.command.install import install


class custom_install(install):
    def run(self):
        install.run(self)

        print "Building MUMmer..."
        cwd = os.getcwd()
        syspath = sys.path[:]
        try:
            syspath.remove(cwd)
        except ValueError:
            pass
        f, pn, d = imp.find_module('strainest', syspath)
        module = imp.load_module('strainest', f, pn, d)

        print "MUMmer path: %s" % module.mummer_path()
        os.chdir(module.mummer_path())
        proc = subprocess.Popen(['make', 'install'],
                                stderr=subprocess.PIPE)
        _, out_stderr = proc.communicate()
        if proc.returncode:
            raise Exception(out_stderr)
        os.chdir(cwd)


setup(
    name='strainest',
    version = '1.2.2',
    packages = find_packages(),
    description='Abundance estimation of strains',
    long_description=open('README.rst').read(),
    url='https://github.com/compmetagen/strainest',
    download_url='https://github.com/compmetagen/strainest/releases',
    license='GPLv3',
    author='Davide Albanese and Claudio Donati',
    author_email='davide.albanese@fmach.it',
    maintainer='Davide Albanese',
    maintainer_email='davide.albanese@fmach.it',
    install_requires=[
        'Click >=5.1',
        'numpy>=1.7.0',
        'scipy',
        'pandas',
        'pysam>=0.12',
        'scikit-learn>=0.16.1',
        'matplotlib>=1.3.0',
        'biopython>=1.50'],
    entry_points='''
            [console_scripts]
            strainest=strainest.scripts.strainest_cmd:cli
        ''',
    use_2to3 = True,
    include_package_data = True,
    cmdclass={
        'install': custom_install,
    }
)
