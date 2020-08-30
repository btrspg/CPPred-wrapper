from setuptools import setup

setup(
    name='CPPred-wrapper',
    version='0.1.0',
    packages=['cppred'],
    scripts=['CPPred.py'],
    include_package_data=True,
    package_data={'': ['data/*.*']},
    url='',
    license='',
    author='CHEN Yuelong',
    author_email='yuelong.chen.btr@gmail.com',
    description='''
This is a wrapper of CPPred algorithm. Download from original source, [http://www.rnabinding.com/CPPred/CPPred/CPPred.tar.gz](http://www.rnabinding.com/CPPred/CPPred/CPPred.tar.gz), would be perfect.

This wrapper is just for myself easy use.

Almost whole codes were from [http://www.rnabinding.com/CPPred/CPPred/CPPred.tar.gz](http://www.rnabinding.com/CPPred/CPPred/CPPred.tar.gz). I just modified the structure and install setup.py for myself easier use.

    '''
)
