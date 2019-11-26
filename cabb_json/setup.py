from setuptools import setup

# This is the analysis-idv Python library.
# Jamie Stevens 2019
# ATCA Senior Systems Scientist
# Jamie.Stevens@csiro.au

setup(name='analysis_idv',
      version='1.1',
      description='Analysis tools for fast variable sources',
      url='https://github.com/ste616/analysis-idv',
      author='Jamie Stevens', author_email='Jamie.Stevens@csiro.au',
      license='MIT', packages=[ 'analysis_idv' ],
      install_requires=[
          'astropy', 'docopt', 'matplotlib'
      ],
      zip_safe=False)

