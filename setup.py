from setuptools import setup, find_packages
setup(
  name = 'parsevasp',
  packages = find_packages(exclude=['test']),
  include_package_data=True,
  version = '0.3.0',
  description = 'A general parser for VASP',
  author = 'Espen Flage-Larsen',
  author_email = 'espen.flage-larsen@sintef.no',
  license = 'MIT',
  url = 'https://github.com/espenfl/parsevasp',
  download_url = 'https://github.com/espenfl/parsevasp/archive/0.3.0.tar.gz',
  keywords = ['VASP', 'parser', 'python', 'xml'],
  classifiers = ['Development Status :: 4 - Beta',
	         'Intended Audience :: Science/Research',
	         'Topic :: Scientific/Engineering :: Physics',
	         'License :: OSI Approved :: MIT License',
	         'Programming Language :: Python :: 2.7',
                 'Programming Language :: Python :: 3.6'],
  install_requires=['numpy', 'lxml', 'six', 'future', 'pyyaml'],
  extras_require={'dev': ['pytest', 'pytest-cov', 'click']}
)
