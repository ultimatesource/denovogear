import ez_setup

try:
    from setuptools import setup
except ImportError:
    ez_setup.use_setuptools()
    from setuptools import setup

setup(
    name='treecall',
    version='0.1.0',
    author='Ni Huang',
    author_email='nihuang@genetics.wustl.edu',
    scripts=['treecall.py'],
    py_modules=['geno', 'utils', 'tree_est'],
    install_requires=['pyvcf', 'ete2'],
    dependency_links=['https://pypi.python.org/packages/20/b6/36bfb1760f6983788d916096193fc14c83cce512c7787c93380e09458c09/PyVCF-0.6.8.tar.gz#md5=3cc70aa59e62dab7b4a85bd5a9f2e714'],
    url='https://github.com/nh3/treecall',
    license='BSD',
    description='Tree-based joint lineage inference and somatic mutation calling',
)
