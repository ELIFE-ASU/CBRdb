from setuptools import setup, find_packages

setup(
    name='CBRdb',
    version='1.4.0',
    author='Louie Slocombe, Camerian Millsaps, Reza Shahjahan, Kamesh Narasimhan, and Sara Walker',
    author_email='louies@hotmail.co.uk',
    description='A curated biochemical database that integrates and refines data from KEGG and ATLAS databases to support precise analyses of biochemical reaction data.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/ELIFE-ASU/CBRdb',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.9',
    install_requires=[
        'numpy',
        'sympy',
        'matplotlib',
        'networkx',
        'pandas',
        'rdkit',
        'chempy',
        'requests',
        'urllib3',
        'chemparse',
        'ase',
        'mace',
        'pymatgen',
        'dimorphite_dl',
        'drfp',
        'equilibrator-api',
        # 'rxnmapper',
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
        ],
    },
)
