from setuptools import setup, find_packages

with open("README.md") as f:
    long_description = f.read()

setup(
    name='prenergy',
    version='0.1.0',
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=[
        'pandas>=1.0.0',
        'numpy>=1.18.0',
        'matplotlib>=3.1.0',
        'seaborn>=0.10.0',
        'statsmodels>=0.11.0',
        'mol2vec>=0.1',
        'rdkit>=2020.03.1',
        'gensim>=3.8.0',
        'scikit-learn>=0.22.0',
        'shap>=0.35.0'
    ],
    author='João V. V. Cassiano',
    author_email='joaocassianox7x@gmail.com',
    description='A package to predict adsorption energies of molecules on surfaces.',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
