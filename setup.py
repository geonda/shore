from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()


    
setup(
    name='shore-pkg-geonda',  # Replace with your package name
    version='0.1.4',          # Initial version
    packages=find_packages(),  # Automatically find packages in the directory
    install_requires=required,
    author='Andrey Geondzhian',       # Your name
    author_email='a.geondzhian@gmail.com',  # Your email
    description='python wrapper for OCEAN-FEFF',
    long_description=open('README.md').read(),  # Read from README file for long description
    long_description_content_type='text/markdown',  # Format of the long description
    url='https://github.com/geonda/shore',  # URL to your project
)
