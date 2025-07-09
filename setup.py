from setuptools import setup

setup(
    name='iinuclear',
    version='0.1',
    author='Sebastian Gomez',
    author_email='sebastian.gomez@austin.utexas.edu',
    description='Functions to determine whether a transient is nuclear.',
    url='https://github.com/gmzsebastian/iinuclear',
    license='MIT',
    python_requires='>=3.6',
    packages=['iinuclear'],
    license_files=["LICENSE"],
    include_package_data=True,
    package_data={'iinuclear': ['ref_data/*']},
    install_requires=[
        'numpy',
        'matplotlib',
        'astropy',
        'scipy',
        'astroquery',
        'emcee',
        'alerce'
    ]
)
