.. image:: images/iinuclear.png
   :alt: iinuclear logo
   :width: 200px
   :align: center
   
Is it nuclear?
==============

``iinuclear`` is a Python package designed to determine whether a transient is nuclear. The code will determine the probability of
whether the location of the transient coincides with the center of its host galaxy. This determination is crucial for various areas
of time-domain astronomy such as TDEs, AGN, or other nuclear transient phenomena.

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   Installation
   Usage
   Methodology

Quick Start
-----------

To determine if a transient is nuclear::

    from iinuclear import isit

    # Using an IAU name
    isit("2018hyz")

    # Using a ZTF name
    isit("ZTF18acpdvos")

    # Using coordinates in degrees
    isit(151.711964138, 1.69279894089)

Key Features
------------

* Automated retrieval of transient positions from ZTF
* Statistical analysis of positional coincidence
* Visualization tools including:
    - PS1 image cutouts
    - Position scatter plots
    - Separation histograms
    - Statistical significance plots
* Quantitative metrics for nuclear classification

Requirements
------------

* Python 3.7 or later
* Access to TNS API (for IAU name queries)
* Having the ``alerce`` API package installed
* Having the ``astroquery`` package installed

See :ref:`installation` for detailed setup instructions.

Citation
--------

If you use ``iinuclear`` in your research, please cite [reference to be added].

License & Attribution
---------------------

Copyright 2024 Sebastian Gomez and `contributors <https://github.com/gmzsebastian/iinuclear/graphs/contributors>`_.
``iinuclear`` is free software made available under the MIT License.