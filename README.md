# Is It Nuclear?

<p align="center">
  <img src="docs/images/iinuclear.png" width="200">
</p>

`iinuclear` is a Python package designed to determine whether a transient is nuclear and whether it coincides with the center of its host galaxy.

* Documentation: [https://iinuclear.readthedocs.io/](https://iinuclear.readthedocs.io/)
* Code: [https://github.com/gmzsebastian/iinuclear](https://github.com/gmzsebastian/iinuclear)
* License: MIT

![Tests](https://github.com/gmzsebastian/iinuclear/actions/workflows/ci_tests.yml/badge.svg)
![License](http://img.shields.io/badge/license-MIT-blue.svg)
[![Coverage Status](https://coveralls.io/repos/github/gmzsebastian/iinuclear/badge.svg?branch=main)](https://coveralls.io/github/gmzsebastian/iinuclear?branch=main)

## Quick Start

The simplest way to use `iinuclear` is:

```python
from iinuclear import isit

# Using an IAU name from TNS
isit("2018hyz")

# Using a ZTF name in Alerce
isit("ZTF18acpdvos")

# Using coordinates in degrees
isit(151.711964138, 1.69279894089)
```

This will create a diagnostic plot showing:
- PanSTARRS image of the host galaxy
- ZTF position measurements
- Statistical analysis of nuclear position
- Confidence assessment

<p align="center">
  <img src="docs/images/2018hyz_iinuclear.png" width="800">
  <br>
  <em>Example output for AT2018hyz showing: (top left) PanSTARRS image with ZTF positions overlaid, 
  (top right) separation histogram, (bottom left) ZTF position density plot, and 
  (bottom right) statistical significance analysis.</em>
</p>

## Installation

Install using pip:
```bash
pip install iinuclear
```

Or install from source:
```bash
git clone https://github.com/gmzsebastian/iinuclear.git
cd iinuclear
pip install -e .
```

## Requirements

* Python 3.7 or later
* Access to TNS API (for IAU name queries)
* Having the ``alerce`` API package installed
* Having the ``astroquery`` package installed

## Citation

If you use `iinuclear` in your research, please cite:

```bibtex
@software{iinuclear,
  author       = {Sebastian Gomez},
  title        = {iinuclear: Nuclear Transient Classifier},
  year         = {2024},
  publisher    = {GitHub},
  url          = {https://github.com/gmzsebastian/iinuclear}
}
```

## License

Copyright 2024 Sebastian Gomez and contributors.

`iinuclear` is free software made available under the MIT License. For details see the LICENSE file.
