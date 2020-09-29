# Coupled Phase-Field Brittle Fracture Simulation

A phase-field brittle fracture mechanics code written in [Nutils](http://www.nutils.org/en/latest/) is coupled to a dummy Python script providing material parameters (fracture toughness, first and second Lamé parameters). The coupling is realized with [preCICE](https://www.precice.org/). The purpose of this code is to showcase how preCICE could be used to volume-couple between an electro-chemistry model (here the dummy script) and a fracture mechanics model.

The original fracture mechanics code was developed by [Clemens Verhoosel](https://www.tue.nl/en/research/researchers/clemens-verhoosel/). The model is documented in this publication:
 
Singh, N., Verhoosel, C.V., De Borst, R. and Van Brummelen, E.H., 2016. A fracture-controlled path-following technique for phase-field modeling of brittle fracture. _Finite Elements in Analysis and Design_, **113**, pp.14-29.

The coupling prototype was developed by [Benjamin Uekermann](https://github.com/uekerman).
Contact Benjamin for questions about this code. Furthermore, the following support channels might be helpful: [Nutils matrix chat room](https://matrix.to/#/#nutils-users:matrix.org), [preCICE forum](https://precice.discourse.group/).

## Dependencies

* Python3 (e.g. 3.6, but any new version should work)
* numpy (e.g. 1.18, but any new version should work)
* matplotlib (e.g. 3.0, but any new version should work)
* [Nutils v6](https://pypi.org/project/nutils/) 
* [treelog](https://github.com/evalf/treelog) (should already be fetched as dependency of nutils)
* [preCICE v2](https://github.com/precice/precice/releases/tag/v2.1.0)
* [pyprecice v2](https://pypi.org/project/pyprecice/) 

## Folder structure

├── `Allclean`: a bash script to clean the case 
├── `dummy.py`: the dummy code, which provides the material parameters
├── `fracture.py`: the fracture mechanics, implemented with Nutils
├── `precice-config.xml`: the preCICE configuration file, describing the coupled setup
└── `README.md`

## How to run

Open two terminals and start one program in each directly from the root directory:
* `python3 fracture.py`
* `python3 dummy.py`

## Visualization

Nutils collects and outputs data in html. Simply open the link that Nutils outputs in your browser.
Something like: `file://///home/uekerman/public_html/fracture.py/log.html`
Click on links there and take notice of the menu on the top right.
