.. highlight:: shell

============
Installation
============

* Compile the shared library via CMake in a separate directory
* Build the bundled gromacs-4.6.5 being sure to define GMX_GIFS with the root of this repository
* The inaqs_build.sh script will do both of the above if you like.
* Add the modified gromacs to your path:
::

   source .../path/to/inaqs/gromacs-4.6.5/install/bin/GMXRC.bash
* Make sure you're using Q-Chem 5.4.1 or later with you $QC variable set correctly
* Try some of the examples in the examples folder; run.sh will get each example going.


PLUMED
------

* If you want to use `PLUMED <www.plumed.org>` with INAQS, install PLUMED as you normally would.
* Run the prep-plumed.sh script in utils/plumed to write patch files for PLUMED.
* In the gromacs directory, run plumed patch --patch --engine gromacs-4.6.5
* Rebuild gromacs.
