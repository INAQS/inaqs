# GIFS
The "Generalized Interface to Surface Hopping" allows one to run gas-phase or solvated (QM/MM with electronic embedding) dynamics using GROMACS for the dynamics and Q-Chem for the electronic structure. We support Born Oppenheimer MD on ground or excited state surfaces and Fewest Switches Surface Hopping with optional decoherence.

For an overview of GIFS and its capabilites, see our [Q-Chem Webinar](https://www.q-chem.com/webinars/45/).

# Building
GIFS is relatively decoupled from Gromacs and is invoked via shared library.

To get started:
* Compile the shared library via make in gifs_src/gifs
* Build the bundled gromacs-4.6.5 being sure to define GMX_GIFS with the root of this repository; see the build.sh or setup.sh in the gromacs folder for details
* Make sure you're using Q-Chem 5.3.2 or later with you $QC variable set correctly
* Try some of the examples in the examples folder; run.sh will get each example going.

# Development
GIFS is developed as a collaboration between the University of Pennsylvania and the University of Groningen. Work in the Netherlands is led by Dr. Maximilian Menger (m.f.s.j.menger@rug.nl) and in the United States by [Dr. Vale Cofer-Shabica](https://vale.science) (valecs@sas.upenn.edu).

# Citation
If you use GIFS, please cite this repository; two papers are in preparation and we will update once they have been submitted.