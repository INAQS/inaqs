=================================================
Generalized Interface for Surface Hopping (GIFS)
=================================================

Idea: Provide a simple api to perform Surface Hopping (SH)
      simulations using an molecular dynamics code (MD) of choice
      e.g. Gromacs. 

Decouple the surface hopping algorithm (active state, electronic propagation, 
qm calculations e.g.) from the nuclear propagation scheme and in case
of Gromacs as well from the pure MM calculations, if QM/MM is used.

To do so, an abstract base class is definded (**Gifs**), which acts as 
a factory to create an concrete **SurfaceHopping** class, containing
the actual implementation.

From the MD side, as only Born-Oppenheimer dynamics is assumed, only
a few interface method need to be exposed, namely:

 * get the energy/gradient at a given nuclear position
 * perform the surface hopping related computations
   (e.g. electronic propagation etc.)
 * rescale velocities after hop

**@Vale**: Can you ask Joe, if this minimalistic interface is enough, or
           if more concrete functions of `Gifs` need to be exposed? I
           personaly think it is enough for our purposes.

.. literalinclude:: ../gifs_src/gifs/gifs/gifs.hpp
   :language: cpp
   :lines: 11-35
   :emphasize-lines: 2-24


Surface Hopping
---------------


The details of the actual Surface Hopping implementation are hidden
in the private instance of the abstract SurfaceHopping class.
The latter is used to give an interface to concrete surface hopping
implementations (e.g. Tullys fewest switches, Landau-Zener, ...), which
can be in that way used without changing the public interface, as long as 
the factory method knows the concrete type.

.. literalinclude:: ../gifs_src/gifs/sh/sh.hpp
   :language: cpp
   :lines: 7-37
   :emphasize-lines: 2-29

The surface hopping routine is responsible to expose the actual implementations
of the get_energy, get_eandg member functions, needed by the **Gifs** class.
As well as, a way to perform qm(/mm) calculations via an **ElectronicStructure**
method.

There are basically two different kinds of electronic structure calculations one, wants
to perform:

  * Mechanical Embedding/No-embedding 
    (MM not present in QM calculation)
  * Electrostatic Embedding           
    (MM included in terms of fixed point charges)

Gromacs supports both cases, and it should make it more extensible, in case someone
wants to support different flavours like polarizable embedding.


QM Calculations
---------------

All electronic structure calculations are performed 
via a concrete implementation of the abstract base class **ElectronicStructure**.

Electronic Structure
~~~~~~~~~~~~~~~~~~~~

The class defines exposed the concrete method to perform the electronic structure 
computations using the **compute** method. 
In the simplest implementation it should just perform energy, gradient computations.
In that way, one can in principle use the framework to perform normal QM/MM calculations
with Gromacs without computational overhead, by providing a fake **SurfaceHopping** class 
that only performs the QM/MM calculations.

**TODO:** Is there a better way to call the electronic structure in a general way, then
          to use polymorphic input/output types?

.. literalinclude:: ../gifs_src/gifs/sh/es.hpp
   :language: cpp
   :lines: 68-86
   :emphasize-lines: 2-17

Input
~~~~~

The input class contains all necessary information to
perform the acutal QM(/MM) computation.
In its most simplest representation it just contains the
coordinates of the system for which the properties should
be computed.


**TODO:** Which properties are always needed, which should be optional?
          

.. literalinclude:: ../gifs_src/gifs/sh/es.hpp
   :language: cpp
   :lines: 11-28
   :emphasize-lines: 9-17

Output
~~~~~~

The output class contains all properties that
are computed by the electronic structure method.

**TODO:** Which properties are always needed, which should be optional?

.. literalinclude:: ../gifs_src/gifs/sh/es.hpp
   :language: cpp
   :lines: 44-57
   :emphasize-lines: 4-13

