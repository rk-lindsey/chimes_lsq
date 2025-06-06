.. chimes_lsq documentation master file, created by
   sphinx-quickstart on Mon Dec 14 14:22:56 2020.

.. Lines starting with two dots are comments if no special commands are found

ChIMES Parameter Generator Documentation
========================================

The **Ch**\ ebyshev **I**\ nteraction **M**\ odel for **E**\ fficient **S**\ imulation (ChIMES) is a machine-learned interatomic potential targeting chemistry in condensed phase systems. ChIMES models are able to approach quantum-accuracy through a systamtically improvable explicitly many-bodied basis comprised of linear combinations of Chebyshev polynomials. Though originally developed to enable description of organic molecular materials, ChIMES has successfuly been applied to systems spanning ambient water to molten carbon, and leveraged as correction for density functional based tight binding simulations. 

The ChIMES least squares (chimes_lsq) utility facilitates generates of ChIMES parameter sets from reference data, and can be used for particle-based simulations and single point calculations through using the `ChIMES Calculator <https://doi.org/10.1021/acs.jctc.7b00867>`_.

The ChIMES Calculator is developed at Lawrence Livermore National Laboratory with funding from the US Department of Energy (DOE), and is open source, distributed freely under the terms of the (specify) License.


For additional information, see:

.. toctree::
   :maxdepth: 1
   
   getting_started
   quick_start
   chimes_overview
   lsq_input_file
   units
   citing
   contributing
   legal



.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

---------------

This work was produced under the auspices of the U.S. Department of Energy by
Lawrence Livermore National Laboratory under Contract DE-AC52-07NA27344.

This work was prepared as an account of work sponsored by an agency of the
United States Government. Neither the United States Government nor Lawrence
Livermore National Security, LLC, nor any of their employees makes any warranty,
expressed or implied, or assumes any legal liability or responsibility for the
accuracy, completeness, or usefulness of any information, apparatus, product, or
process disclosed, or represents that its use would not infringe privately owned
rights. Reference herein to any specific commercial product, process, or service
by trade name, trademark, manufacturer, or otherwise does not necessarily
constitute or imply its endorsement, recommendation, or favoring by the United
States Government or Lawrence Livermore National Security, LLC. The views and
opinions of authors expressed herein do not necessarily state or reflect those
of the United States Government or Lawrence Livermore National Security, LLC,
and shall not be used for advertising or product endorsement purposes.
