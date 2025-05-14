.. _page-quick_start:

Quick Start Guide
=============================================

This page provides a step-by-step guide for generating a 2+3-body ChIMES model for a 1:1 C:O mixture at 9350 K and 2.5 g/cc.

.. note::  
    We *highly* recommend users read the :ref:`ChIMES Overview <page-overview>` and :ref:`Generating a ChIMES model <page-running>` before attempting to generate a simulation-ready model. 

---------------


Step 1: Obtain the tutorial files
*********************************************

Follow instructions on the :ref:`Getting Started <page-getting_started>` page to obtain a copy of the ChIMES LSQ repository and compile the code. Next, in a terminal, create a new folder called ``tutorial`` in your desired location, and copy the tutorial files there:

.. code-block:: bash

    cp /path/to/chimes_lsq/test_suite-lsq/special3b/fm_setup.in                   /my/tutorial_folder
    cp /path/to/chimes_lsq/test_suite-lsq/special3b/9350.combine-scs-2-a.fit.xyzf /my/tutorial_folder
    
    
Step 2: Inspect the input
*********************************************    

Aside from the ``chimes_lsq`` code, a minimum of two input files are needed to generate a ChIMES model: A training trajectory and a ``chimes_lsq`` configuration file (``9350.combine-scs-2-a.fit.xyzf`` and ``fm_setup.in`` in this example, respectively). By running ``head 9350*`` in the tutorial folder, one can see that the trajectory file is in an extended xyz format, where for a given frame, the first line provides the number of atoms, the next line, the three box lengths (Å), and the remaining lines provide an atom type, *x*, *y*, and *z* coordinates (Å), and the corresponding *x*, *y*, and *z* components of the forces acting on those atoms, in Hartree (H) per Bohr (B). Frame information is repeated in this format :math:`n_{\mathrm{frames}}` times. 

.. note::

    In this example, we're only concerned with generating a ChIMES model based on forces, but in practice, one can also fit to energies and/or stresses as well. Similarly, in this example we're concerned with an orthorhombic system, but triclinic cells can be used as well. For additional details, see the :ref:`Input Files and Options <sec-input_files>` section of :ref:`Generating a ChIMES model <page-running>`.
    
Next, open the configuration file, ``fm_setup.in``, in your text editor of choice. Add the lines that are highlighted below.

The "control variables" section is as given below. Note that ``# VARNAME #`` are ``chimes_lsq`` key words and should not be edited; key word settings are specified on the following line.

.. tip::

    ChIMES configuration file comments can either be made as new lines prepended with a hash ("#"), or at the end of a line, prepended with an exclamation mark (!).

.. code-block:: bash
    :emphasize-lines: 19, 20
    
    ####### CONTROL VARIABLES #######

    # TRJFILE #
            9350.combine-scs-2-a.fit.xyzf   ! Training trajectory file
    # NFRAMES # 
            100                             ! Number of frames in training trajectory file
    # NLAYERS #
            1                               ! Number of replicate layers
    # FITCOUL #
            false                           ! Fit charges?
    # FITSTRS # 
            false                           ! Fit stresses? If so, how?
    # FITENER #
            false                           ! Fit energies? If so, how?
    # PAIRTYP # 
            CHEBYSHEV 12 5 0                ! Polynomial orders - See note below!!!
    # CHBTYPE #
            MORSE                           ! Coordinate transformation style
    # SPLITFI #
            false                           ! Should the design matrix be split into multiple files?


Next comes the "topology variables", which specify atom-type and interaction-cluster-type-specific aspects of the model:

.. code-block:: bash
 
    ####### TOPOLOGY VARIABLES #######

    # Specify the total number of atom types being fit

    # NATMTYP # 
            2

    # For each atom type, specify the type index, atom type, fixed charge, and mass.

    # TYPEIDX #     # ATM_TYP #     # ATMCHRG #     # ATMMASS # 
    1               C                0.0            12.011 
    2               O                0.0            15.999

    # Next, figure out how many pair types can be generated. Specify the pair type index, constituent atom types,
    # inner cutoff (s_minim), outer cutoff (s_maxim), and morse transformation variable. Note that s_delta is a 
    # vestigial field that can use any floating point number as a place holder.

    # PAIRIDX #     # ATM_TY1 #     # ATM_TY2 #     # S_MINIM #     # S_MAXIM #     # S_DELTA #     # MORSE_LAMBDA #
    1               C               C               0.988           8.00            0.01            1.35            
    2               O               O               0.911           8.00            0.01            2.55            
    3               C               O               0.876           8.00            0.01            1.15            


Model efficiency can be increased by using many-body cutoffs that are shorter than those used for two-body interactions (specified above). The ``fm_setup.in`` section below shows how to do so for 3-body interactions. Here, we want to do so for each of the 4 possible 3-body interaction types (specified by ``SPECIAL 3B S_MAXIM: SPECIFIC 4``). The following line shows how this is done for a 3-body interaction comprised entirely of carbon atoms. For 3 carbon atoms, C, we have 3 atom pairs, CC, CC, CC, which yields a cluster "name" of CCCCCC. For each of these constituent atom pairs, we want to use 3-body outer cutoffs of 4.4 Å. For a C O O cluster, we have pairs CO, CO, OO, for which we would like to use respective 3-body outer cutoffs of 4.0, 4.0, and 6.5.

.. code-block:: bash
 
    SPECIAL 3B S_MAXIM: SPECIFIC 4
    CCCCCC CC CC CC 4.4 4.4 4.4
    COCOCC CC CO CO 4.4 4.0 4.0
    OOCOCO CO CO OO 4.0 4.0 6.5
    OOOOOO OO OO OO 6.5 6.5 6.5

Finally, we specify that we will use a cubic smoothing function, and then specify the end of the configuration file:

.. code-block:: bash

    # FCUTTYP #
            CUBIC   

    # ENDFILE #
    
.. Warning::

    The polynomial orders listed in the input file for ``# PAIRTYP # CHEBYSHEV`` should be read as <O2B> <O3B+1> <O4B+1>, when O3B and O4B are greater than zero. For example, ``CHEBYSHEV 12 5 0`` actually indicates a model that is 12th order in 2-body interactions, and 4th order in 3-body interactions. 
    
.. Tip::

    Note: All of the following information is provided in published ChIMES manuscripts.
    
    * Inner cutoffs are set to the lowest sampled distance (:math:`r_{\mathrm{samp,min}}`) in for each given pair type the training trajectory, or slightly less (:math:`r_{\mathrm{samp,min}}-0.02`).

    * 2-body outer cutoffs are usually set to encompass at least the 2nd non-bonded solvation shell, and usually set to about 8 Å.

    * 3-body outer cutoffs are usually set to encompass the 1st non-bonded solvation shell.

    * 4-body outer cutoffs are usually set to between the first or second RDF minimum; choice depends on the nature of the system (i.e. do extended 4-body interactions like molecular torsions need to be described, or is the system characterized by short-ranged interactions like in metallic clusters? 

    * The Morse lambda is usually set to distance corresponding to the the first (bonding) peak in the in the pair RDF.

    * There is  no hard-and-fast rule for setting the 2/3/4-body polynomial orders, but 12/7/3 (i.e., ``CHEBYSHEV 12 8 4``) is usually a good starting point.

    * Weighting is often needed when energies and stresses tensors are included in the fit. Typical respective weighting factors are 0.1 – 5.0 and 200 – 500. This will be discussed in greater detail in the following section.

    * If a TERSOFF smoothing function is used (i.e., ``# FCUTTYP # TERSOFF <var>``) The TERSOFF var determines when the smooth step cutoff function kicks in and is typically set to between 0.25 and 0.50, where larger values modify the interaction at shorter distances, but can make it easier to obtain a smooth interaction.
    
    * Be mindful of outer cutoff distances with respect to box lengths. 

    * If you don’t want to determine all pair/triplet/quadruplet interaction types and inner cutoffs by hand, setup a dummy fm_setup.in file with: 
    
        * ``CHEBYSHEV 2 2 2``
        * 2-body inner cutoffs = 0.0
        * 2-body outer cutoffs = 0.5
        * Do not explicitly set 3- or 4-body outer cutoffs (the code will use the 2-body value for them)
        * Run the chimes_lsq code using this file
        * Check the bottom section of the output, which lists all this information for you
    
Step 3: Generate the design matrix
*********************************************   

As described in :ref:`ChIMES Overview <page-overview>` and :ref:`Generating a ChIMES model <page-running>`, ChIMES is parameterically linear, meaning the fitting problem can be recast as a matrix equation of the form :math:`\mathbf{Ax=b}`. The purpose of ``chimes_lsq`` is to generate :math:`\mathbf{A}` and :math:`\mathbf{b}` based on the user-defined ChIMES hyperparameters (i.e., contents of ``fm_setup.in``), and training trajectory (e.g., ``9350.combine-scs-2-a.fit.xyzf``). To do this, simply run the following from your tutorial folder:

.. code-block::
    
    /path/to/chimes_lsq/build/chimes_lsq fm_setup.in | tee fm_setup.log
    
Running ``ls`` will show generation of ``A.txt`` and ``b.txt`` files, which we will use to determine our model parameters (i.e. :math:`\mathbf{x}` in the equation :math:`\mathbf{Ax=b}`). Several other files are produced including:

* dim.txt: A-matrix dimension <# columns> <# rows>

* natoms.txt: For each line in b.txt, provides the number of atoms that were in the corresponding frame

* b-labeled.txt: Indicates the atom type for each force component, and if applicable, labels stress tensor components (prepended with "s\_"), and system energies (labeled "+1").

* params.header: Formatted output containing ChIMES hyperparameters

* ff_groups.map: Assigns cluster type indices for cluster name permutations


Step 4: Solve for ChIMES parameters
*********************************************   

We will use principal component analysis (PCA) based on the singular value decomposition (SVD) of the force derivative matrix with the (default) regularization of 1E-5 to solve for ChIMES parameters:

.. code-block::

    python /path/to/chimes_lsq/build/chimes_lsq.py > params.txt
    
Note, this step can take a few minutes to run.

Once complete, the script produces two files, ``params.txt``, the ChIMES parameter file, and ``force.txt``, the forces (and optionally stresses and energies) predicted for the training trajectory based on the presently developed ChIMES model. To inspect the RMS error in predicted forces (kcal/mol/Å), run:

.. code-block::
    
    grep -F "RMS force error" params.txt
    
and to visualize model peformance, run:

.. code-block::

    paste b.txt force.txt > compare.txt
    xmgrace compare.txt
    
    
Step 5: Prepare for MD simulation
*********************************************       
    
Before using resulting parameters in MD simulations, be sure to add penalty function parameters. To do so, locate the containing ``FCUT TYPE: CUBIC``, and enter add the highlighted text:

.. code-block::
    :emphasize-lines: 8, 9

    # PAIRIDX #     # ATM_TY1 #     # ATM_TY1 #     # S_MINIM #     # S_MAXIM #     # CHBDIST #     # MORSE_LAMBDA #
            0               C               C               0.988           8               MORSE           1.35
            1               O               O               0.911           8               MORSE           2.55
            2               C               O               0.876           8               MORSE           1.15
    
    FCUT TYPE: CUBIC
    
    PAIR CHEBYSHEV PENALTY DIST:    0.02
    PAIR CHEBYSHEV PENALTY SCALING: 1E6
    
    SPECIAL 3B S_MAXIM: SPECIFIC 4
    CCCCCC CC CC CC 4.40000 4.40000 4.40000
    CCCOCO CC CO CO 4.40000 4.00000 4.00000
    COCOOO CO CO OO 4.00000 4.00000 6.50000
    OOOOOO OO OO OO 6.50000 6.50000 6.50000

Simulations using the resulting model can now be run with ``chimes_md`` (ancillary support), or `the ChIMES Calculator <https://github.com/rk-lindsey/chimes_calculator>`_

.. Warning::

    This tutorial is only intended to familiarize the user with ``chimes_lsq``. Generating an accurate and stable machine learned model requires careful consideration of the target application, relevant training data, and suitable model hyperparameters. The :ref:`ChIMES Overview <page-overview>` and :ref:`Generating a ChIMES model <page-running>` pages provide a more in-depth discussion of some of these topics and their relation to ``chimes_lsq``. We also recommend reviewing the :ref:`ChIMES literature <page-citing>` to see how models for complex problems are developed. 