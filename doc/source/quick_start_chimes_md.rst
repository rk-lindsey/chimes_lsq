.. _page-quick_start:

Quick Start Guide
=============================================
 
This provides a step-by-step guide for generating a 2+3-body ChIMES MD simulation for H2O at 300 K and 1.0 g/cc. 

.. note::   

    We *highly* recommend users read the :ref:`ChIMES Overview <page-overview>` and :ref:`Generating a ChIMES model <page-running>` before attempting to run a simulation using ChIMES MD software.

Step 1: Obtain the tutorial files
*********************************************

Follow instructions on the :ref:`Getting Started <page-getting_started>` page to obtain a copy of the ChIMES LSQ repository and compile the code. Next, in a terminal, create a new folder called ``tutorial`` in your desired location, and copy the tutorial files there:

.. code-block:: bash

    cp /path/to/chimes_lsq/test_suite-lsq/test_suite-md/h2o-2bcheby/params.cheby.txt          /my/tutorial_folder
    cp /path/to/chimes_lsq/test_suite-lsq/test_suite-md/h2o-2bcheby/run_md.in          /my/tutorial_folder
    cp /path/to/chimes_lsq/test_suite-lsq/test_suite-md/h2o-2bcheby/input.xyz          /my/tutorial_folder
    
    
Step 2: Inspect the input
*********************************************    

Aside from the ``chimes_md`` code, a minimum of three input files are needed to run the MD code: A trajectory file, parameters file and a ``chimes_md`` configuration file (``input.xyz``, ``params.cheby.txt`` and ``run_md.in`` in this example, respectively). By running ``head input.xyz`` in the tutorial folder, one can see that the trajectory file is in an extended xyz format, where for a given frame, the first line provides the number of atoms, the next line, the three box lengths (Å), and the remaining lines provide an atom type, *x*, *y*, and *z* coordinates (Å), and the corresponding *x*, *y*, and *z* components of the forces acting on those atoms, in Hartree (H) per Bohr (B). Frame information is repeated in this format :math:`n_{\mathrm{frames}}` times. 

.. note::

    In this example, we're only concerned with generating a ChIMES MD files. For additional details, see the :ref:`Input Files and Options <sec-input_files>` section of :ref:`Generating a ChIMES model <page-running>`.
    
Next, open the configuration file, ``run_md.in``, in your text editor of choice. Add the lines that are highlighted below.


.. tip::

    ChIMES configuration file comments can either be made as new lines prepended with a hash ("#"), or at the end of a line, prepended with an exclamation mark (!).

.. code-block:: bash
   :emphasize-lines: 28
    
    ## Notes: Compare with "cheby_md.in // params.cheby.txt in non-generalized version of the code's
    ##        h2o_md example folder.
 
    ###################################
    #### GENERAL CONTROL VARIABLES ####
    ###################################
 
 
    # RNDSEED # ! Seed. If not specified, default value 123457 is used
    12357
    # TEMPERA # ! In K
    2000.0
    # CMPRFRC # ! Compare computed forces against a set of input forces? ...If true, provide name of the file containing the forces for comparison
    false
    # TIMESTP # ! In fs
       0.025
    # N_MDSTP # ! Total number of MD steps
     16000
    # NLAYERS # ! x,y, and z supercells.. small unit cell should have >= 1
      2
    # USENEIG # 
     true 
    # PRMFILE # ! Parameter file (i.e. params.txt)
      params.cheby.txt
    # CRDFILE # ! Coordinate file (.xyz) or force file (.xyzf)
      input.xyz 
    # TRAJEXT # ! coordinate file type 
      XYZ
 
    ###################################
    ####    SIMULATION  OPTIONS    ####
    ###################################
 
    # VELINIT # (options are READ or GEN)
      READ
    # CONSRNT # (options are HOOVER <hoover time> or VELSCALE <scale freq>)
      NVT-MTK HOOVER 50
    !# CONSRNT # (options are HOOVER <hoover time> or VELSCALE <scale freq>)
    !      NVE
    # PRSCALC # (options are ANALYTICAL or NUMERICAL)
      ANALYTICAL
    ###################################
    ####      OUTPUT  CONTROL      ####
    ################################### 
     # WRPCRDS # 
      false 
    # FRQDFTB # ! Frequency to output the DFTB gen file can also be called FRQTRAJ
      20
    # FRQENER # ! Frequency to output energies
      10 
    # PRNTFRC # ! Print computed forces? Forces are printed to force_out.txt 
      false

    # ENDFILE #
      
.. Note::

    * The box dimensions are read from the comment line of the xyz.
    * The box dimsnsion can be int two different format; Non-orthogonal which is denpted as ``NON_ORTHO <xx> <xy> <xz> <yx> <yy> <yz> <zx> <zy> <zz>``, orthogonal which is ``<x> <y> <z>``. The default is orthogonal.
    
Step 3: Run the MD simulation
*********************************************   

ChIMES MD can be ran parallel or serial. To run the MD simulation in serial, run the following command:
.. code-block::
    
    /path/to/chimes_lsq/build/chimes_md-serial run_md.in > run_md.log

To run the MD simulation in parallel, run the following command:
.. code-block::
    
    /path/to/chimes_lsq/build/chimes_md run_md.in > run_md.log

This will generate run_md.out, traj_bad_r.lt.rin.xyz, traj_bad_r.lt.rin+dp.xyz, traj_bad_r.ge.rin+dp_dftbfrq.xyz, traj.xyz, run_md.out, restart.xyzv ,restart.bak ,output.xyz , and md_statistics.out, but the most important files are ``traj.xyz`` , ``run_md.out`` and ``md_statistics.out``.
``traj.xyz`` is printed based on  ``FRQDFTB`` aka ``FRQTRAJ`` which is the trajectory output frequency, ``run_md.out`` is the main output of MD simulation, and ``md_statistics.out`` is the statistics of the MD simulation printed every ``FRQENER``.

.. tip::

    To check the conserved quantity, you can plot the 2nd vs 8th column of ``md_statistics.out``.
