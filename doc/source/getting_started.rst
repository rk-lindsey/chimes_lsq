.. _page-getting_started:

Getting Started
=============================================

* :ref:`Obtaining the code     <sec-obtaining>`
* :ref:`Compiling and running  <sec-compiling>`

---------------

.. _sec-obtaining:


Obtaining the code (Non-LLNL Users):
****************************************

The ``chimes_lsq`` code is stored in two separate GitHub repositories, one maintined by the Lindsey Lab at UM, and the other maintained by collaborators at LLNL. The repository links are given below:

* `UM-maintained   <https://github.com/LindseyLab-umich/chimes_lsq-LLfork>`_
* `LLNL-maintained <https://github.com/rk-lindsey/chimes_lsq>`_

To use these codes, we recommend the following procedure:

1. Fork the codes 
^^^^^^^^^^^^^^^^^

For each repository (e.g., chimes_lsq-myLLfork):

* Click the repository name
* At the top, next to the repository name, click the fork button
* In the window that appears, select your github user name as the owner, and add "-fork" to the repository name
* Click fork

2. Cloning the codes
^^^^^^^^^^^^^^^^^^^^

* Execute mkdir -p ~/Codes in a terminal window
* For each of your newly forked repositories, (e.g., chimes_lsq-myLLfork):
  * Go to the github website (e.g., https://github.com/RKLindsey/chimes_lsq-myLLfork)
  
  * Click the green "Code" button
  * Copy the address starting with "git@github.com:"
  * Execute cd ~/Codes in a terminal window
  * Execute git clone <copied git@github.com:> .
  * Don't forget the "." at the end!
  
---------------

Obtaining the code (For LLNL Employees):
****************************************

Note that ``chimes_lsq`` is stored in a `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_lsq/browse>`_. The following steps must be completed to obtain/access ``chimes_lsq``:

* Send an e-mail to lindsey11 titled "ChIMES LSQ BB Acess Request" to be granted access. 
* If the code is being copied to the LC: Set up SSH RSA keys

    * `Generate your keys <https://www.ssh.com/ssh/keygen/>`_. Use the default location and do not leave the pass-phrase blank.
    * Copy the contents of your key to the LLNL-hosted Bitbucket and add them to your account `(accessible here) <https://mybitbucket.llnl.gov/plugins/servlet/ssh/account/keys>`_

Once access is granted and, if applicable, keys have been configured, ``chimes_lsq`` it can be obtained three ways, which are described below:

1. Forking (recommended for developers/contributors)
2. Cloning 
3. Downloading as an archive

1. Forking
^^^^^^^^^^

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_lsq/browse>`_ page, click the fork button (fourth button down from the top left, which looks like a two-prong trident). This will create a copy of ``chimes_lsq`` in your LLNL Bitbucket account. You can clone, edit, and commit to this repository as you see fit, and any contributions to the original "parent" reposoitory via the pull request mechanism. For additional details, see :ref:`Contributing <page-contributing>`.


2. Cloning
^^^^^^^^^^

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_lsq/browse>`_ page, click the clone button (first button from the top left, which looks like a monitor with an arrow). Copy the ``chimes_lsq`` repository url that appears. In a terminal, navigate to the desired download location, and type: ``git clone <the full copied url> .``

.. note::

    To clone on a LLNL LC machine, execute:
    
    ``git clone ssh://git@mybitbucket.llnl.gov:7999/chms/<repo name> .``
    
    From any other computer, use:
    
    ``git clone https://mybitbucket.llnl.gov/scm//chms/<repo name> .``
    
    Allow 30 minutes to an hour for the clone to complete.
    

3. Archive Download
^^^^^^^^^^^^^^^^^^^

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_lsq/browse>`_ page, click the three dots next to "master" under the heading "Source," and select download from the drop down menu.  
  
  
  
  
  
  
---------------



Compiling and running the code
****************************************

Navigate to the root project directory (e.g, chimes_lsq-myLLfork), and then execute ``ls modfiles``. You will see a list of files named like UM-ARC.mod. The name of these files correspond to different high-performance computing systems (HPC) and contain the modules necessary to compile the code on those platforms. e.g., UM-ARC is for the UM ARC Great Lakes HPC. If you are on one of these computers, execute, e.g.,  ``export hosttype=UM-ARC`` prior to installation. To complete installation, execute ``./install.sh`` in the root project directory. 

.. note::

      If you do not see your platform represented in the modfiles, we suggest you create one that contains modules necessary for compiling C, C++, Fortran, and MPI code. You may need to include a module for cmake, python. You will also need to update the install.sh script accordingly. See, e.g., UM-ARC.mod for reference. As described in greater detail in :ref:`Generating a ChIMES model <page-running>`, chimes_lsq.py depends on native `numpy <https://numpy.org>`_, `scipy <https://www.scipy.org>`_, and `sklearn <https://scikit-learn.org/stable/>`_ installations for `python3.x <https://www.python.org>`_. 

Parameter set generation requires two steps, i.e. generation and solution of the design matrix. Provided suitable input (e.g. fm_setup.in, a ``chimes_lsq`` input file and a reference data file), this can be as simple as:

.. code-block:: bash
    
    /path/to/repo/build/chimes_lsq fm_setup.in > fm_setup.log
    python3 /path/to/repo/build/chimes_lsq.py > params.txt
    
For a more detailed description of how to use ``chimes_lsq``, see: :ref:`Generating a ChIMES model <page-running>`.



