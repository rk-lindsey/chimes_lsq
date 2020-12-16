Getting Started
=============================================

* :ref:`Obtaining the code     <sec-obtaining>`
* :ref:`Compiling and running  <sec-compiling>`

---------------

.. _sec-obtaining:

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

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_lsq/browse>`_ page, click the fork button (fourth button down from the top left, which looks like a two-prong trident). This will create a copy of ``chimes_lsq`` in your LLNL Bitbucket account. You can clone, edit, and commit to this repository as you see fit, and any contribibutions to the original "parent" reposoitory via the pull request mechanism. For additional details, see :ref:`Contributing <page-contributing>`.


2. Cloning
^^^^^^^^^^

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_lsq/browse>`_ page, click the clone button (first button from the top left, which looks like a monitor with an arrow). Copy the ``chimes_lsq`` repository url that appears. In a terminal, navigate to the desired download location, an type: ``git clone <the full copied url> .``


2. Archive Download
^^^^^^^^^^^^^^^^^^^

On the `LLNL-hosted Bitbucket repository <https://mybitbucket.llnl.gov/projects/CHMS/repos/chimes_lsq/browse>`_ page, click the three dots next to "master" under the heading "Source," and select download from the drop down menu


---------------


.. _sec-compiling:

Compiling and running the code
****************************************

To compile, simply navigate to the ``src`` folder and type ``make chimes_lsq``. Parameter set generation requires two steps, i.e. generation and solution of the design matrix. Provided suitable input (e.g. fm_setup.in, a ``chimes_lsq`` input file and a reference data file), this can be as simple as:

.. code-block:: bash
    
    /path/to/repo/src/chimes_lsq fm_setup.in > fm_setup.log
    /path/to/repo/src/chimes_lsq.py > params.txt
    
    
Note: As described in greater detail in :ref:`Generating a ChIMES model <page-running>`, chimes_lsq.py depends on native `numpy <https://numpy.org>`_, `scipy <https://www.scipy.org>`_, and `sklearn <https://scikit-learn.org/stable/>`_ installations for `python2.x <https://www.python.org>`_ **[THIS NEEDS TO BE UPGRADED TO 3.X]**.

For a more detailed description of how to use ``chimes_lsq``, see: :ref:`Generating a ChIMES model <page-running>`.
