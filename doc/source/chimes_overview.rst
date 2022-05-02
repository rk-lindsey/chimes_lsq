.. _page-overview:

ChIMES Overview
=============================================

This page provides an overview of key ChIMES equations. For further details, see the complete set of ChIMES references available at :ref:`Citing ChIMES <page-citing>`.


The ChIMES equations
********

ChIMES energy
^^^^^^^^^^^^^



ChIMES is a reactive explicitly many-bodied machine learned interatomic potential (ML-IAP) for which interactions are 
computed on the basis of atom clusters. For example, the total ChIMES energy is given as:

.. math::
   :nowrap:

   \begin{eqnarray}
      E_{n_\mathrm{B}} 
      = \sum_{i_1}        ^{n_a} {}^{1}\!E_{i_1}
      + \sum_{i_1>i_2}    ^{n_a} {}^{2}\!E_{i_1i_2} 
      + \sum_{i_1>i_2>i_3}^{n_a} {}^{3}\!E_{i_1i_2i_3} 
      + \dots
      + \sum_{i_1>i_2\dots >i_{n\mathrm{B}-1}>i_{n\mathrm{B}}}^{n_a} {}^{n}\!E_{i_1i_2\dots i_n}
   
   \end{eqnarray}

   
where :math:`E_{n_\mathrm{B}}` is the total ChIMES system energy, :math:`n_{\mathrm{B}}` is the maximum bodiedness, 
:math:`{}^{n}\!E_{i_1i_2\dots i_n}` is the :math:`n`-body ChIMES energy for  a given set of :math:`n` atoms with indices :math:`i = {i_1, i_2, \dots , i_n}`, and :math:`n_a` is the total number of atoms in the system. In the ChIMES framework, single-body energies are constant values and :math:`n`-body energies are constructed from the product of polynomials of transformed atom pair distances. Thus, a 2-body interaction would involve a single pair distance, :math:`ij` and corresponding polynomial :math:`T_n(s_{ij})`, where :math:`s_{ij}` is the pair distance transformed to a range of :math:`[-1:1]`. Similarly, a  three-body interaction would involve :math:`3\choose 2` (i.e. three) pairs, :math:`ij, ik,` and :math:`jk`, and a corresponding three-body polynomial, :math:`T_\alpha(s_{ij}) T_\beta(s_{ik}) T_\gamma(s_{jk})`. A 4-body interaction would involve :math:`4\choose 2 = 6` pairs, 
and so on. Currently, the ChIMES parameter generator supports up to 4-body interactions.


ChIMES atom cluster energies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. Taking a 3-body interaction as an example, we define the following: :math:`\mathbf{A}=\{i,j,k\}` is the index over atoms within an atom interaction cluster, with the corresponding set of pairs given by :math:`\mathbf{P}=\{ij,ik,jk\}`, their element pair types by :math:`\mathbf{E}=\{e_ie_j, e_ie_k, e_je_k\}`, and the polynomial order for each pair given by :math:`\mathbf{O}=\{\alpha , \beta , \gamma\}`. 

.. Two mapping functions are used to relate pair indices :math:`\mathbf{P}` to the three aforementioned pair properties: :math:`m_1 = \mathbf{P}\to \mathbf{E}` and :math:`m_2 = \mathbf{P}\to \mathbf{O}`, where an index :math:`y` refers to a particular component of :math:`\mathbf{P}`, defining an interaction pair.

.. Using these definitions, the generalized ChIMES energy for a cluster of :math:`n` atoms is written:

.. .. math::
   :nowrap:
   
   \begin{equation}
   ^{n}\!E = \prod _{y \in \mathbf{P}} f_s^{m_1(y)} (r_y)    \times \sum_{\mathbf{O}}^{\mathcal{O}^*}c_{\mathbf{O}}^{\mathbf{E}} \prod _{y \in \mathbf{P}} T_{m_2(y)}(s_y^{m_1(y)}),
   
   \end{equation}
   
.. where the :math:`\sum_{\mathbf{O}}` notation indicates a multiple sum for which there are :math:`\binom{n}{2}` distinct indices, :math:`\mathcal{O}^*` is the maximum polynomial order for a :math:`n`-body interaction, and the asterisk indicates a sufficient number of non-zero terms exist such that the graph formed by the edges of interacting atoms connects all :math:`n` atoms, which guarantees a true :math:`n`-body interaction. :math:`T_{m_2(y)}(s_y^{m_1(y)})` is a Chebyshev polynomial of order :math:`m_2(y)` that depends on pair distance :math:`s_y^{m_1(y)}` for pair :math:`y` of atom types :math:`m_1(y)` that has been transformed from :math:`r_y` to ensure it existis in the [-1,1] domain over which Chebyshev polynomials are defined, and :math:`f_s^{m_1(y)} (r_y)` is a cutoff function that ensures smooth behavior at the outer cutoff.


Using the above definition of :math:`^{n}\!E`, one arrives at the following for a two-body interaction:

.. math::
   :nowrap:
   
   \begin{equation}
   ^{2}\!E = f_s^{e_ie_j} (r_{ij}) \sum_{\alpha = 1} ^{\mathcal{O}_{2B}} c_{\alpha}^{e_ie_j} T_{\alpha}(s_{ij}^{e_ie_j}).
   
   \end{equation}
   
Here, :math:`\mathcal{O}_{2B}` is the maximum two-body polynomial order, :math:`c_{\alpha}^{e_ie_j}` is the permutationally invariant linear coefficient of polynomial of order :math:`\alpha`, :math:`e_ie_j` are the element pairs of atoms :math:`i` and :math:`j`, and :math:`T_{\alpha}(s_{ij}^{e_ie_j})` is the Chebyshev polynomial of the first kind computed for the transformed distance :math:`s_{ij}^{e_ie_j}`. The function :math:`f_s^{e_ie_j}(r_{ij})` is a smooth cutoff function, discussed below (:ref:`smooth cutoff function  <sec-cutoff>`). Note that an additional penalty term, :math:`f_p^{e_ie_j}` is added to each two-body interaction, described in greater detail below.

For a higher-body interaction, e.g. between three atoms, one obtains the following:

.. math::
   :nowrap:
   
   \begin{equation}
   ^{3}\!E  = f_s^{e_ie_j}(r_{ij})f_s^{e_ie_k}(r_{ik})f_s^{e_je_k}(r_{jk})  \sum_{\alpha = 1} ^{\mathcal{O}_{3B}} \sum_{\beta  = 1} ^{\mathcal{O}_{3B}} \sum_{\gamma = 1} ^{\mathcal{O}_{3B}}   c_{\alpha , \beta , \gamma}^{\mathbf{E}} T_{\alpha}(s_{ij}^{e_ie_j})  T_{\beta }(s_{ik}^{e_ie_k}) T_{\gamma}(s_{jk}^{e_je_k}).
   
   \end{equation}

The ChIMES penalty function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Two-body ChIMES interactions include an additional penalty function, which prevents spurious close contacts and discourages sampling of interatomic distances smaller than model inner cutoffs (i.e. where the Chebyshev polynomials are undefined). The penalty function takes on the following form:

.. math::
   :nowrap:
   
   \begin{equation}
   f_p^{m_1(y)} (r_y) =  
   \begin{cases}
   A_{\mathrm{p}}^{m_1(y)}(r_{\mathrm{c,in}}^{m_1(y)}+d^{m_1(y)}_{\mathrm{p}} -r_y)^3,& \text{if } r_y < r_{\mathrm{c,in}}^{m_1(y)}+d^{m_1(y)}_{\mathrm{p}}\\
   0,              & \text{otherwise}
   \end{cases}
   \end{equation}  

where :math:`r_{\mathrm{c,in}}`, :math:`A_{\mathrm{p}}`, and :math:`d^{m_1(y)}_{\mathrm{p}}` are inner cutoffs, penalty prefactors, and penalty initation distances, respectively. Note that optimal choice of :math:`A_{\mathrm{p}}`, and :math:`d^{m_1(y)}_{\mathrm{p}}` will depend on target conditions and model quality, but respective values of :math:`10^5` kcal/(mol Å :math:`^3`) and 0.01 Å are generally reasonable initial values.

ChIMES smoothing functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _sec-cutoff:

ChIMES models include smoothing functions for each constituent pair interaction, ensuring smoothness as outer cutoffs, :math:`r_{\mathrm{c,out}}^{m_1(y)}` are approached. Currently, two smoothing function forms are supported, "cubic" and "Tersoff". The former has been shown to work reasonably well for models including up to 3-body interactions, but the Tersoff form is best for models including higher-bodied interactions since it can be tuned to minimized the smoothing function cutoff on the overall interactions (see `this paper <https://doi.org/10.1063/5.0021965>`_ for further details).

The cubic cutoff function
""""""""""""""""""""""""""""""

The cubic cutoff function is given by:

.. math::
   :nowrap:

     \begin{equation}
    f_s^{m_1(y)} (r_y) =  \left( 1- \frac{r_y}{r_{\mathrm{c,out}}^{m_1(y)}} \right)^3.
    \end{equation}
    
The Tersoff cutoff function
""""""""""""""""""""""""""""""

The Tersoff cutoff function takes the form:

.. math::
   :nowrap:
   
   \begin{equation}
   f_s^{m_1(y)} (r_y) =  
   \begin{cases}
   0, & \text{if } r_y > r_{\mathrm{c,out}}^{m_1(y)} \\
   1, & \text{if } r_y < d_{\mathrm{t}} \\
   \frac{1}{2}+\frac{1}{2}\sin\left( \pi \left[ \frac{r_y-d_{\mathrm{t}}}{r_{\mathrm{c,out}}^{m_1(y)} - d_{\mathrm{t}} } \right] + \frac{\pi}{2} \right), & \text{otherwise} 
   \end{cases}
   \end{equation} 

where :math:`d_{\mathrm{t}}` is the threshold distance, given by :math:`d_{\mathrm{t}} = r_{\mathrm{c,out}}^{m_1(y)}(1-f_\mathrm{O})`, and :math:`f_\mathrm{O}` is a value in [0,1] generally taken to be 0.5 to 0.75. This form exhibits a smooth step, allowing :math:`^{n}\!E` to remain unmodified by the smoothing function for all :math:`r_{m_1(y)}<d_{\mathrm{t}}`, which is particularly useful for many-body interactions of large :math:`n`, where the product of :math:`\binom{n}{2}f_s^{m_1(y)}(r_y)` factors is used, and can otherwise severely reduce :math:`^{n}\!E` contributions to the total energy.


On permutational invariance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ChIMES explicitly enforces permuational invariance. 

In particular, we require :math:`^n\!E_A = ^n\!E\Pi _A`, where :math:`\Pi` is a permutation operator acting on the :math:`n` atoms in the cluster. This leads to equality conditions among the coefficients, i.e. :math:`c_\mathbf{O}^\mathbf{E} =  c_{\prod \mathbf{O}}^{\prod \mathbf{E}}`. As an example, for a 4-body interaction, we have:

.. math::
   :nowrap:
   
   \begin{equation}
    c_{\alpha \beta \gamma \delta \epsilon \zeta}^{e_ie_j,e_ie_k,e_ie_l,e_je_k,e_je_l,e_ke_l} = c_{\alpha \delta \epsilon \beta \gamma \zeta}^{e_je_i,e_je_k,e_je_l,e_ie_k,e_ie_l,e_ke_l},
    \end{equation}  
    
which is derived by permuting atoms :math:`i` and :math:`j`. We note that in ChIMES, permutational invariance is enforced by treating permutationally related coefficients as the same unique fitting variable.
