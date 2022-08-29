.. currentmodule:: ldgm
.. _sec_python_api:

==========
Python API
==========

The ``ldgm`` Python API provides methods for inferring an LDGM from a tree sequence. See the :ref:`MATLAB API page<sec_matlab_api>` for documentation on inferring LDGM precision matrices from LDGMs.


Creating an LDGM
================

.. autofunction:: ldgm.make_ldgm


Intermediate steps to creating an LDGM
======================================

Creating a "bricked tree sequence"
----------------------------------

.. autofunction:: ldgm.brick_ts


Creating a brick graph from a bricked tree sequence
---------------------------------------------------

.. autofunction:: ldgm.brick_haplo_graph


Creating a brick graph from a bricked tree sequence
---------------------------------------------------

.. autofunction:: ldgm.reduce_graph


Utility functions
======================================

.. autofunction:: ldgm.prune_sites

   
.. autofunction:: ldgm.make_snplist
