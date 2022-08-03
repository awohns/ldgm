.. currentmodule:: ldgm
.. _sec_python_api:

==========
Python API
==========

The ``ldgm`` Python API provides methods for inferring graphical models of SNPs from a given tree sequence. See :ref:`MATLAB API <sec_matlab_api>` for documentation on inferring LDGM precision matrices from LDGMs.


Creating an LDGM
================

.. autofunction:: ldgm.reduce


Intermediate steps to creating an LDGM
======================================

Creating a "bricked tree sequence"
----------------------------------

.. autofunction:: ldgm.brick_ts


Creating a brick graph from a bricked tree sequence
---------------------------------------------------

.. autofunction:: ldgm.brick_graph


Creating a brick graph from a bricked tree sequence
---------------------------------------------------

.. autofunction:: ldgm.reduce_graph


Utility functions
======================================

.. autofunction:: ldgm.prune_sites

   
.. autofunction:: ldgm.return_site_list
