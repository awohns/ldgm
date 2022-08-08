.. _sec_introduction:

============
Introduction
============

LDGMs are "linkage disequilibrium graphical models": graphical models created from `tree sequences <https://tskit.dev/tutorials/what_is.html>`_ which in turn allow the inference of sparse LDGM precision matrices (LD correlation matrix inverse). The ldgm software package supports (1) :ref:`the creation of LDGMs from tree sequences <sec_python_api>`, (2) :ref:`the inference of LDGM precision matrices <sec_matlab_api>`, and (3) matrix algebra operations that use LDGM precision matrices to massively speed up statistical genetic computations. Please see our forthcoming preprint for more details on LDGMs.

We provide LDGMs for the five continential ancestries represented in the `1000 Genomes Project <http://www.internationalgenome.org>`_. The LDGMs can be downloaded from (link) and the LDGM precision matrices from (link). This documentation contains tutorials and guidelines for using these files in your applications. Most users will be using these precomputed LDGMs in their statistical genetics methods, rather than creating new LDGMs.
