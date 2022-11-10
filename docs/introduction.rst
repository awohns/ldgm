.. _sec_introduction:

============
Introduction
============

LDGMs are "linkage disequilibrium graphical models": graphical models created from `tree sequences <https://tskit.dev/tutorials/what_is.html>`_ which in turn allow the inference of sparse LDGM precision matrices (LD correlation matrix inverse). The ``ldgm`` software package supports matrix algebra operations that use LDGM precision matrices to massively speed up statistical genetic computations. This code is currently implemented in MATLAB (a Python implementation is forthcoming).

We provide LDGMs for the five continential ancestries represented in the `1000 Genomes Project <http://www.internationalgenome.org>`_. The LDGMs and the LDGM precision matrices can be `downloaded from here <https://www.dropbox.com/sh/1huaxgad2bjjv9a/AACI6pRZpr6gmK5E1mo4fck6a?dl=0>`_. This documentation contains tutorials and guidelines for using these files in your applications. Most users will be using these precomputed LDGMs in their statistical genetics methods, rather than creating new LDGMs.

This package also supports :ref:`the creation of LDGMs from tree sequences <sec_python_api>` (written in Python) and :ref:`the inference of LDGM precision matrices <sec_matlab_api>` (written in MATLAB).

Please see our forthcoming preprint for more details on LDGMs.


`Link text <link URL>`_
