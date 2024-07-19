Introduction
============

Overview
--------

`OP-PIC <https://github.com/OP-DSL/OP-PIC>`_ is a high-level embedded domain specific language (eDSL) for writing **unstructured mesh** Particle in Cell algorithms with automatic parellelisation on multi-core and many-core architectures. The API is embedded in C/C++.

The current OP-PIC eDSL supports generating code targeting multi-core CPUs with OpenMP threading, many-core GPUs with CUDA offloading, and distributed memory cluster variants of these using MPI. 
Particularly parallelizations are, sequential, OpenMP, MPI, CUDA, HIP and distributed memory CUDA+MPI, HIP+MPI.

These pages provide detailed documentation on using OP-PIC, including an installation guide, an overview of the C++ API, a walkthrough of the development of an example application, and developer documentation.

Licencing
---------

OP-PIC is released as an open-source project under the BSD 3-Clause License. See the `LICENSE <https://github.com/OP-DSL/OP-PIC/blob/main/LICENSE>`_ file for more information.
