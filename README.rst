======================================================================
Super-simple recipes for getting productive with ScaLAPACK on Windows 
======================================================================
Be sure to check the samples out on `MSDN`_!

About
-----
I wrote these sample programs many years ago to help me benchmark
different implementations of ScaLAPACK on supercomputers running Linux
(vanilla Beowulf clusters as well as SGI Altices running MPT). Over
time, this collection of hacked-together benchmarks evolved into a set
of self-contained teaching samples and validation tests that I used to
help folks in my team come up to speed with ScaLAPACK as well as to
ensure that clusters and interconnects on our customers' site were
property configured.

Last year over the Christmas break, I ported my samples to run on
Windows HPC Server using MSMPI (the Microsoft implementation of the
MPI standard). I also simplified and streamlined a lot of the code
using a few handy features in the C++11 standard. My friends in the
Microsoft HPC team then made the examples available on `MSDN`_ shortly
afterwards.

Compiling the samples
~~~~~~~~~~~~~~~~~~~~~
Although ScaLAPACK is a Fortran 77 library, building and running the
samples does not require you to know Fortran (you don't even need a
Fortran compiler or to even know what Fortran is) - all you need is a
decent C++ compiler, ScaLAPACK and MPI. If you're just starting out, I
would suggest using Microsoft Visual Studio 2010, Intel MKL and
Microsoft MPI respectively since they're the easiest to set up.

For detailed step-by-step instructions on setting up prerequisites,
compiling the samples and running them on a HPC cluster, please refer
to the `document`_ that accompanies the examples on MSDN.

License
-------
The examples are released by Microsoft under the `Apache license`_, version 2.0:

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
   
   http://www.apache.org/licenses/LICENSE-2.0
   
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

.. _`MSDN`: http://code.msdn.microsoft.com/Using-ScaLAPACK-on-Windows-d16a5e76
.. _`document`: http://code.msdn.microsoft.com/Using-ScaLAPACK-on-Windows-d16a5e76/file/55347/1/ScaLAPACKExample.docx
.. _`Apache license`: http://www.apache.org/licenses/LICENSE-2.0.html
