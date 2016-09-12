lwe-frodo parameter scripts
===========================

**lwe-frodo** is a C cryptographic library for post-quantum key exchange based on the learning with errors (LWE) problem.  It is based on the following research paper:

- Joppe Bos, Craig Costello, Léo Ducas, Ilya Mironov, Michael Naehrig, Valeria Nikolaenko, Ananth Raghunathan, Douglas Stebila.  **Frodo: Take off the ring!  Practical, quantum-secure key exchange from LWE**.  In *ACM Conference on Computer and Communications Security (CCS) 2016*, ACM, October, 2016.  DOI:[10.1145/2976749.2978425](http://dx.doi.org/10.1145/2976749.2978425), Eprint [http://eprint.iacr.org/2016/659](http://eprint.iacr.org/2016/659).

The lwe-frodo software is available in the [lwe-frodo/lwe-frodo](https://github.com/lwe-frodo/lwe-frodo) GitHub repository.

This directory contains **Python scripts for selecting parameters**.

You can [download the PDF](https://github.com/lwe-frodo/lwe-frodo/blob/master/LWE-Frodo-full-version.pdf) of the paper from the lwe-frodo GitHub repository.

Files
-----

The main standalone scripts in the directory are:

- **search_params.py** - searches the parameter space (takes about an hour)
- **print_tables.py** - prints out tables for the paper and C code (takes a few minutes)

There are some additional support scripts in the directory:

- pqsec.py - functions for estimating hardness of various LWE assumptions
- renyi.py - functions for computing Renyi divergence
- discrete_distr.py - functions for manipulating discrete distributions
- failure_prob.py - functions for estimating failure probability of key agreement
- approx_distr.py - functions for finding an optimal approximation to rounded Gaussian

License
-------

This software is licensed under the MIT License.  For details, see [LICENSE.txt](https://github.com/lwe-frodo/parameter-selection/blob/master/LICENSE.txt).

Acknowledgements
----------------

JB and LD were supported in part by the Commission of the European Communities through the Horizon 2020 program under project number 645622 (PQCRYPTO).  DS was supported in part by Australian Research Council (ARC) Discovery Project grant DP130104304 and a Natural Sciences and Engineering Research Council of Canada (NSERC) Discovery Grant.  The authors would like to thank Adam Langley, Eric Grosse, and Úlfar Erlingsson for their inputs. A large part of this work was done when VN was an intern with the Google Security and Privacy Research team.
