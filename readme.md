README - Phase Retrieval code v20150728

https://github.com/JeffFessler/PhaseRetrievalWithOutliers

Created by Daniel Weller, University of Virginia, and
           Jeffrey Fessler, University of Michigan

Please direct questions and comments to Daniel Weller 
at d.s.weller@ieee.org.

The contents of this program are subject to copyright 
as set forth in the file LICENSE.TXT accompanying this
software. If you did not receive a copy of LICENSE.TXT, 
please email d.s.weller@ieee.org.

## I. Introduction

The enclosed program includes implementations of the 
phase retrieval algorithm proposed in the paper 

D. S. Weller, A. Pnueli, G. Divon, O. Radzyner, 
Y. C. Eldar, and J. A. Fessler,
["Undersampled Phase Retrieval with Outliers,"
IEEE Transactions on Computational Imaging, 1(4):247-58, Dec. 2015.](http://doi.org/10.1109/TCI.2015.2498402)

A draft version of the manuscript appeared on
[arXiv: 1411.2183 [stat.AP]](http://arxiv.org/abs/1411.2183).

We ask that any user please cite this work in any 
publications that present results based on this software.

## II. Installation Instructions

This code has been tested with MATLAB (The Mathworks, 
Inc.) versions R2014a and R2015a. Please note that 
some major changes in the graphics subsystem between 
R2014a and R2015a means that some of the plotting code 
no longer functions equivalently across platforms.
In addition, other changes in more recent releases 
may no longer be compatible with this code.

To use this code, you need to install the
Michigan Image Reconstruction Toolbox (MIRT) provided from 
http://web.eecs.umich.edu/~fessler/code/index.html.

Installation instructions for that toolbox can be
found on the website or in the `readme` file included 
with the toolbox. Once the toolbox is installed, 
you should either include the toolbox directories 
on your MATLAB path at startup, or be sure to call 
the `setup_IRT.m` function in this software before 
running the included scripts.

Once MIRT is installed, you should be almost ready to go.
Clone this repo
and redirect the MATLAB working directory to the `src/` folder.
You should now be able to run any of the
`run_*.m` or `plot_*.m` scripts therein
after following the directions in that script file. 
Note that some scripts require that you have run others first.

A quick note about the "results" directory: the 
scripts will not attempt to create a results 
directory to store simulation data. You must 
create a readable/writable directory and update 
the path "outfolder" in your scripts to point 
to this folder.


## III. Running the Code

A number of scripts are included with this software,
mainly with the goal of recreating figures from the paper above.
To generate figures, run in order:

* `surrogate_examples.m`
* `run_phase_retrieval_padmm_mus.m`
* `run_phase_retrieval_betas.m`
* `plot_phase_retrieval_betas.m`
* `run_phase_retrieval_MC.m`
* `run_sparse_fienup_MC.m`
* `plot_phase_retrieval_MC.m`
* `run_phase_retrieval_MC_noise.m` (with settings varying outliers)
* `run_sparse_fienup_MC_noise.m` (with settings varying outliers)
* `plot_phase_retrieval_MC_outliers.m`
* `run_phase_retrieval_MC_noise.m` (with settings varying noise levels)
* `run_sparse_fienup_MC_noise.m` (with settings varying noise levels)
* `plot_phase_retrieval_MC_noise.m`
* `run_phase_retrieval_SOD.m`
* `run_sparse_fienup_SOD.m`
* `plot_phase_retrieval_SOD.m`

Each script contains a description of its purpose
and which figures it can help recreate.
Note that many of the figures in the paper 
contain results from competing methods such 
as GESPAR and PR-GAMP that are not included with this code.
Please contact the authors of those methods
to obtain software to test them directly.


## IV. Adapting this Code

The main functionality for the proposed method 
requires these files: 

* `phase_retrieval_p1_q2.m`
* `ox_fessler_p1_q2.m`
* `PADMM_1split_p1_q2.m`
* `surrogate_p1_q2.m`
* `shrink.m`
* `compute_beta_normalization.m`

The first three files contain the three levels 
of the algorithm (initializations, MM, and ADMM). 
The last three are utility functions
used to evaluate the majorizer, perform soft-thresholding, 
and compute the normalization factor for beta, 
respectively.
We hope the scripts described in Part III above
provide adequate examples
of how to use these functions for your own code. 

Experimental code for the quadratic data fit case is also included.
Those files have `p2_q2` in their names.

Enjoy!

Daniel Weller
Assistant Professor
University of Virginia
d.s.weller@ieee.org

July 30, 2015

This github version
created by Jeff Fessler on 2022-08-16,
using code previously downloaded on 2016-05-03
from http://people.virginia.edu/~dsw8c/sw.html.
