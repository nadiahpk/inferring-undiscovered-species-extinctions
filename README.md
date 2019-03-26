# Inferring undiscovered species' extinctions

## Overview

This repository stores data and code (mostly Python3, with some interface to R) that was used to infer the number of vascular plant species that went extinct in Singapore since colonialisation, including those that went extinct before they could be discovered, using the SEUX model (Chisholm et al. 2016). 

**See also**: ['seux' R package](https://github.com/nadiahpk/seux "seux")

## About the manuscript

Kristensen, N.P., ... Chishom, R.A. (in prep) Estimating extinctions of undiscovered plant species in Singapore

**Abstract**

The highly developed island city-state of Singapore has one of the best-documented tropical floras in the world.  We collated museum and herbarium records of plant collections in Singapore since 1796, with the goal of estimating the total extinction rate since 1822.  Our database includes 2090 native species, of which 501 (24%) are now apparently locally extinct.  We used a statistical model and two inference techniques -- classical and Bayesian -- to estimate the number of undetected extinctions.  We estimate that between 319 and 544 additional species went extinct before they were ever discovered, suggesting a total extinction rate of 34--37%.  These numbers serve to illustrate the magnitudes of local species loss that can be expected in the tropics as extensive development continues.

## About the repository

The diagram below gives an overview of the directory structure, and highlights important data and results files. The Singapore plants database can be found in `merged.csv`, and the input to the model in `first last detns final.csv`. The method and results for classical and Bayesian inference can be found in the `classical` and `mcmc` directories, respectively.

![repository structure](https://raw.githubusercontent.com/nadiahpk/inferring-undiscovered-species-extinctions/master/repo_structure.png)

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or distribute this software, either in source code form or as a compiled binary, for any purpose, commercial or non-commercial, and by any means.

In jurisdictions that recognize copyright laws, the author or authors of this software dedicate any and all copyright interest in the software to the public domain. We make this dedication for the benefit of the public at large and to the detriment of our heirs and successors. We intend this dedication to be an overt act of relinquishment in perpetuity of all present and future rights to this software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
