# Inferring undiscovered species' extinctions

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3733920.svg)](https://doi.org/10.5281/zenodo.3733920)

## Overview

This repository stores data and code (mostly Python3, with some interface to R) that was used to infer the number of vascular plant species that went extinct in Singapore since colonialisation, including those that went extinct before they could be discovered, using the SEUX model (Chisholm et al. 2016). 

**See also**: ['seux' R package](https://github.com/nadiahpk/seux "seux")

## About the manuscript

Kristensen, N.P., Seah, W.W., Chong, K.Y., Yeoh, Y.S., Fung, T., Berman, L.M., Tan, H.Z., Chishom, R.A. (2020) Extinction rate of discovered and undiscovered plants in Singapore, Conservation Biology

**Abstract**

Extinction is a key issue in the assessment of global biodiversity.  However, many extinction-rate measures do not account for species that went extinct before they could be discovered.  The highly developed island city-state of Singapore has one of the best-documented tropical floras in the world. This allowed us to estimate the total rate of local floristic extinctions in Singapore since 1822 after accounting for sampling effort and cryptoextinctions by collating herbaria records.  Our database comprised 34224 specimens from 2076 native species, of which 464 species (22%) were considered nationally extinct.  We assumed that undiscovered species had the same annual per-species extinction rates as discovered species and that no undiscovered species remained extant.  With classical and Bayesian algorithms, we estimated that, respectively, 304 (95% CI 213-414) and 412 (95% credible interval 313-534) additional species went extinct before they could be discovered; thus, corresponding total extinction rate estimates were 32% and 35% (range 30–38%).  We detected violations of our two assumptions that could cause our extinction estimates, particularly the absolute numbers, to be biased downward.  Thus, our estimates should be treated as lower bounds.  Our results illustrate the possible magnitudes of plant extirpations that can be expected in the tropics as development continues.

## About the repository

The diagram below gives an overview of the directory structure, and highlights important data and results files. The Singapore plants database can be found in `merged.csv`, and the input to the model in `first last detns final.csv`. The method and results for classical and Bayesian inference can be found in the `classical` and `mcmc` directories, respectively.

![repository structure](https://raw.githubusercontent.com/nadiahpk/inferring-undiscovered-species-extinctions/master/repo_structure.png)

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or distribute this software, either in source code form or as a compiled binary, for any purpose, commercial or non-commercial, and by any means.

In jurisdictions that recognize copyright laws, the author or authors of this software dedicate any and all copyright interest in the software to the public domain. We make this dedication for the benefit of the public at large and to the detriment of our heirs and successors. We intend this dedication to be an overt act of relinquishment in perpetuity of all present and future rights to this software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
