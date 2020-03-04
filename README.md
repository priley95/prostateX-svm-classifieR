# Multi-parametric classification of prostate cancer through mpMRI

This repository realtes to the MSc project. This is an implementation in R. There is another which is in Python but utilises shared code, and forms the output of a conference paper.

**Note** this was on image patch data, extracted from the mpMRI. Some associated patient-level metadata was also analysed - the age and weight of the patient & the prostate zone in which the lesion is located.

The .pdf of the MSc project is also provided.

_Pre-process_

A snippet of the code is provided; as various MRI measurements were taken, it's quite a lot! Much work was required to show how the various scan parameters needed to be combined as appropriate. 

_Methods_

This script shows an SVM implementation using the **mlr** package in R.
