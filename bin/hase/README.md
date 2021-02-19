# HASE
Framework for efficient high-dimensional association analyses. This is a fork of the [original repository](https://github.com/roshchupkin/hase).

## Installation HASE

Navigate to directory where you want to install HASE and clone this repository:
     ```
     git clone https://github.com/roshchupkin/hase.git
     ```

## requirements

1. Install the [HDF5](https://www.hdfgroup.org/downloads/hdf5/) software.
   The HDF5_DIR has to be added to the environment variables.
    ```
    export HDF5_DIR=~/hdf5<version-number>/hdf5/
    ```
2. HASE uses [python](https://www.python.org/downloads/) version 2.7
3. To install the required packages, use `pip install -r requirements.txt`. Where `requirements.txt` is the file located within the root folder of this repository.
 
## User Guide
Use [this wiki](https://github.com/roshchupkin/hase/wiki) for the, original, upstream wiki. Or [this wiki](https://github.com/CAWarmerdam/hase/wiki/Running-HASE-meta-analysis-in-the-example-study) for a guide on using hase with the example data.

### Usage notes
- Variants that do not have any variance across the samples are to be removed as this will cause an exception.
- Tests have made apparent that converting VCF and PLINK files not always result in correct HDF5 files. It is recommended to check this prior to proceeding analysis.

## Changes from the upstream repository
- Fixed bug causing an exception when more than 1000 individuals were used.
- Resolved bug causing the `--intercept` option having no effect.
- Made version numbers of pip packages explicit.
- Added commentary to code in places.
- Lines 355-357 of hase.py were commented out because this caused pipeline to crash when >1 datasets were added.
- Line 355 of /hdgwas/data.py were changed `self.chunk_size=10000 --> self.chunk_size=20000`.
- For eQTLGen pipelines: removed folders with unit tests and test data, in order to keep it lightweight.

### Interactions development branch
- Implemented the possibility for using interaction terms.
- Started implementing tests for both using, and not using interaction terms.

## Citation 
If you use HASE framework, please cite:

[Roshchupkin, G. V. et al. HASE: Framework for efficient high-dimensional association analyses. Sci. Rep. 6, 36076; doi: 10.1038/srep36076 (2016)](http://www.nature.com/articles/srep36076) 

## Licence
This project is licensed under GNU GPL v3.

## Authors
Gennady V. Roshchupkin (Department of Epidemiology, Radiology and Medical Informatics, Erasmus MC, Rotterdam, Netherlands)

Hieab H. Adams (Department of Epidemiology, Erasmus MC, Rotterdam, Netherlands) 

## Contacts

If you have any questions/suggestions/comments or problems do not hesitate to contact us!

* g.roshchupkin@erasmusmc.nl
* h.adams@erasmusmc.nl
 
