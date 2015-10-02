Installing and running alpha-CENTAURI
===========================================

:Contributors: 

Volkan Sevim, Jason Chin, Ali Bashir, Karen Miga

:Version: 

0.2 of 2014/12/10

alpha-CENTAURI is a Pyhton package for mining alpha satellites and their higher-order structures in sequence data. It requires an initial consensus sequence (a sample is provided in the package). This consensus sequence is used to build an HMM model, which is employed to detect alpha-satellite monomers in the sequence data.

alpha-CENTAURI can in principle run on any sequence data, however, its performance will increase greatly if the reads are filtered based on similarity to the initial monomer set. 

Installation
--------------

Prerequisites (you might need superuser privilages to install these packages):

* g++ compiler
* python2.7
* python-dev
* git ( http://git-scm.com/ for installing the dependencies from github )
* python pip 
* virtualenv

Make sure you are using python2.7. First create a clean virtualenv and activate it:

    $ export CENT_HOME=/some/path/to/your/CENT_ENV
    $ virtualenv -p /usr/bin/python2.7 $CENT_HOME
    $ cd $CENT_HOME
    $ . bin/activate

Install ``numpy``:

    $pip install numpy

Install ``networkx``

    $ pip install networkx

Install PacBio package ``Falcon`` (Install the commit specified below to avoid a future incompatibility.)

    $ pip install git+https://github.com/PacificBiosciences/FALCON.git@96230ec9d6027e465deaccdb6fe3c045e5b820a3#falcon
 
Install ``HMMer`` (Instructions are intended for the 32-bit version. If you have a 64-bit system, locate the corresponding file on the HMMER website, and install as explained below.)

    $ wget http://selab.janelia.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-ia32.tar.gz
    $ tar -xvf hmmer-3.1b1-linux-intel-ia32.gz.tar
    $ cd hmmer-3.1b1-linux-intel-ia32
    $ ./configure
    $ make
    $ make check
    $ make install

Clone ``alpha-CENTAURI``:

    $ cd $CENT_HOME
    $ git clone https://github.com/volkansevim/alpha-CENTAURI.git


Workflow Example
--------------
For this example, we will use the dataset under ``example`` folder provided in the package. There are 5 files in this folder:

* pread_HuPac_example.fa: A filtered set of reads that contain sequences similar to the provided monomers.
* alpha.sto, alpha.rc.sto : Human alpha-satellite repeat consensus sequence aligned to itself, and its reverse complement. 
The files below are the outputs of the step 1. of the workflow. They are provided for convenience.

* alpha.hmm, alpha.rc.hmm: HMM built using the multiple sequence alignment of the consensus and its reverse complement.

### Workflow Steps

For this workflow, you will need a consensus sequence for the monomers in your repeats. The HMM needs two files: consensus sequence and its reverse complement aligned to themselves. If you have your own consensus sequence, you can just modify the .sto files provided in the package using that sequence.

First, build an HMM based on the alignment.

    $ hmmbuild alpha.hmm alpha.sto
    $ hmmbuild alpha.rc.hmm alpha.rc.sto
 
Infer monomers from sequence data using the HMM, write them into inferred_monomers.fa. 

    $ python ../src/chop_to_monomers.py pread_HuPac_example.fa alpha.hmm alpha.rc.hmm 

(Here minimum monomer length is assumed 150bp, and shorter inferred monomers are discarded. Use the -l flag to modify this number in order to analyze repeats other than alpha satellites. Use -h flag for help.)

Analyze the higher order structures in the sequence data. 

    $ python ../src/monomer_graph_analysis.py pread_HuPac_example.fa inferred_monomers.fa
    
This script is pre-tuned for analyzing alpha-satellite repeats. Use the command-line arguments below to modify the analysis parameters. (Use -h flag for help.)

  -l: Average length of a monomer.
  -d: Maximum allowed head-to-tail distance between two adjacent monomers.
  -s: Minimum allowed read length.
  -t: Specifies a clustering threshold. Multiple -t allowed. Values sorted and tested in descending order.

Default clustering threshold list is 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.89, 0.88. Values are tested in descending order, until an HOR is detected. 
In order to specify a different (set of) threshold(s) use the -t flag. For example,  

    $ python ../src/monomer_graph_analysis.py pread_HuPac_example.fa inferred_monomers.fa -t 0.95 -t 0.93 -t 0.90

would test threshold vales 0.95, 0.93, and 0.90 in that order (i.e., specified list is tested in descending order).
	
Higher Order Repeat (HOR) Analysis Output 
-------------------	
monomer\_graph\_analysis.py creates four FASTAs and three text files:

* **regularHORs.fa**: Regular HORs extracted from reads. (Leading and trailing partial monomeric sequence in the each read is excluded.) See below for definition of regularity. 
* **irregularHORs.fa**: Irregular HORs.
* **no_HOR_reads.fa**: Reads that contain monomeric sequences but no detectable HOR.
* **too_short_reads.fa**: Reads that are too short to be considered for analysis. Default threshold is 2Kbases.
* **inversions.fa**: Reads that contain an inversion. 
* **missing_monomer.fa**: Reads with a missing monomer. See KNOWN ISSUES for details. 
* **regularHORs_pattern.txt**: Symbolic repeat pattern on each regular HOR, e.g., ABCDABCDABCD. Follows the same order as BASE_regularHORs.fa.
* **irregularHORs_pattern.txt**: Symbolic repeat pattern on each irregular HOR. Not meaningful for irregularities caused by non-monomeric insertions.
* **inversions_pattern.txt**: Symbolic repeat pattern on each regular HOR, e.g., ABCDABCDABCD. Follows the same order as BASE_regularHORs.fa.
* **stats.txt**: HOR statistics on reads. See below for detailed description.  

**Read ID format in regularHORs.fa and irregularHORs.fa:** 

_OriginalID_ \_\_\_ _length_ \_\_ _start_ \_ _end_ \_\_HOR _n_

Here, _OriginalID_ indicates the ID of the original read, _n_ indicates the period, and _length_ indicates the length of complete HOR structure. _start_ and _end_ indicate the first and last base positions of the HOR structure in the original read. 

Example: `6ed935a_20072_0___8646__102_8369__HOR8`

_OriginalID_=6ed935a\_20072\_0, _length_=8646, _start_=102, _end_=8369, _HOR period_=8

### KNOWN ISSUES 
As one lowers the clustering threshold, repeat structure tends to converge onto 
a dimer, i.e., HOR2. Thus, some HOR2s classified as regular can be irregular 
repeats. We recommend a visual inspection of all regular HOR2 predictions. An 
alternative is to raise the lowest clustering threshold, however, that could 
result in missing some regular repeats with periods longer that 2.

Some reads that are classified as irregular are in fact regular. Ther are two 
reason for this misclassification:

(a) One or more monomers in the read are not recognized by the HMM.

(b) HOR unit contains more multiple instances of a certain monomer.

Currently, cases in (a) are reported seperately in missing_monomer.fa. The 
file contains both regular and irregular reads with an unrecognized monomer. 
Also, some reads with TE insertions could be potentially found in this file 
as the current algorithm cannot distinguish a TE from an unrecognized monomer.

We will improve the workflow in the next version of the software to correctly 
classify such reads.
### Stats.txt Content 
Please see the publication for details about the algorithm.

**RID**: Read ID

**Regularity**: R=regular, I=irregular, N=No HOR detected, V=Inversion

**Read Len**: Read length

**Thresh**: Threshold used for detecting monomer identity. Algorithm starts from a high threshold and gradually reduces it.

**#all monomers**: Number of all detected monomers, clustered or not clustered. 

**#mono in a cluster**: Total number of clustered monomers. 

**Isolates**: Number of monomers that do not belong to a cluster.  

**Clustered monomer fraction in read**: Fraction of monomers that belong to a cluster.

**#total clusters**: Total number of clusters, i.e., number of distinct monomers in HOR.

**mean identity within clusters**: Mean identity for monomers of the same kind, averaged over all clusters.

**mean identity between clusters**: Mean identity for monomers of the different kind.

**min, max, median monomeric period**: Median monomeric period is the median distance between monomers of the same kind.

**Normalized min, max monomeric period**: Minimum (maximum) distance between monomers of the same kind normalized by the median distance. 

**min, max, median head to tail interval**: Number of bases between two adjacent monomers. This is generally 0. A large number indicates an insertion.  

