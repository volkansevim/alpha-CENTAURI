Installing and running alpha-CENTAURI
===========================================

:Authors: 
    Volkan Sevim & Jason Chin

:Version: 0.1.0 of 2014/08/05

alpha-CENTAURI is a Pyhton package for mining alpha satellites and their higher-order structures in sequence data. It requires an initial set of "monomer" sequences (a set for human genome is provided in the package). This initial set is used to build an HMM model, which is employed to detect alpha-satellite monomers in the sequence data.

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

Compile ``libhdf5`` (
http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.9/src/ ) and install it
in under virtualenv:

    $ wget http://www.hdfgroup.org/ftp/HDF5/prev-releases/hdf5-1.8.9/src/hdf5-1.8.9.tar.gz
    $ tar zxvf hdf5-1.8.9.tar.gz
    $ cd hdf5-1.8.9
    $ ./configure --prefix=$HBAR_HOME --enable-cxx
    $ make install
    $ cd ..
    $export HDF5INCLUDEDIR=/home/CENT_ENV/include/
    $export HDF5LIBDIR=/home/CENT_ENV/lib/

Install ``numpy``:

    $pip install numpy

Install ``h5py`` ( http://h5py.googlecode.com/files/h5py-2.0.1.tar.gz )::

    $ wget http://h5py.googlecode.com/files/h5py-2.0.1.tar.gz
    $ tar zxvf h5py-2.0.1.tar.gz
    $ cd h5py-2.0.1
    $ python setup.py build --hdf5=$HBAR_HOME
    $ python setup.py install

Install ``networkx``

    $ pip install networkx

Install PAcBio packages ``PBcore`` and ``Falcon``

    $ pip install git+https://github.com/PacificBiosciences/pbcore.git#pbcore
    $ pip install git+https://github.com/PacificBiosciences/FALCON.git#falcon
 
Install ``HMMer``

    $ wget http://selab.janelia.org/software/hmmer3/3.1b1/hmmer-3.1b1-linux-intel-ia32.gz.tar
    $ tar -xvf hmmer-3.1b1-linux-intel-ia32.gz.tar
    $ ./configure
    $ make
    $ make check
    $ make install

Install ``ClustalW``

    $ wget http://www.clustal.org/download/current/clustalw-2.1.tar.gz
    $ tar -xvzf clustalw-2.1.tar.gz
    $ cd clustalw-2.1
    $ ./configure
    $ make

Clone ``alpha-CENTAURI``:

    $ cd $CENT_HOME
    $ git clone https://github.com/volkansevim/alpha-CENTAURI.git


Workflow Example
--------------
For this example, we will use the dataset under ``example`` folder provided in the package. There are 9 files in this folder:

* MigaKH.HigherOrderRptMon.fa: Initial set of monomers from the human genome.
* pread_HuPac_example.fa: A filtered set of reads that contain sequences similar to the provided monomers.

The files below are the outputs of the steps 1. and 2. of the workflow. They are provided for convenience.

* MigaKH.HigherOrderRptMon.aln and .dnd: Multiple sequence alignments for the initial set of monomers.
* monomers.hmm, .h3m, .h3i, h3f, h3p: HMM built using the multiple sequence alignment of the initial set of monomers.

### Workflow Steps

Generate multiple sequence alignments on the initial set of monomers 

    $ cd alpha-CENTAURI/example/
    $ $CENT_HOME/clustal-2.1/src/clustalW MigaKH.HigherOrderRptMon.fa

Build an HMM based on the alignment.

    $ hmmbuild monomers.hmm MigaKH.HigherOrderRptMon.aln
    $ hmmpress monomers.hmm

Infer monomers from sequence data using the HMM, write them into HuPac_monomers.fa. 

    $  python ../src/chop_to_monomers.py -f pread_HuPac_example.fa -m monomers.hmm > HuPac_monomers.fa

(Here minimum monomer length is assumed 160bp. Use the -l flag to modify it to analyze repeats other than alpha satellites. Use -h flag for help.)

Analyze the higher order structures in the sequence data.

    $ python ../src/monomer_graph_analysis.py -f pread_HuPac_example.fa -m HuPac_monomers.fa

This script is pre-tuned for analyzing alpha-satellite repeats. Use the command-line arguments below to modify the analysis parameters. (Use -h flag for help.)

  -l: Average length of a monomer.
  -d: Maximum allowed head-to-tail distance between two adjacent monomers.
  -s: Minimum allowed read length.
  -t: Specifies a clustering threshold. Multiple -t allowed. Values sorted and tested in descending order.

Default clustering threshold list is 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.89, 0.88. Values are tested in descending order, until an HOR is detected. 
In order to specify a different (set of) threshold(s) use the -t flag. For example,  

    $ python ../src/monomer_graph_analysis.py -f pread_HuPac_example.fa -m HuPac_monomers.fa -t 0.95 -t 0.93 -t 0.90

would test threshold vales 0.95, 0.93, and 0.90 in that order.
	
Higher Order Repeat (HOR) Analysis Output 
-------------------	
monomer\_graph\_analysis.py creates four FASTAs and three text files:

* **regularHORs.fa**: Regular HORs extracted from reads. (Leading and trailing partial monomeric sequence in the each read is excluded.) See below for definition of regularity. 
* **irregularHORs.fa**: Irregular HORs.
* **no_HOR_reads.fa**: Reads that contain monomeric sequences but no detectable HOR.
* **too_short_reads.fa**: Reads that are too short to be considered for analysis. Default threshold is 2Kbases.
* **inversions.fa**: Reads that contain an inversion. 
* **regularHORs_pattern.txt**: Symbolic repeat pattern on each regular HOR, e.g., ABCDABCDABCD. Follows the same order as BASE_regularHORs.fa.
* **irregularHORs_pattern.txt**: Symbolic repeat pattern on each irregular HOR. Not meaningful for irregularities caused by non-monomeric insertions.
* **inversions_pattern.txt**: Symbolic repeat pattern on each regular HOR, e.g., ABCDABCDABCD. Follows the same order as BASE_regularHORs.fa.
* **stats.txt**: HOR statistics on reads. See below for detailed description.  

**Read ID format in regularHORs.fa and irregularHORs.fa:** 

_OriginalID_ \_\_\_ _length_ \_\_ _start_ \_ _end_ \_\_HOR _n_

Here, _OriginalID_ indicates the ID of the original read, _n_ indicates the period, and _length_ indicates the length of complete HOR structure. _start_ and _end_ indicate the first and last base positions of the HOR structure in the original read. 

Example: `6ed935a_20072_0___8646__102_8369__HOR8`

_OriginalID_=6ed935a\_20072\_0, _length_=8646, _start_=102, _end_=8369, _HOR period_=8

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

**median, min, max #monomers in HOR**: Median, min, and max of interval/average_monomer_length. 
Here 'interval' is an array containing distances between the monomers within the same cluster. 
Min, max, and median are identical in a regular HOR. (Also see 'min, max, median monomeric period' below.)

**median, min, max monomer len**:  Lengths of detected monomers can be different. Large deviations from 171bp points to an irregularity.

**mean identity within clusters**: Mean identity for monomers of the same kind, averaged over all clusters.

**mean identity between clusters**: Mean identity for monomers of the different kind.

**min, max, median monomeric period**: Median monomeric period is the median distance between monomers of the same kind.

**Normalized min, max monomeric period**: Minimum (maximum) distance between monomers of the same kind normalized by the median distance. 

**min, max, median head to tail interval**: Head to tail distance is the number of bases between two adjacent monomers. This is generally 0. A large number indicates an insertion.  

