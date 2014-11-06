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

    $wget http://www.clustal.org/download/current/clustalx-2.1-linux-i686-libcppstatic.tar.gz
    $cd clustalw-2.1
    $./configure
    $make

Clone ``alpha-CENTAURI``:

    $cd $CENT_HOME
    $git clone https://github.com/volkansevim/alpha-CENTAURI.git


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

    $ python ../src/chop_to_monomers.py monomers.hmm preads_HuPac.fa > HuPac_monomers.fa

Analyze the higher order structures in the sequence data.

    $ python ../src/monomer_graph_analysis.py preads_HuPac.fa HuPac_monomers.fa 
	
Higher Order Repeat (HOR) Analysis Output 
-------------------	
monomer\_graph\_analysis.py creates four FASTAs and three text files:

* **regularHORs.fa**: Regular HORs extracted from reads. (Leading and trailing partial monomeric sequence in the each read is excluded.) See below for definition of regularity. 
* *irregularHORs.fa**: Irregular HORs.
* **no_HOR_reads.fa**: Reads that contain monomeric sequences but no detectable HOR.
* **too_short_reads.fa**: Reads that are too short to be considered for analysis. Default threshold is 2Kbases.

* **regularHORs_pattern.txt**: Symbolic repeat pattern on each regular HOR, e.g., ABCDABCDABCD. Follows the same order as BASE_regularHORs.fa.
* **regularHORs_pattern.txt**: Symbolic repeat pattern on each irregular HOR. Not meaningful for irregularities caused by non-monomeric insertions.
* **stats.txt**: HOR statistics on reads. See below for detailed description.  

** Read ID format in regularHORs and irregularHORs.fa: **
>__Original\_Read\_ID__\_\_\__length of HOR_\_\_Start base locus in read\_End base locus in read\_\_HOR__n__
__n__ indicates the period of HOR.


