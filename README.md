Installing and running alpha-CENTAURI
===========================================

:Authors: 
    Volkan Sevim & Jason Chin

:Version: 0.1.0 of 2014/08/05


Prerequisites:

* python2.7
* git ( http://git-scm.com/ for installing the dependencies from github )

Install Python
--------------

Make sure you are using python2.7. We can create a clean virtualenv and
activate it::

    $ export CENT_HOME=/some/path/to/your/CENT_ENV
    $ virtualenv -p /usr/bin/python2.7 $CENT_HOME
    $ cd $CENT_HOME
    $ . bin/activate

Next you want to install ``pbcore`` ( http://www.numpy.org ) library and its
dependencies. First install ``numpy``::

    $ pip install numpy

Install Python Libraries
-----------------------------

Next, install the PacBio python libraries::
    

    $ pip install git+https://github.com/PacificBiosciences/pbcore.git#pbcore
  
If you do a ``pip freeze``, this is what you will see::

    $ pip freeze
    h5py==2.0.1
    html5lib==0.95
    isodate==0.4.9
    matplotlib==1.2.0
    numpy==1.6.2
    pbcore==0.6.0
    pbtools.hbar-dtk==0.1.0
    pbtools.pbdagcon==0.2.0
    pbtools.pbh5tools==0.75.0
    pyparsing==1.5.7
    pypeflow==0.1.0
    rdfextras==0.4
    rdflib==3.4.0
    wsgiref==0.1.2

Install Other Dependencies
--------------------------------

We need BLASR for the pre-assembly mapping. BLASR is included in the SMRT(R)
Analysis installation and is also available on github. You need to copy a blasr
binary into your ``$CENT_HOME/bin``::

    $ cp blasr $CENT_HOME/bin

Last, we need a copy of Celera Assembler for the assembly itself::

    $ wget http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-7.0/wgs-7.0-PacBio-Linux-amd64.tar.bz2
    $ tar jxvf wgs-7.0-PacBio-Linux-amd64.tar.bz2 -C $CENT_HOME/bin/
    $ ln -sf $CENT_HOME/bin/wgs-7.0/Linux-amd64/bin/* $CENT_HOME/bin/
 
For the final step of polishing the assembly using Quiver, we need a SMRT
Analysis installation, plus Quiver from github. Quiver is available via a link
from PacBio DevNet or directly on github. Please follow the installation
instructions there.

Running CENT_WF.py
=================

Warning
---------

- While the general strategy of HBAR will work for larger genome in principle.
  Special consideration should be taken to do the distributed computing
  efficiently.

Set up the environment
-----------------------

Make sure you have clean UNIX shell environment. (Please be sure you do not
have ``PYTHON_PATH`` environment variable and other random non-standard paths
in your ``PATH`` environment variable.) If your shell environment is clean, do::

    $ export PATH_TO_CENT_ENV=/the_full_path_to_your_installation
    $ source $PATH_TO_CENT_ENV/bin/activate

You can "deactivate" the ``CENT_ENV`` by::
 
    $ deactivate

Prepare data, set up the configuration and run
----------------------------------------------

Prepare a working directory and create a file ``input.fofn`` that points to the
base files (``bas.h5`` files) for assembly. Let call this directory
``my_assembly``.  You also need to make sure the paths in the ``input.fofn``
file are absolute and not relative paths.

Here is an example of the ``input.fofn`` files::

    /mnt/data/m120803_022519_42141_c100388772550000001523034210251234_s1_p0.bas.h5
    /mnt/data/m120803_041200_42141_c100388772550000001523034210251235_s1_p0.bas.h5
    /mnt/data/m120803_055858_42141_c100388772550000001523034210251236_s1_p0.bas.h5
    /mnt/data/m120803_074648_42141_c100388772550000001523034210251237_s1_p0.bas.h5

Copy the example configuration to the working directory::

    $ cd my_assembly
    $ cp $PATH_TO_CENT_ENV/etc/HBAR.cfg .

Here is the content of ``HBAR.cfg``::

    [General]
    # list of files of the initial bas.h5 files
    input_fofn = input.fofn

    # The length cutoff used for seed reads used for initial mapping
    length_cutoff = 4500

    # The length cutoff used for seed reads usef for pre-assembly
    length_cutoff_pr = 4500

    # The read quality cutoff used for seed reads
    RQ_threshold = 0.75

    # SGE job option for distributed mapping 
    sge_option_dm = -pe smp 8 -q fas

    # SGE job option for m4 filtering
    sge_option_mf = -pe smp 4 -q fas

    # SGE job option for pre-assembly
    sge_option_pa = -pe smp 16 -q fas

    # SGE job option for CA 
    sge_option_ca = -pe smp 4 -q fas

    # SGE job option for Quiver
    sge_option_qv = -pe smp 16 -q fas

    # SGE job option for "qsub -sync y" to sync jobs in the different stages
    sge_option_ck = -pe smp 1 -q fas 

    # blasr for initial read-read mapping for each chunck (do not specific the "-out" option). 
    # One might need to tune the bestn parameter to match the number of distributed chunks to get more optimized results 
    blasr_opt = -nCandidates 50 -minMatch 12 -maxLCPLength 15 -bestn 4 -minPctIdentity 70.0 -maxScore -1000 -nproc 4 -noSplitSubreads

    #This is used for running quiver
    SEYMOUR_HOME = /mnt/secondary/Smrtpipe/builds/Assembly_Mainline_Nightly_Archive/build470-116466/

    #The number of best alignment hits used for pre-assembly
    #It should be about the same as the final PLR coverage, slight higher might be OK.
    bestn = 36

    # target choices are "pre_assembly", "draft_assembly", "all"
    # "pre_assembly" : generate pre_assembly for any long read assembler to use
    # "draft_assembly": automatic submit CA assembly job when pre-assembly is done
    # "all" : submit job for using Quiver to do final polish
    target = draft_assembly

    # number of chunks for distributed mapping
    preassembly_num_chunk = 8 

    # number of chunks for pre-assembly. 
    # One might want to use bigger chunk data sizes (smaller dist_map_num_chunk) to 
    # take the advantage of the suffix array index used by blasr
    # It would be great to use ramdisk for this. Set tmpdir to a NFS mount will probably have very bad performance.
    dist_map_num_chunk = 4

    # "tmpdir" is for preassembly. A lot of small files are created and deleted during this process. 
    tmpdir = /tmp

    # "big_tmpdir" is for quiver, better in a big disk
    big_tmpdir = /tmp
    
    # various trimming parameters
    min_cov = 8
    max_cov = 64
    trim_align = 50
    trim_plr = 50

    # number of processes used by by blasr during the preassembly process
    q_nproc = 16 

Please change the various ``sge_option_*`` to the proper SGE queue for the SGE
cluster to run the code.

You should estimate the overall coverage and length distribution for putting in
the correct options in the configuration file.  You will need to decide a
length cutoff for the seeding reads. The optimum cutoff length will depend on
the distribution of the sequencing read lengths, the genome size and the
overall yield. The general guideline is the coverage of the seeding sequences
should be above 20x of the genome and the overall coverage should be at least
3x of the coverage of the seeding sequences. Start the Hierarchical Genome
Assembly Process b the assembly process by::

    $ CENT_WF.py HBAR.cfg  

If you want to kill the jobs, you should kill the python process using
``kill`` command and using ``qdel`` for the SGE jobs submitted by the python
process. 

The spec file used by the Celera Assembler is at ``$CENT_HOME/etc/asm.spec``.
In the future, this will be configurable using the configuration file.
