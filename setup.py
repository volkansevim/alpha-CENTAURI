from distutils.core import setup

setup(
        name = 'alpha-CENTAURI',
        version = '0.1.0',
        author = 'Volkan Sevim & Jason Chin',
        author_email = 'vsevim@pacificbiosciences.com',
        packages = ['alpha-CENTAURI'],
	package_dir = {'alpha-CENTAURI': 'src'},
        scripts = [ 'calc_alignment_scores_for_inferred_monomers.py', 
		   'chop_to_monomers.py', 
		   'cluster_monomers.py', 
		   'count_found_monomers_in_fa.py', 
		   'get_reads_with_monomers.py' ],
	
        url = 'https://github.com/volkansevim/alpha-CENTAURI',
        download_url = 'https://github.com/volkansevim/alpha-CENTAURI/',
	license = 'LICENSE.txt',
)
