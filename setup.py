from distutils.core import setup

setup(
        name = 'alpha-CENTAURI',
        version = '0.1.0',
        author = 'Volkan Sevim & Jason Chin',
        author_email = 'vsevim@pacificbiosciences.com',
        packages = ['alpha-CENTAURI'],
	package_dir = {'alpha-CENTAURI': 'src'},
        scripts = ['chop_to_monomers.py',  
		   'count_found_monomers_in_fa.py', 
		   'monomer_graph_analysis.py' ],
	i
        url = 'https://github.com/volkansevim/alpha-CENTAURI',
        download_url = 'https://github.com/volkansevim/alpha-CENTAURI/',
	license = 'LICENSE.txt',
)
