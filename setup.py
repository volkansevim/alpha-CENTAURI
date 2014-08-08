from distutils.core import setup

setup(
        name = 'alpha-CENTAURI',
        version = '0.1.0',
        author = 'Volkan Sevim & Jason Chin',
        author_email = 'vsevim@pacificbiosciences.com',
        packages = ['alpha-CENTAURI'],
	package_dir = {'alpha-CENTAURI': 'src'},
        scripts = ['src/chop_to_monomers.py',  
		   'src/count_found_monomers_in_fa.py', 
		   'src/monomer_graph_analysis.py' ],
	
        url = 'https://github.com/volkansevim/alpha-CENTAURI',
        download_url = 'https://github.com/volkansevim/alpha-CENTAURI/',
	license = 'LICENSE.txt',
	install_requires = [ "pbcore >= 0.6.3", "falcon >= 0.1", "networkx >= 1.7", "numpy >= 1.7" ],
	include_package_data = True,
	package_data = {'': ['examples/*'] }
)
