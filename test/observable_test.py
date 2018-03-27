import unittest
from subprocess import Popen
class ObservableTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from os import mkdir, chdir,devnull
        self.dir_name='./observable_test/'
        from os.path import exists
        if exists(self.dir_name):
            rmtree(self.dir_name)
        mkdir(self.dir_name)
        chdir(self.dir_name)
        self.oblivion=open(devnull, 'w')
    @classmethod
    def tearDownClass(self):
        from os import chdir
        from shutil import rmtree
        chdir('..')
        rmtree(self.dir_name)
    def tearDown(self):
        import subprocess
        subprocess.call('rm parm*',shell=True)

    def call(self,parm, output=False):
        import pyalps
        import subprocess
        input_file = pyalps.writeInputFiles('parm',parm)
        if not output:
            return subprocess.call('mc++ parm.in.xml', stdout=self.oblivion, stderr=subprocess.STDOUT, shell=True)
        else:
            return subprocess.call('mc++ parm.in.xml', shell=True)
    def should_crash(self,parm,output=False):
        self.assertNotEqual(self.call(parm,output),0)
    def should_work(self,parm, output=False):
        self.assertEqual(self.call(parm,output),0)
    def check_data_length(self,string,should_be, output=False):
        import pyalps
        l=len(pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0][0].y)
        self.assertEqual(l,should_be)
    def check_has_observable(self,string, wait=False):
        import pyalps
        if wait:
            a=input('wait')
        self.assertGreater(len(pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0]), 0)

    def test_basic_observables(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'L'              : 4,
            }]
        self.should_work(parm)
        self.check_has_observable('M')
        self.check_has_observable('M^2')
        self.check_has_observable('M^4')
        self.check_has_observable('M staggered')
        self.check_has_observable('M staggered^2')
        self.check_has_observable('M staggered^4')
    def test_component_observables(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 20,
                 'component_observables' : True,
                 "ALGORITHM"      : "xy",
                 'L'              : 4,
            }]
        self.should_work(parm)
        self.check_has_observable('M striped')
        self.check_has_observable('M microvortex')
    def test_last_configuration(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'measure last configuration' : True,
                 'L'              : 4,
            }]
        self.should_work(parm)
        self.check_has_observable('Last Configuration')
        self.check_data_length('Last Configuration', 4*4)
    def test_structure_factor(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'structure_factor' : True,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'L'              : 4,
            }]
        self.should_work(parm)
        self.check_has_observable('|Structure Factor|^2')
        self.check_data_length('|Structure Factor|^2', 4*4)
    def test_MCRG_nothing_loaded(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 1000,
                 'mcrg_iteration_depth' : 2,
                 "ALGORITHM"      : "xy",
                 'L'              : 32,
            }]
        self.should_crash(parm)
    def test_MCRG_no_interactions(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 1000,
                 'mcrg_iteration_depth' : 2,
                 'MCRG Reduction Technique': 'Blockspin',
                 "ALGORITHM"      : "xy",
                 'L'              : 32,
            }]
        self.should_crash(parm)
    def test_MCRG_no_reduction(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 1000,
                 'mcrg_iteration_depth' : 2,
                 'MCRG Interactions': 'small',
                 "ALGORITHM"      : "xy",
                 'L'              : 32,
            }]
        self.should_crash(parm)
    def test_MCRG_working(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 1.,
                 'J'              : 1.,
                 'Initialization' : "Random",
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 1000,
                 'MCRG Interactions': 'small',
                 'MCRG Reduction Technique': 'Blockspin',
                 'mcrg_iteration_depth' : 1,
                 "ALGORITHM"      : "xy",
                 'L'              : 32,
            }]
        self.should_work(parm)
        self.check_has_observable('MCRGe S_alpha0')
        self.check_has_observable('MCRGe S_alpha0 S_beta0')
        self.check_has_observable('MCRGe S_alpha0 S_beta1')
        self.check_has_observable('MCRGe S_alpha1')
        self.check_has_observable('MCRGe S_alpha1 S_beta1')
        self.check_has_observable('MCRGo S_alpha0')
        self.check_has_observable('MCRGo S_alpha0 S_beta0')
        self.check_has_observable('MCRGo S_alpha0 S_beta1')
        self.check_has_observable('MCRGo S_alpha1')
        self.check_has_observable('MCRGo S_alpha1 S_beta1')
    def test_LLG_basic_observables(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 'llg'            : True,
                 'LLG measures per measure' : 5,
                 "ALGORITHM"      : "xy",
                 'cutoff_distance': 3.,
                 'L'              : 10,
            }]
        self.should_work(parm)
        self.check_has_observable('LLG M')
        self.check_has_observable('LLG Ms')
        self.check_has_observable('LLG E')
    def test_LLG_Muon(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 'llg'            : True,
                 'LLG measures per measure' : 5,
                 'LLG Measure Muon': True,
                 "ALGORITHM"      : "xy",
                 'cutoff_distance': 3.,
                 'L'              : 10,
            }]
        self.should_work(parm)
        self.check_has_observable('LLG M')
        self.check_has_observable('LLG Muon Depolarization')
    def test_spin_autocorrelation(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'Spin autocorrelation analysis length' : 5,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'L'              : 4,
            }]
        self.should_work(parm)
        self.check_has_observable('Spin Autocorrelation')
    def test_field_histogram(self):
        nbins=64
        parm=[{
		 'LATTICE'        : "square lattice",
		 'T'              : 0.,
		 'D'              : 1.,
		 'THERMALIZATION' : 10,
		 'SWEEPS'         : 200,
		 'Field Histogram': True,
		 'Field Histogram fixed z': False,
		 'Field Histogram diameter split': 3,
	         'Field Histogram n_bins' : nbins,
                 "ALGORITHM"      : "xy",
		 'cutoff_distance': 3.,
		 'L'              : 8,
            }]
        self.should_work(parm)
        self.check_has_observable('Field Histogram Absolute Value')
        self.check_has_observable('Field Histogram xy'            )
        self.check_has_observable('Field Histogram z'             )
        self.check_data_length('Field Histogram Absolute Value', nbins)
        self.check_data_length('Field Histogram xy'            , nbins)
        self.check_data_length('Field Histogram z'             , nbins)
    def test_field_histogram_log(self):
        nbins=64
        parm=[{
		 'LATTICE'        : "square lattice",
		 'T'              : 0.,
		 'D'              : 1.,
		 'THERMALIZATION' : 10,
		 'SWEEPS'         : 20,
		 'Field Histogram': True,
		 'Field Histogram Log Scale' : True,
                 'Field Histogram fixed z': False,
		 'Field Histogram diameter split': 3,
	         'Field Histogram n_bins' : nbins,
                 "ALGORITHM"      : "xy",
		 'cutoff_distance': 3.,
		 'L'              : 8,
            }]
        self.should_work(parm)
        self.check_has_observable('Field Histogram Absolute Value')
        self.check_has_observable('Field Histogram xy'            )
        self.check_has_observable('Field Histogram z'             )
        self.check_data_length('Field Histogram Absolute Value', nbins)
        self.check_data_length('Field Histogram xy'            , nbins)
        self.check_data_length('Field Histogram z'             , nbins)
