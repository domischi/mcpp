import unittest
from subprocess import Popen
class ObservableTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from os import mkdir, chdir,devnull
        self.dir_name='./observable_test/'
        mkdir(self.dir_name)
        chdir(self.dir_name)
        self.oblivion=open(devnull, 'w')
    @classmethod
    def tearDownClass(self):
        from os import chdir
        from shutil import rmtree
        chdir('..')
        rmtree(self.dir_name)

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
    def check_data_length(self,string,should_be,error_msg, output=False, EQUAL=True, GREATER=False):
        import pyalps
        l=len(pyalps.collectXY(pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string]),x='T',y=string)[0].y)
        if EQUAL:
            self.assertEqual(l,should_be,error_msg)
        if GREATER:
            self.assertGreater(l,should_be,error_msg)
    def check_has_observable(self,string,error_msg, output=False):
        import pyalps
        self.assertGreater(len(pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0]), 0, error_msg)

    def test_basic_observables(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 2,
                 'UPDATE'         : "ssf",
                 'L'              : 4,
            }]
        self.should_work(parm)
    def test_last_configuration(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 2,
                 'UPDATE'         : "ssf",
                 'measure last configuration' : True,
                 'L'              : 4,
            }]
        self.should_work(parm)
        self.check_data_length('Last Configuration', 4*4, 'mc++ returned a wrong last configuration (data length)')
    def test_structure_factor(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'structure_factor' : True,
                 'SWEEPS'         : 2,
                 'UPDATE'         : "ssf",
                 'L'              : 4,
            }]
        self.should_work(parm)
        self.check_has_observable('|Structure Factor|^2','mc++ returnd a wrong structure factor (data length)')
        self.check_has_observable('|Structure Factor|^2','mc++ returnd a wrong structure factor (data length)')
    def test_MCRG_nothing_loaded(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 2,
                 'mcrg_iteration_depth' : 2,
                 'UPDATE'         : "ssf",
                 'L'              : 32,
            }]
        self.should_crash(parm)
    def test_MCRG_no_interactions(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 2,
                 'mcrg_iteration_depth' : 2,
                 'MCRG Reduction Technique': 'Blockspin',
                 'UPDATE'         : "ssf",
                 'L'              : 32,
            }]
        self.should_crash(parm)
    def test_MCRG_no_reduction(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 2,
                 'mcrg_iteration_depth' : 2,
                 'MCRG Interactions': 'small',
                 'UPDATE'         : "ssf",
                 'L'              : 32,
            }]
        self.should_crash(parm)
    def test_MCRG_working(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 2,
                 'MCRG Interactions': 'small',
                 'MCRG Reduction Technique': 'Blockspin',
                 'mcrg_iteration_depth' : 2,
                 'UPDATE'         : "ssf",
                 'L'              : 32,
            }]
        self.should_work(parm)
        self.check_data_length('MCRGe S_alpha0',0,'mc++ returnd a wrong MCRG matrix (data length)', EQUAL=False, GREATER=True)
        self.check_data_length('MCRGo S_alpha0',0,'mc++ returnd a wrong MCRG matrix (data length)', EQUAL=False, GREATER=True)
        self.check_data_length('MCRGe S_alpha1',0,'mc++ returnd a wrong MCRG matrix (data length)', EQUAL=False, GREATER=True)
    def test_LLG_basic_observables(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 2,
                 'llg'            : True,
                 'UPDATE'         : "ssf",
                 'cutoff_distance': 3.,
                 'L'              : 10,
            }]
        self.should_work(parm)
    def test_LLG_Muon(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 2,
                 'llg'            : True,
                 'LLG Measure Muon': True,
                 'UPDATE'         : "ssf",
                 'cutoff_distance': 3.,
                 'L'              : 10,
            }]
        self.should_work(parm)
    def test_spin_autocorrelation(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'Positional Disorder': 0.5,
                 'Spin autocorrelation analysis length' : 5,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 2,
                 'UPDATE'         : "ssf",
                 'cutoff_distance': 3.,
                 'L'              : 4,
            }]
        self.should_work(parm)
