import unittest
from subprocess import Popen
class DeepTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from os import mkdir, chdir,devnull, listdir, path
        self.dir_name='./deep_test/'
        if not path.exists(self.dir_name):
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
    def check_has_observable(self,string, output=False):
        import pyalps
        self.assertGreater(len(pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0]), 0)
    def check_value_scalar(self, string, should_be, significance=8, debug=False):
        import pyalps
        val=pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0][0].y[0].mean
        err=pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0][0].y[0].error
        if debug:
            print('\n')
            print(val)
            print(err)
            print(abs(should_be-val)/max([err,1e-6]))
            print('\n')
        self.assertLess(abs(should_be-val)/max([err,1e-6]),significance)
    def check_value_vector(self, string, index, should_be, significance=8):
        import pyalps
        val=pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0][0].y.mean[index]
        err=pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0][0].y.error[index]
        self.assertLess(abs(should_be-val)/max(err,1e-6),significance)
    def check_no_double_values(self,string, output=False):
        import pyalps
        from numpy import diff, sort
        val=pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'),[string])[0][0].y
        if len(val)>1:
            self.assertGreater(min(diff(sort(abs(val)))),1e-10)

    def test_GS_exchange(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'Initialization' : 'Ferro',
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 100,
                 'SWEEPS'         : 20,
                 'UPDATE'         : "ssf",
                 'L'              : 4,
            }]
        self.should_work(parm)
        self.check_value_scalar('M',1.)
    def test_GS_dipolar_striped(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'Initialization' : 'GS',
                 'D'              : 1,
                 'THERMALIZATION' : 100,
                 'SWEEPS'         : 200,
                 'Spin autocorrelation analysis length':5,
                 'UPDATE'         : "ssf",
                 'L'              : 8,
            }]
        self.should_work(parm)
        self.check_value_vector('Spin Autocorrelation',-1, 1.)
    def test_GS_dipolar_vortex(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'Initialization' : 'Vortex',
                 'D'              : 1,
                 'THERMALIZATION' : 100,
                 'SWEEPS'         : 200,
                 'Spin autocorrelation analysis length':5,
                 'UPDATE'         : "ssf",
                 'L'              : 8,
            }]
        self.should_work(parm)
        self.check_value_vector('Spin Autocorrelation',-1, 1.)
    def test_finite_temperature_ising_exchange(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'Initialization' : 'Ferro',
                 'Ising'          : True,
                 'T'              : 10.,
                 'J'              : 1.,
                 'THERMALIZATION' : 10000,
                 'SWEEPS'         : 50000,
                 'UPDATE'         : "ssf",
                 'L'              : 6,
            }]
        self.should_work(parm)
        self.check_value_scalar('M',0.165) #value from ALPS/spinmc
    #in python 3.4 they added subTest, which would make the following a lot nicer
    #TODO implement this test also for the other interaction sets
    def test_MCRG_no_double_entries_very_small(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 10.,
                 'Initialization' : 'Random',
                 'J'              : 1.,
                 'THERMALIZATION' : 100,
                 'SWEEPS'         : 4,
                 'MCRG Interactions': 'very small',
                 'MCRG Reduction Technique': 'Blockspin',
                 'mcrg_iteration_depth' : 1,
                 'UPDATE'         : "ssf",
                 'L'              : 8
            }]
        self.should_work(parm)
        self.check_has_observable('MCRGe S_alpha0')
        self.check_has_observable('MCRGo S_alpha0')
        self.check_no_double_values('MCRGe S_alpha0')
        self.check_no_double_values('MCRGo S_alpha0')
    def test_MCRG_no_double_entries_small(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 10.,
                 'Initialization' : 'Random',
                 'J'              : 1.,
                 'THERMALIZATION' : 100,
                 'SWEEPS'         : 4,
                 'MCRG Interactions': 'small',
                 'MCRG Reduction Technique': 'Blockspin',
                 'mcrg_iteration_depth' : 1,
                 'UPDATE'         : "ssf",
                 'L'              : 8
            }]
        self.should_work(parm)
        self.check_has_observable('MCRGe S_alpha0')
        self.check_has_observable('MCRGo S_alpha0')
        self.check_no_double_values('MCRGe S_alpha0')
        self.check_no_double_values('MCRGo S_alpha0')
    def test_MCRG_no_double_entries_very_small(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 10.,
                 'Initialization' : 'Random',
                 'J'              : 1.,
                 'THERMALIZATION' : 100,
                 'SWEEPS'         : 4,
                 'MCRG Interactions': 'ising',
                 'MCRG Reduction Technique': 'Blockspin',
                 'mcrg_iteration_depth' : 1,
                 'UPDATE'         : "ssf",
                 'L'              : 8
            }]
        self.should_work(parm)
        self.check_has_observable('MCRGe S_alpha0')
        self.check_has_observable('MCRGo S_alpha0')
        self.check_no_double_values('MCRGe S_alpha0')
        self.check_no_double_values('MCRGo S_alpha0')
