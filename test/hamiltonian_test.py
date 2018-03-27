import unittest
from subprocess import Popen
class HamiltonianTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from os import mkdir, chdir,devnull, getcwd
        self.dir_name='./hamiltonian_test/'
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

    def test_no_hamiltonian(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'basic_observables': False,
                 'L'              : 6,
            }]
        self.should_crash(parm)
    def test_NN_run(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 10,
                 'basic_observables': False,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'L'              : 6,
            }]
        self.should_work(parm)
    def test_mpi_run(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1.,
                 'THERMALIZATION' : 10,
                 'basic_observables': False,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'L'              : 6,
            }]
        self.should_work(parm)
    def test_pure_dipolar_w_honeycomb_lattice(self):
        parm=[{
                 'LATTICE'        : "honeycomb lattice",
                 'basic_observables': False,
                 'T'              : 0.,
                 'D'              : 1.,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'cutoff_distance': 3.,
                 'L'              : 10,
            }]
        self.should_work(parm)
    def test_pure_dipolar_run(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'basic_observables': False,
                 'T'              : 0.,
                 'D'              : 1.,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'cutoff_distance': 3.,
                 'L'              : 10,
            }]
        self.should_work(parm)
    def test_dipolar_crash_for_too_small_cutoff(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'basic_observables': False,
                 'T'              : 0.,
                 'D'              : 1.,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'cutoff_distance': 3.,
                 'L'              : 4,
            }]
        self.should_crash(parm)
    def test_positional_disorder_dipolar_run(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'basic_observables': False,
                 'T'              : 0.,
                 'D'              : 1.,
                 'Positional Disorder': 0.5,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'cutoff_distance': 3.,
                 'L'              : 10,
            }]
        self.should_work(parm)
    def test_diluted_dipolar_run(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'basic_observables': False,
                 'T'              : 0.,
                 'D'              : 1.,
                 'Dilution Rate'  : 0.5,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'cutoff_distance': 3.,
                 'L'              : 10,
            }]
        self.should_work(parm)
    def test_shape_anisotropy_run(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'Shape Anisotropy Strength' : 1.,
                 'basic_observables': False,
                 'Shape Anisotropy p' : 2,
                 'THERMALIZATION' : 10,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'cutoff_distance': 3.,
                 'L'              : 4,
            }]
        self.should_work(parm)
