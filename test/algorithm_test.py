import unittest
from subprocess import Popen
class AlgorithmTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from os import mkdir, chdir,devnull, getcwd
        self.dir_name='./algorithm_test/'
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

    def test_algorithm_xy(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "xy",
                 'basic_observables': False,
                 'L'              : 4
            }]
        self.should_work(parm)
    def test_algorithm_exmc(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T_MIN'          : 0.1,
                 'T_MAX'          : 0.2,
                 'NUM_REPLICAS'   : 2,
                 'OPTIMIZATION_ITERATIONS':2,
                 'OPTIMIZE_TEMPERATURE': True,
                 'J'              : 1,
                 'THERMALIZATION' : 1000,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "exmc",
                 'basic_observables': False,
                 'L'              : 4
            }]
        self.should_work(parm)
    def test_algorithm_not_detected(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 200,
                 "ALGORITHM"      : "something else",
                 'basic_observables': False,
                 'L'              : 4
            }]
        self.should_crash(parm)
