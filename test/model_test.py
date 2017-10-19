import unittest
from subprocess import Popen
class ModelTest(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        from os import mkdir, chdir,devnull, getcwd
        self.dir_name='./model_test/'
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

    def test_not_enough_sweeps(self):
        parm=[{
                 'LATTICE'        : "square lattice",
                 'T'              : 0.,
                 'J'              : 1,
                 'THERMALIZATION' : 1,
                 'SWEEPS'         : 3,
                 "ALGORITHM"      : "xy",
                 'basic_observables': False,
                 'L'              : 4,
                 'Each_Measurement' : 5
            }]
        self.should_crash(parm)
