import unittest
import options
class CompileTest(unittest.TestCase): #TODO setup the class such that it works with a pure git commit
    @classmethod
    def setUpClass(self):
        from os import mkdir, chdir,devnull, getcwd
        self.old_dir=getcwd()
        chdir(options.MCPP_TEST_OPTIONS_BUILD_DIR)
        self.oblivion=open(devnull, 'w')
    @classmethod
    def tearDownClass(self):
        from os import chdir
        chdir(self.old_dir)
    def test_cmake(self):
        from subprocess import call
        OUT=call('cmake ../', shell=True, stdout=self.oblivion)
        self.assertEqual(OUT,0)
    def test_make(self):
        from subprocess import call
        OUT=call('make -j', shell=True, stdout=self.oblivion)
        self.assertEqual(OUT,0)
