import unittest

class ImportTest(unittest.TestCase):
    def test_import_numpy(self):
        import numpy as np
    def test_import_scipy(self):
        import scipy as np
    def test_import_pyplot(self):
        import matplotlib.pyplot as plt
    def test_import_alpspython(self):
        import pyalps
        import pyalps.alea
        import pyalps.plot
    def test_import_os(self):
        from os import mkdir, path, listdir, chdir
