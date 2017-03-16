#!/opt/alps/bin/alpspython
import results
import hamiltonian_test
import import_test
import compile_test
import observable_test
import deep_test
import unittest
from sys import exit

res = results.my_test_results()
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(compile_test.CompileTest))
suite.addTest(unittest.makeSuite(import_test.ImportTest))
suite.addTest(unittest.makeSuite(hamiltonian_test.HamiltonianTest))
suite.addTest(unittest.makeSuite(observable_test.ObservableTest))
suite.addTest(unittest.makeSuite(deep_test.DeepTest))
suite.run(res)
print(res.get_string())
exit(res.get_number_errors())
