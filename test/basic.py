#!/opt/alps/bin/alpspython
import results
import algorithm_test
import compile_test
import hamiltonian_test
import import_test
import model_test
import observable_test
import unittest
from sys import exit

res = results.my_test_results()
suite = unittest.TestSuite()
suite.addTest(unittest.makeSuite(compile_test.CompileTest))
suite.addTest(unittest.makeSuite(algorithm_test.AlgorithmTest))
suite.addTest(unittest.makeSuite(model_test.ModelTest))
suite.addTest(unittest.makeSuite(import_test.ImportTest))
suite.addTest(unittest.makeSuite(hamiltonian_test.HamiltonianTest))
suite.addTest(unittest.makeSuite(observable_test.ObservableTest))
suite.run(res)
print(res.get_string())
exit(res.get_number_errors())
