from __future__ import print_function
import unittest
import re
class my_test_results(unittest.TestResult):
    def __init__(self):
        super(my_test_results, self).__init__(self)
        self.Suite =list()
        self.Test  =list()
        self.Result=list()
        self.success_msg='SUCCESS'
        self.fail_msg   ='FAIL'
        self.err_msg    ='ERROR'
    def get_delimiter(self):
        return '-'*65+'\n'
    def get_summary(self):
        n_err =len([x for x in self.Result if x == self.err_msg])
        n_suc =len([x for x in self.Result if x == self.success_msg])
        n_fail=len([x for x in self.Result if x == self.fail_msg])
        n_tests=len(self.Suite)
        overview_string=''
        if(n_err==0 and n_fail==0 and n_suc==n_tests):
            overview_string+='ALL TESTS OK\n'
        overview_string+=self.success_msg.ljust(15)+str(n_suc ).rjust(4)+'/'+str(n_tests)+'\n' 
        overview_string+=self.err_msg    .ljust(15)+str(n_err ).rjust(4)+'/'+str(n_tests)+'\n'
        overview_string+=self.fail_msg   .ljust(15)+str(n_fail).rjust(4)+'/'+str(n_tests)+'\n'
        return overview_string
    def get_string(self):
        len_suite  =max([len(x) for x in self.Suite ])+1
        len_test   =max([len(x) for x in self.Test  ])+2
        len_result =max([len(x) for x in self.Result])
        z=zip(self.Suite,self.Test,self.Result)
        delimiter=self.get_delimiter()
        overview_string=self.get_summary()
        s='\n'+delimiter+overview_string+delimiter
        for r in z:
            s+=r[0].ljust(len_suite)+': '+r[1].replace('_',' ').ljust(len_test)+'['+r[2].center(len_result)+']\n'
        s+=delimiter+overview_string+delimiter
        return s
    def analyze_output(self, test, result):
        s=str(test)
        import string
        m=re.match(r'test_(\w+) \(\w+\.(\w+)\)',s)
        self.Suite .append(m.group(2))
        self.Test  .append(m.group(1))
        self.Result.append(result)
    def addSuccess(self, test):
        import sys
        super(my_test_results, self).addSuccess(test) 
        print('.',end='')
        sys.stdout.flush()
        self.analyze_output(test,self.success_msg)
    def addError(self, test, err):
        import sys
        print('E',end='')
        sys.stdout.flush()
        super(my_test_results, self).addError(test, err) 
        self.analyze_output(test,self.err_msg)
    def addFailure(self, test, err):
        import sys
        print('F',end='')
        sys.stdout.flush()
        super(my_test_results, self).addFailure(test, err) 
        self.analyze_output(test,self.fail_msg)
    def get_number_errors(self):
        n_err =len([x for x in self.Result if x == self.err_msg])
        n_fail=len([x for x in self.Result if x == self.fail_msg])
        return n_err+n_fail 
