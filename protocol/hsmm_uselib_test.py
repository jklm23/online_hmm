import unittest
from hsmm_uselib import *
from main import *

class fpTest(unittest.TestCase):
    def setUp(self) -> None:
        packets = readhexfile('modbus_query.txt')
        self.lines=[]
        for packet in packets:
            for num in packet:
                self.lines+=str(num)
            break

    def test_apriori(self):
        self.fp=hsmm_uselib(self.lines,'frequent_result.txt',4)
        self.fp.fit()

if __name__ == '__main__':
    suite=unittest.TestSuite()
    suite.addTest(fpTest('test_apriori'))

    runner=unittest.TextTestRunner()
    runner.run(suite)
