import unittest
from hsmm import *
from main import *

class fpTest(unittest.TestCase):
    def setUp(self) -> None:
        packets = readhexfile('modbus_query.txt')
        self.lines=[]
        for packet in packets:
            for num in packet:
                self.lines.append(str(num))
            break
        

    def test_apriori(self):
        self.fp=hsmm(self.lines,'frequent_result.txt',4)
        for i in range(10):
            self.fp.train(self.fp.cur_seq)

if __name__ == '__main__':
    suite=unittest.TestSuite()
    suite.addTest(fpTest('test_apriori'))

    runner=unittest.TextTestRunner()
    runner.run(suite)
