import unittest
from hsmmv2 import *
from main import *

class fpTest(unittest.TestCase):
    def setUp(self) -> None:
        packets = readhexfile('modbus_query.txt')
        self.maxlen=-1
        self.lines=[]
        count=0
        for packet in packets:
            
            if count==3:break
            line=[]
            for num in packet:
                line.append(str(num))
            if len(line)>self.maxlen:
                self.maxlen=len(line)
            self.lines.append(line)
            count+=1
        

    def test_apriori(self):
        self.fp=hsmm(self.lines,10,'frequent_result.txt',self.maxlen)
        self.fp.train()

if __name__ == '__main__':
    suite=unittest.TestSuite()
    suite.addTest(fpTest('test_apriori'))

    runner=unittest.TextTestRunner()
    runner.run(suite)
