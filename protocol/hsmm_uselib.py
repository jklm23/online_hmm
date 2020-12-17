import numpy as np
import math as mt
from hsmmlearn.hsmm import MultinomialHSMM

Sigma = [i for i in range(256)]


class hsmm_uselib:
    def __init__(self, sequence, filename, D) -> None:
        super().__init__()
        self.cur_seq = sequence
        self.len = len(sequence)
        self.filename = filename  # 频繁模式结果集
        self.M = len(self.cur_seq)+1  # 状态个数
        self.keys = []
        self.D = D  # 状态的最大持续时间D
        self.readfile()

        
        self.DD=[]  # 状态i持续时间为D的概率
        self.pi_i = []  # 初始状态分布概率
        self.A = []  # 状态转移矩阵
        self.b_i_c_arr=[]   # 状态i观察到c的概率
        
        

        self.initialize()

    def initialize(self):

        # 初始化状态转移矩阵为等概率
        for i in range(self.M):
            tmp = []
            for j in range(self.M):
                tmp.append(1 / self.M)
            self.A.append(tmp)

        self.p_i_d_init()  # 初始化p_i(d)，状态i停留时间为d的概率

        for i in range(self.M):
            # 初始化初始状态分布
            self.pi_i.append(1 / self.M)
        
        # 初始化b_i_c_arr
        for i in range(self.M):
            self.b_i_c_arr.append([])
            for v in range(256):
                self.b_i_c_arr[i].append([])
                self.b_i_c_arr[i][v]=self.b_i_c(i,str(v))
        print('initialize finish')

    def b_i_c(self, i, c):
        if i in range(0, self.N):
            if c in self.keys[i]:
                tmp = self.keys[i].split(',')
                slen = len(tmp)
                # if len(tmp) == 2 and tmp[1] == '':
                #     slen -= 1
                return mt.exp(-slen / 10)
            else:
                return 0
        else:
            return mt.exp(-20)  # -20

    def readfile(self, ):
        file = open(self.filename, 'r')
        lines = file.readlines()
        res=[]
        self.N = len(lines)
        self.frec=0
        for line in lines:
            line = line.split(' ')[0]
            # 每个频繁项集看作一个状态
            self.keys.append(line)
            # for c in line:
            #     if c not in res:
            #         res.append(c)
            #         self.frec+=1

    def b_i_c(self, i, c):
        if i in range(0, self.N):
            if c in self.keys[i]:
                tmp = self.keys[i].split(',')
                slen = len(tmp)
                # if len(tmp) == 2 and tmp[1] == '':
                #     slen -= 1
                # return mt.exp(-slen / 10)
                return 1/slen
            else:
                return 0
        else:
            return mt.exp(-20)  # -20

    def p_i_d_init(self):
        for i in range(self.M):
            self.DD.append([])
            for d in range(0,self.D):
                self.DD[i].append(0)

        allsum = 0
        for i in range(1, self.D + 1):
            allsum += i * i
        for i in range(self.M):
            for d in range(0, self.D):
                dd=d+1
                self.DD[i][d] = (dd * dd) / allsum


    def fit(self,):
        hsmm_learn=MultinomialHSMM(np.asarray(self.b_i_c_arr),np.asarray(self.DD),np.asarray(self.A),np.asarray(self.pi_i))
        numseq=[]
        count=0
        for i in self.cur_seq:
            numseq.append(int(i))
        hsmm_learn.fit(numseq)
        print(hsmm_learn.durations)
        print(hsmm_learn.tmat)
