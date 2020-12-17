import numpy as np
import math

class hsmm:
    def __init__(self, sequences, M, filename, seqmax):
        self.seqmax = seqmax
        self.keys = []
        self.filename = filename
        self.klMax = -1  # key的最大长度
        self.readKeys()
        self.sequences = sequences  # 报文序列
        self.br_list = []  # 各个报文的起始位置
        self.br_len = []  # 报文的各个长度
        self.T = 0  # 总报文长度
        self.M = M  # 设定状态个数
        # self.Dmax=Dmax #状态的最大持续时间（包括data部分）
        self.seq = ['']  # 报文连接在一起，添加一个占位符

        for i in range(len(self.sequences)):
            self.br_list.append(self.T+1)  # 记录报文的开始位置
            for j in self.sequences[i]:
                self.seq.append(j)
                self.T += 1

        for i in range(len(self.sequences)-1):
            self.br_len.append(self.br_list[i+1]-self.br_list[i])
        self.br_len.append(self.T-self.br_list[-1]+1)

        # HSMM四个参数的初始化
        self.A = np.full((M, M), 1/M, dtype=float)
        self.pi = np.full(M, 1/M, dtype=float)
        self.L_arr = np.full((M, len(self.keys)+1, self.seqmax+1),
                             1/self.seqmax, dtype=float)    # [:,:,0]空出来
        self.k_j_key = np.full((M, len(self.keys)+1),
                               1/(len(self.keys)+1), dtype=float)  # 还要考虑非key序列
        self.ct = []
        # self.b_j_d_o_arr=np.zeros((self.M,self.Dmax,self.T,self.T))
       
        # self.yita_arr={}

        self.initialize()

    def readKeys(self, ):
        file = open(self.filename, 'r')
        lines = file.readlines()
        self.N = len(lines)
        for line in lines:
            line = line.split(' ')[0]
            if len(line) > self.klMax:
                self.klMax = len(line.split(','))
            # 每个频繁项集看作一个状态
            self.keys.append(line)

    def find_key_index(self, key):
        i = 0
        while i < len(self.keys):
            if self.keys[i].split(',') == key:
                break
            i += 1
        return i

    def l_j_key(self, j, key, d):
        # print(d)
        i = self.find_key_index(key)
        if i < len(self.keys):
            print('找到'+key+'位置：%d'%i)
            return self.L_arr[j, i, d]
        else:
            return self.L_arr[j, -1, d]

    def k_j(self, j, o):
        i = self.find_key_index(o)
        if i < len(self.keys):
            return self.k_j_key[j, i]
        else:
            return self.k_j_key[j, -1]  # non-key

    def b_j_d_o(self, j, d, start):
        # print(d)
        res = 0
        for kl in range(1, min(d, self.klMax)+1):
            res += self.k_j(j, self.seq[start:start+kl]) * \
                self.l_j_key(j, self.seq[start:start+kl], d)
        return res

    def alpha(self, t, j):

        if self.alpha_t_j_arr[t, j] > -1:
            return self.alpha_t_j_arr[t, j]
        if t == 0:
            self.alpha_t_j_arr[t, j] = self.pi[j]
            return self.pi[j]
        startindex = 0
        while startindex < len(self.br_list) and t >= self.br_list[startindex]:
            startindex += 1

        # print('%d %d'%(t,self.br_list[startindex-1]))
        # if startindex>0:
        Dt = t-self.br_list[startindex-1]+1  # startindex-1为当前报文

        res = 0
        for i in range(self.M):
            for d in range(1, Dt+1):
                res += self.alpha(t-d, i) * \
                    self.A[i, j]*self.b_j_d_o(j, d, t-d+1)
        self.alpha_t_j_arr[t, j] = res
        # print('t:%d j:%d Dt:%d res:%lf'%(t,j,Dt,res))
        return res

    def Lkh_calc(self):
        for j in range(self.M):
            self.Lkh += self.alpha(self.T-1, j)

    def Et_(self, t):
        # startindex指向下一个报文
        if t==0:
            return self.br_len[0]
        startindex = 0
        while startindex < len(self.br_list) and t >= self.br_list[startindex]:
            startindex += 1

        if t < self.br_list[startindex-1]+self.br_len[startindex-1]-1:
            # 如果t不在这个报文的末尾
            Et = self.br_list[startindex-1]+self.br_len[startindex-1]-1-t
        else:
            # 如果t在这个报文的末尾
            # print('startindex:%d,len(self.br_list):%d'%(startindex,len(self.br_list)))
            if startindex-1 < len(self.br_list)-1:
                # 如果不是最后一个报文
                Et = self.br_len[startindex]
            else:
                # 如果是最后一个报文
                Et = 0
        # print('t:%d,Et:%d'%(t,Et))
        return Et

    def beta(self, t, i):
        if self.beta_t_i_arr[t, i] > -1:
            return self.beta_t_i_arr[t, i]
        if t == self.T:
            self.beta_t_i_arr[t, i] = 1
            return 1
        hEt = self.Et_(t)
        # print('Et:%d' % hEt)
        res = 0
        for j in range(self.M):
            for d in range(1, hEt+1):
                res += self.A[i, j]*self.b_j_d_o(j, d, t+1)*self.beta(t+d, j)
        self.beta_t_i_arr[t, i] = res
        # print('t:%d i:%d Et:%d res:%lf'%(t,i,Et,res))
        return res

    # 定义中间变量
    def kesi_t_i_j(self, t, i, j):
        if self.kesi_arr[t, i, j] > -1:
            return self.kesi_arr[t, i, j]
        res = 0
        for d in range(1, self.Et_(t)+1):
            res += self.A[i, j]*self.b_j_d_o(j, d, t+1)*self.beta(t+d, j)
        self.kesi_arr[t, i, j] = res*self.alpha(t, i)
        return self.kesi_arr[t, i, j]

    def pisi_(self, t, j, len):
        if self.psi_arr[t, j, len] > -1:
            return self.psi_arr[t, j, len]
        res = 0
        for i in range(self.M):
            for d in range(len, self.Et_(t)+1):
                res += self.alpha(t, i)*self.A[i, j]*self.k_j(j, self.seq[t+1:t+len+1])*self.l_j_key(
                    j, self.seq[t+1:t+len+1], d)*self.beta(t+d, j)
        self.psi_arr[t, j, len] = res
        return res

    def yita_(self, t, j, key, d):
        res = 0
        for i in range(self.M):
            if self.seq[t-d+1:t-d+len(self.keys[key].split(','))+1] == self.keys[key].split(','):
                res += self.alpha(t-d, i) * \
                    self.A[i, j]*self.k_j(j, key)*self.beta(t, j)
        return res

    def train(self):
        self.initialize()
        A=self.A.copy()
        k_j_key=self.k_j_key.copy()
        l_j_key=self.L_arr.copy()
        pi_i=self.pi.copy()

        for i in range(self.M):
            for j in range(self.M):
                tmp=0
                for t in range(0,self.T):
                    tmp+=self.kesi_t_i_j(t,i,j)
                A[i,j]=tmp*(1/self.Lkh)
                print('finish A(%d,%d)'%(i,j))

        for j in range(self.M):
            for key in range(len(self.keys)):
                tmp=0
                keystr=self.keys[key].split(',')
                for t in range(0,self.T-len(keystr)+1):
                    if self.seq[t+1:t+len(keystr)+1] == keystr or key==len(self.keys)-1:
                        tmp+=self.pisi_(t,j,len(keystr))
                k_j_key[j,key]=tmp*(1/self.Lkh)
                print('finish k(%d,%d)'%(j,key))

                for d in range(1,self.seqmax+1):
                    tmp=0
                    for t in range(d,self.T+1):
                        tmp+=self.yita_(t,j,key,d)
                    l_j_key[j,key,d]=tmp*(1/self.Lkh)
                    print('finish l(%d,%d,%d)'%(j,key,d))

        for i in range(self.M):
            tmp=0
            for j in range(self.M):
                tmp+=self.kesi_t_i_j(0,i,j)
            pi_i[i]=(1/self.Lkh)*tmp
            print('finish pi[%d]'%i)
        
        self.A=A
        self.pi=pi_i
        self.k_j_key=k_j_key
        self.L_arr=l_j_key
        self.Lkh_calc()
        print('train finish')


    def initialize(self):
        # initialize l_j_key_d
        self.alpha_t_j_arr = np.full((self.T+1, self.M), -1, dtype=float)
        self.Lkh = 0
        self.beta_t_i_arr = np.full((self.T+1, self.M), -1, dtype=float)

        # 中间变量
        self.kesi_arr = np.full((self.T+1, self.M, self.M), -1, dtype=float)
        self.psi_arr = np.full((self.T+1, self.M, self.T+1), -1, dtype=float)

        for j in range(self.M):
            for key in range(len(self.keys)):
                tmp=np.zeros((self.seqmax+1),dtype=float)
                keystr=self.keys[key]
                d=0
                while d!=len(keystr.split(',')):
                    d+=1
                tmp[d:]=1/(self.seqmax-d+1)
                self.L_arr[j,key]=tmp




        self.Lkh_calc()

        # #计算比例因子
        # for t in range(self.T+1):
        #     res=0
        #     for i in range(self.M):
        #         res+=self.alpha(t,i)
        #     self.ct.append(1/res)

        for t in range(1,self.T+1):
            for i in range(self.M):
                self.alpha(t, i)
                self.beta(t, i)
        
