import numpy as np
import math as mt

Sigma = [i for i in range(256)]


class hsmm:
    def __init__(self, sequence, filename, D):
        # self.sequences = sequences  # 观测序列
        self.cur_seq=sequence   #当前的观测序列
        # self.len = len(sequences)
        self.filename = filename

        self.keys = []  # 频繁项集
        self.keys_int=[]
        
        
        
        self.readfile()
        self.initial_keys()
        self.N = len(self.keys)  # 频繁模式个数

        self.D = D  # 状态的最大持续时间D
        self.DD = {}  # {状态i,停留时间d:概率}
        self.A = []  # 状态转移矩阵
        self.M = self.N+1  # 状态个数
        self.pi_i = []  # 初始状态分布概率
        self.b_i_c_arr=[]

        self.initialize()  # ininitialize
        self.init_forward_backward()

    # def set_cur_seq(self,seq):


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
                self.b_i_c_arr[i].append(self.b_i_c(i,str(v)))

    def init_forward_backward(self):
        
        self.alpha = []  # forward procedure,the length is max(len(message))::
        self.beta = []  # backward procedure
        self.gamma = []  # {(t,i}:possibility}
        self.yita_arr={}

        for i in range(len(self.cur_seq)):
            self.alpha.append({})
            self.beta.append({})

        for i in range(len(self.cur_seq)):
            self.gamma.append([])

        

        # 初始化前向和后向矩阵
        # TODO: underflow_arr=[] #防止下溢，在beta上乘以一个比例因子
        for t in range(len(self.cur_seq)):
            for i in range(self.M):
                for d in range(1,self.D+1):
                    self.alpha[t][(i,d)]=self.__forward(t,i,d)

        for i in range(len(self.cur_seq)):
            self.gamma.append([])
            for j in range(self.M):
                self.gamma[i].append(0)

        for t in range(len(self.cur_seq)-1,-1,-1):
            for i in range(self.M):
                self.gamma[t][i]=self.__gamma_t(t,i)
                for d in range(1,self.D+1):
                    self.beta[t][(i,d)]=self.__backward(t,i,d)

        for t in range(len(self.cur_seq)):
            for i in range(self.M):
                for d in range(1,self.D+1):
                    self.yita_arr[(t,i,d)]=self.yita(t,i,d)
        
        
        print('initialize finish')
        

    def kesi(self, t, i, j):
        # TODO:t=0?
        if t==0:return 1
        tmp = 0
        for d in range(1, self.D+1):
            tmp += self.DD[(j, d)]*self.beta[t][(j, d)]
        return self.alpha[t-1][(i, 1)]*self.A[i][j]*self.b_i_c_arr[j][int(self.cur_seq[t])]*tmp

    def yita(self, t, i, d):
        # if i==self.M-1:
        #     if 
            
        tmp = 0
        for j in range(self.M):
            if j == i:
                continue
            tmp += self.alpha[t-1][(j, 1)]*self.A[j][i]
        return tmp*self.b_i_c_arr[i][int(self.cur_seq[t])]*self.DD[(i, d)]*self.beta[t][(i, d)]


    def p_i_d_init(self):

        allsum = 0
        for i in range(1, self.D + 1):
            allsum += i * i
        for i in range(self.M):
            for d in range(1, self.D + 1):
                self.DD[(i, d)] = (d * d) / allsum
    # def k_j_key_init(self):
    #     self.b_i_c={}
    #     for i in self.keys.keys():
    #         self.k_j_key[i]=1/self.keylen

    def readfile(self, ):
        file = open(self.filename, 'r')
        lines = file.readlines()

        self.N = len(lines)
        for line in lines:
            line = line.split(' ')[0]
            # 每个频繁项集看作一个状态
            self.keys.append(line)

    def train(self,o):
        self.cur_seq=o
        tmpA=self.A.copy()
        tmpPi=self.pi_i.copy()
        tmpb_i_c=self.b_i_c_arr.copy()
        tmpDD=self.DD.copy()

        for i in range(self.M):
            sigma_gamma=0
            for j in range(self.M):
                sigma_gamma+=self.gamma[0][j]
            tmpPi[i]=self.gamma[0][i]/sigma_gamma
        
            for j in range(self.M):
                fenzi=0
                for t in range(0,len(self.cur_seq)):
                    fenzi+=self.kesi(t,i,j)
                fenmu=0
                for k in range(self.M):
                    if k==i:continue
                    for t in range(0,len(self.cur_seq)):
                        fenmu+=self.kesi(t,i,k)


                tmpA[i][j]=fenzi/fenmu
            
            for d in range(1,self.D+1):
                fenzi=0
                for t in range(len(self.cur_seq)):
                    fenzi+=self.yita_arr[(t,i,d)]
                fenmu=0
                for dd in range(1,self.D+1):
                    for t in range(len(self.cur_seq)):
                        fenmu+=self.yita_arr[(t,i,dd)]

                tmpDD[(i,d)]=fenzi/fenmu

            for v in range(256):
                fenzi=0
                for t in range(len(self.cur_seq)):
                    if self.cur_seq[t]==str(v):
                        continue
                    fenzi+=self.gamma[t][i]
                fenmu=0
                for k in range(256):
                    
                    for t in range(len(self.cur_seq)):
                        if self.cur_seq[t]==str(k):continue
                        fenmu+=self.gamma[t][i]
                tmpb_i_c[i][v]=fenzi/fenmu
            print('%d finish'%i)
        print('train finish')
        self.A=tmpA
        self.pi_i=tmpPi
        self.DD=tmpDD
        self.b_i_c_arr=tmpb_i_c
        self.init_forward_backward()

        T=len(self.cur_seq)-1
        result=[]
        while T>0:
            maxres=-1
            maxj=maxd=0
            for j in range(self.M):
                for d in range(1,self.D):
                    if self.yita_arr[(T,j,d)]>maxres:
                        maxres=self.yita_arr[(T,j,d)]
                        maxj=j
                        maxd = d
            T-=maxd
            result.append((maxj,maxd))
        print(result)

    def initial_keys(self):
        for key in self.keys:
            tmp=key.split(',')
            self.keys_int.append([int(t) for t in tmp])


    def b_i_c(self, i, c):
        if i in range(0, self.N):
            if int(c) in self.keys_int[i]:
                
                slen = len(self.keys_int[i])
                # if len(tmp) == 2 and tmp[1] == '':
                #     slen -= 1
                return 1/slen
            else:
                return mt.exp(-20)
        else:
            return mt.exp(-20)  # -20



    def __forward(self,t,i,d):
        if t==0:
            # t=1
            # return self.pi_i[i]*self.b_i_c_arr[i][int(self.cur_seq[0])]*self.DD[(i,d)]
            return self.pi_i[i]
        forward_last = self.alpha[t - 1]
        if d!=self.D:
            tmp1=forward_last[(i,d+1)]*self.b_i_c_arr[i][int(self.cur_seq[t])]
        else:
            tmp1=0
        tmp2=0
        for j in range(self.M):
            if j == i:
                continue
            tmp2 += forward_last[(j, 1)] * self.A[j][i]
        return tmp1 + tmp2 * \
            self.b_i_c_arr[i][int(self.cur_seq[t])] * self.DD[(i, d)]

    def __backward(self,t,i,d):
        if t==len(self.cur_seq)-1:
            return 1
        
        beta_next = self.beta[t + 1]
        if d>1:
            return self.b_i_c_arr[i][int(self.cur_seq[t+1])]*beta_next[(i,d-1)]
        else:
            res = 0
            for j in range(self.M):
                if j==i:
                    continue
                tmp=0
                for d in range(1,self.D+1):
                    tmp+=self.DD[(j,d)]*beta_next[(j,d)]
                res+=self.A[i][j]*self.b_i_c_arr[j][int(self.cur_seq[t+1])]*tmp
            return res

    def __gamma_t(self,t,i):
        res=0
        if t==len(self.cur_seq)-1:
            
            for d in range(1,self.D+1):
                # print('%d,%d'%(i,d))
                res+=self.alpha[t][(i,d)]
            return res
        else:

            for j in range(self.M):
                if j==i:
                    continue
                res+=self.kesi(t+1,i,j)-self.kesi(t+1,j,i)
            return res+self.gamma[t+1][i]