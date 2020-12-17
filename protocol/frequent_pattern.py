import matplotlib.pyplot as plt

class frequent_pattern:
    def __init__(self,seqs,threshold):
        self.seqs=seqs
        self.threshold=threshold    #最小支持度阈值

    def freq(self,query):
        '''
        :param query:tuple,待查询的串
        :return: 计算query在每个观察序列中出现的次数，超过阈值则返回True
        '''
        cont=0
        query_str=str(query)
        query_str=query_str[1:len(query_str)-1]

        for seq in self.seqs:
            seq=str(seq)
            seq=seq[1:len(seq)-1]
            if query_str in seq:
                cont += 1

        if cont>=self.threshold:
            return True,cont
        else:
            return False,0

    # def byteinseq(self,target,seq):
    #     slen=len(seq)
    #     for i in range(0,slen,2):
    #         if target==seq[i:i+len(target)]:
    #             return True
    #     return False

    def apriori(self):
        cachelist=[]#已经计算过的模式串
        collects=[{}]
        collect1=collects[0]
        for i in range(256):
            # hexstr=str(hex(i)).split('0x')[1]
            # if len(hexstr)!=2:
            #     hexstr='0'+hexstr
            cont=0
            for seq in self.seqs:
                if i in seq:
                    cont+=1
            if cont>=self.threshold:
                tmp=(i,)
                collect1[tmp]=cont


        i=0
        while i<3 and len(collects[i])!=0:
            keys=list(collects[i].keys())   #频繁i项集[(),(),...]
            keys_len=len(keys)
            # if i!=0:
            #     tmpdict={}
            #     for j in range(keys_len):
            #
            #         flag,sup=self.freq(keys[j])
            #         if flag:
            #             tmpdict[keys[j]]=sup
            #     collects[i]=tmpdict

            newcollect={}
            for j in range(keys_len):
                for k in range(j+1,keys_len):

                    tmp1=keys[j]
                    tmp2=keys[k]
                    # print(i>=1 and tmp1[:i*2]==tmp2[2:])
                    if i>=1:
                        if tmp1[:i]==tmp2[1:]:
                            newstr=tmp2+(tmp1[i],)
                        elif tmp2[:i]==tmp1[1:]:
                            newstr = tmp1 + (tmp2[i],)
                    elif i==0:
                        newstr=tmp1+tmp2
                        newstr2=tmp2+tmp1
                        if newstr2 not in cachelist:
                            # print(newstr2)
                            flag, sup = self.freq(newstr2)
                            cachelist.append(newstr2)
                            if flag:
                                newcollect[newstr2] = sup
                    if newstr not in cachelist:
                        # print(newstr)
                        flag, sup = self.freq(newstr)
                        cachelist.append(newstr)
                        if flag:
                            newcollect[newstr] = sup
            collects.append(newcollect)
            i+=1


        # x=collect1.keys()
        # y=collect1.values()
        # print(y)
        # plt.plot(x,y)
        # plt.show()

        print(collects)

        i=0
        while i<len(collects)-1 and len(collects[i+1])!=0:
            keys1=list(collects[i].keys())
            keys2=list(collects[i+1].keys())
            for key1 in keys1:
                strkey1=str(key1)
                if i==0:
                    strkey1=strkey1[1:len(strkey1)-1]
                else:
                    strkey1=strkey1[1:len(strkey1)-2]   #去掉(0,)中的,
                for key2 in keys2:
                    strkey2=str(key2)
                    strkey2=strkey2[1:len(strkey2)-1]
                    if strkey1 in strkey2:
                        # print(key1)
                        del collects[i][key1]   #如果i+1项集的key包含i项集的key，则将i项集的key删掉
                        break
            i+=1
        print(collects)

        result_file=open('frequent_result.txt','w')
        for collect in collects:
            if len(collect)==0:
                continue
            for keys,value in collect.items():
                toKey=''
                for key in keys:
                    toKey+=str(key)+','
                toKey=toKey[0:len(toKey)-1]
                result_file.write(str(toKey)+' '+str(value)+'\n')
            # result_file.write(str(collect.items())+'\n')
        result_file.close()