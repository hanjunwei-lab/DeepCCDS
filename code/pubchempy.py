# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 09:22:25 2023

@author: Administrator
"""

import pubchempy
import pandas as pd
import numpy as np


#批量获取pubchem药物信息
with open('D:\\1work\\细胞系数据\\GDSC\\官网\\nonmatch.txt','r',encoding='utf-8-sig') as file1:

    file_lines = file1.readlines()

    name_list=[]

    a=[]

    cc=[]
    
    d1 = []

    d=[]

    e=[]

    f=[]

    g=[]

#readlines读取的每行是字符串格式，采用以下代码将其转换成列表格式

    for i in file_lines:

        j=i.strip() #去掉每行头尾空白

        name_list.append(str(j))

    for k in  name_list:

        results = pubchempy.get_compounds(k, 'name')

        for l in results:

            try:

                #print('CID: {}\tMass: {}\tName: {}\tMolfor: {}\tSmi: {}\tSyn: {}'.format(l.cid,l.exact_mass,l.iupac_name,l.molecular_formula,l.isomeric_smiles,l.synonyms))

                print(l)
                MFs=l.molecular_formula

                MWs=l.molecular_weight
                
                Cas = l.canonical_smiles

                ISs=l.isomeric_smiles

                Sys=l.synonyms

                Cis=l.cid

                a.append(k)

                cc.append(MFs)
                
                d1.append(Cas)

                d.append(ISs)

                e.append(Sys)

                f.append(Cis)

                g.append(MWs)

            except (pubchempy.BadRequestError,TimeoutError,urllib.error.URLError,ValueError):

                pass

            dataframe=pd.DataFrame({'name':a,'molecular_formula':cc,'molecular_weight':g,'Canonical_smiles':d1,'Isomeric_smiles':d,'synonyms':e,'cid':f})

            dataframe.to_csv ("D:\\1work\\深度学习预测药物敏感性\\数据\\化学结构\\GDSC_drug_pubchem1.csv",index=False,sep=',')





###检索特定化合物
from pubchempy import get_compounds, Compound

for compound in get_compounds('Myriocin-12-en', 'name'):

    b1 = compound.cid

    c1 = compound.isomeric_smiles

    d1 = compound.molecular_formula

    e1 = compound.molecular_weight

    f1 = compound.iupac_name

import pandas as pd

dataframe = pd.DataFrame({'molecular_weight': e1,

                          'molecular_formula': d1,

                          'isomeric_smile': c1,

                          'iupac_name': f1,

                          'cid': b1}, index=[0])

dataframe.to_csv("D://1.csv", index=False, sep=',')

#显示所有列

pd.set_option('display.max_columns', None)

#显示所有行

pd.set_option('display.max_rows', None)

#设置value的显示长度为100，默认为50

print(dataframe)


