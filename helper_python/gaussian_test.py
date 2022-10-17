import numpy as np
global insert_min, insert_max, insert_skew
insert_min=40
insert_max=850
maxvals=100000
meanvalue=200
standard_deviation=10

def distribution1(minval,maxval):
    distribution=[0 for i in range(0,maxval)]
    return distribution

def distribution2count(minval,maxval,factor):
    distribution_count=[(i+int(minval/factor))*factor for i in range(0,int((maxval-minval)/factor)+1)]
    return distribution_count

valcounts10=distribution1(int(insert_min/10),int(insert_max/10))
valcounts10_count=distribution2count(insert_min,insert_max,10)

valcounts1=distribution1(int(insert_min),int(insert_max))
valcounts1_count=distribution2count(insert_min,insert_max,1)

high10=int((insert_max)/10)
low10 =int(insert_min/10)
maxpos10=len(valcounts10)

high1=insert_max
low1 =insert_min
maxpos1=len(valcounts1)


print("high10 %s low10 %s"%(high10,low10))

for i in range(0,maxvals):
    #num=int((modevalue*np.random.rayleigh(scale=1.0, size=None))/10)+insert_min
    #num=int((np.random.normal(loc = meanvalue , scale=standard_deviation, size=None))/10)+insert_min
    num1=int(np.random.normal(loc = meanvalue , scale=standard_deviation, size=None))
    num10=int(num1/10)

    #print("num: %s; pos %s "%(num,pos))
    if num10>=0 and num10 <maxpos10:
        valcounts10[num10]+=1
        #print("num: %s; pos %s "%(num,pos))
    if num1>=0 and num1 <maxpos1:
        valcounts1[num1]+=1
        #print("num: %s; pos %s "%(num,pos))
print("valcounts5")
print(valcounts10)
print(valcounts10_count)
print(valcounts1)
print(valcounts1_count)
