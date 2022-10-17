import numpy as np



global insert_min, insert_max, insert_skew
insert_min=40
insert_max=850
#insert_min=40
#insert_max=100
insert_kurtosis=3
maxvals=1000000
meanvalue=300
modevalue = np.sqrt(2 / np.pi) * meanvalue


seqlen=20000


def distribution1(minval,maxval):
    distribution=[0 for i in range(0,maxval)]
    return distribution

def distribution1count(minval,maxval):
    distribution_count=[i+minval for i in range(0,maxval)]
    return distribution_count

def distribution2(minval,maxval,factor):
    distribution=[0 for i in range(0,int(maxval/factor))]
    return distribution

def distribution2count(minval,maxval,factor):
    distribution_count=[(i+int(minval/factor))*factor for i in range(0,int((maxval-minval)/factor)+1)]
    return distribution_count


'''
valcounts=distribution2(insert_min,insert_max)
for i in range(0,maxvals):
    insert=np.random.randint(insert_min,insert_max)
    #print(insert)
    valcounts[insert-insert_min]+=1
print(valcounts)
'''

    
    

#x = np.random.rayleigh(scale=2, size=(1, 100))
#print(x)

#y=np.random.rayleigh(scale=2,size=100)
#print(y)
'''
valcounts2=distribution1(int(insert_min/10),int(insert_max/10))
meanvalue=300
modevalue = np.sqrt(2 / np.pi) * meanvalue
print("modevalue %s"%modevalue)
#s = np.random.rayleigh(modevalue, 100)
#print(s)
high=int((insert_max)/10)
low =int(insert_min/10)
print("high %s low %s"%(high,low))
for i in range(0,maxvals):
    num=np.random.rayleigh(modevalue, 1)
    insert=int(num/10)
    #print("num: %s; insert: %s"%(num,insert))
    if low <= insert <= high:
        valcounts2[insert-insert_min]+=1
#print(valcounts2)

valcounts3=distribution1(int(insert_min/10),int(insert_max/10))

high=int((insert_max)/10)
low =int(insert_min/10)
maxpos=len(valcounts3)-1
print("high %s low %s"%(high,low))
for i in range(0,maxvals):
    num=int((meanvalue*np.random.rayleigh(scale=1.0, size=None))/10)
    #insert=int(num/10)
    pos=num-insert_min
    #print("num: %s; pos %s "%(num,pos))
    if pos>=0 and pos <maxpos:
        valcounts3[pos]+=1
        #print("num: %s; pos %s "%(num,pos))
#print(valcounts3)

valcounts4=distribution1(int(insert_min/10),int(insert_max/10))
high=int((insert_max)/10)
low =int(insert_min/10)
maxpos=len(valcounts4)
print("high %s low %s"%(high,low))
for i in range(0,maxvals):
    num=int((meanvalue*np.random.rayleigh(scale=1.0, size=None))/10)+insert_min
    #insert=int(num/10)
    pos=num-insert_min
    #print("num: %s; pos %s "%(num,pos))
    if pos>=0 and pos <maxpos:
        valcounts4[pos]+=1
        #print("num: %s; pos %s "%(num,pos))
#print(valcounts4)
'''

valcounts5=distribution1(int(insert_min/10),int(insert_max/10))
valcounts5_count=distribution2count(insert_min,insert_max,10)
high=int((insert_max)/10)
low =int(insert_min/10)
maxpos=len(valcounts5)
print("high %s low %s"%(high,low))
for i in range(0,maxvals):
    num=int((modevalue*np.random.rayleigh(scale=1.0, size=None))/10)+insert_min
    #insert=int(num/10)
    pos=num-insert_min
    #print("num: %s; pos %s "%(num,pos))
    if pos>=0 and pos <maxpos:
        valcounts5[pos]+=1
        #print("num: %s; pos %s "%(num,pos))
print("valcounts5")
print(valcounts5)
print(valcounts5_count)

def get_insert_length(minlen,maxlen):
    barf=True
    while barf:
        num=int((modevalue*np.random.rayleigh(scale=1.0, size=None))/10)+insert_min
        if minlen <= num <=maxlen:
            barf=False
    #print("num %s; minlen %s; maxlen %s"%(num,minlen,maxlen))
    return num

def get_insert_length2(minlen,maxlen):
    num=int((modevalue*np.random.rayleigh(scale=1.0, size=None))/10)+insert_min
     #print("num %s; minlen %s; maxlen %s"%(num,minlen,maxlen))
    return num


def get_insert_length3(minlen,maxlen):
    num=int(modevalue*np.random.rayleigh(scale=1.0, size=None))+insert_min
     #print("num %s; minlen %s; maxlen %s"%(num,minlen,maxlen))
    return num

def get_insert_length4(minlen,maxlen):
    barf=True
    while barf:
        num=int(modevalue*np.random.rayleigh(scale=1.0, size=None))
        if minlen <= num <=maxlen:
            barf=False
        
     #print("num %s; minlen %s; maxlen %s"%(num,minlen,maxlen))
    return num

'''
valcounts6=distribution1(int(insert_min/10),int(insert_max/10))
high=int((insert_max)/10)
low =int(insert_min/10)
maxpos=len(valcounts6)
print("high %s low %s"%(high,low))
for i in range(0,maxvals):
    num2=get_insert_length(insert_min,insert_max)
    pos=num2-insert_min
    #print("num2: %s; pos %s "%(num2,pos))
    if pos>=0 and pos <maxpos:
        valcounts6[pos]+=1
        #print("num: %s; pos %s "%(num,pos))
#print("valcounts6")
#print(valcounts6)

valcounts6a=distribution1(int(insert_min/10),int(insert_max/10))
high=int((insert_max)/10)
low =int(insert_min/10)
maxpos=len(valcounts6a)
print("high %s low %s"%(high,low))
for i in range(0,maxvals):
    num2=get_insert_length3(insert_min,insert_max)
    pos=int((num2-insert_min)/10)
    #print("num2: %s; pos %s "%(num2,pos))
    if pos>=0 and pos <maxpos:
        valcounts6a[pos]+=1
        #print("num: %s; pos %s "%(num,pos))
#print("valcounts6a")
#print(valcounts6a)

valcounts6b=distribution1(int(insert_min/10),int(insert_max/10))
high=int((insert_max)/10)
low =int(insert_min/10)
maxpos=len(valcounts6b)
print("high %s low %s"%(high,low))
for i in range(0,maxvals):
    num2=get_insert_length4(insert_min,insert_max)
    pos=int((num2)/10)
    if num2< insert_min:
        print("num2: %s; pos %s "%(num2,pos))
    if pos>=0 and pos <maxpos:
        valcounts6b[pos]+=1
        #print("num: %s; pos %s "%(num,pos))
print("valcounts6b")
print(valcounts6b)


valcounts7=distribution2(insert_min,insert_max,10)
valcounts7_count=distribution2count(insert_min,insert_max,10)
high=int((insert_max)/10)
low =int(insert_min/10)
maxpos=len(valcounts7)
print("high %s low %s"%(high,low))
totvals=0
while totvals<= maxvals:
    num=int((modevalue*np.random.rayleigh(scale=1.0, size=None))/10)
    #insert=int(num/10)
    pos=num-int(insert_min)
    #print("num: %s; pos %s "%(num,pos))
    if pos>=0 and pos <maxpos:
        valcounts7[pos]+=1
        #print("num: %s; pos %s "%(num,pos))
        totvals+=1
#print(valcounts7)
#print(valcounts7_count)
'''
