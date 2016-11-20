import glob
import decimal

def get(line):

    rank = int(line.split(' ')[0])
    matrsize = int(line.split(' ')[1])
    time = line.split(' ')[2][:-1]
    time = float(decimal.Decimal(time))
    return [rank,matrsize,time]

mas=[]
for filename in glob.glob("*.txt"):
    f = open(filename,"r")
    for line in f:
        if line[0]=='_':
            continue
        mas.append(get(line))

for k in mas:
    if k[1]==17:
        print k
