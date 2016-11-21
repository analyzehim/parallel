import glob
import decimal

def get_z_ef(y):
    import math
    return float((2+2*y)*(math.pow(2,y)*math.pow(2,y))) 


def save_file(mas):
    mas1 = sorted(mas,reverse=True)
    f = open("data.txt","w")
    proc = 0
    f.write('# x1 x2 y\n')
    for value in mas1:
        if value[0]!=proc:
            f.write('\n')
        proc = value[0]
        time = float(value[2])
        op = float(get_z_ef(value[1]))
        ef = op/(time*1000000000)
        st = str(value[0])+" "+str(value[1])+" "+str(ef)+'\n'
        f.write(st)
    f.close()
    
def get(line):
    rank = int(line.split(' ')[0])
    matrsize = int(line.split(' ')[1])
    time = line.split(' ')[2][:-1]
    time = float(decimal.Decimal(time))
    return [rank,matrsize,time]

def get_mas()
    mas=[]
    for filename in glob.glob("data/*.txt"):
        f = open(filename,"r")
        for line in f:
            if line[0]=='_':
                continue
            try:
                mas.append(get(line))
            except:
                print filename, line
    return mas


 
save_file(get_mas())
