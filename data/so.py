import glob
import decimal

def save_file(mas):
	f = open("data.txt","w")
	for value in mas:
		st = str(value[0])+" "+str(value[1])+" "+str(value[2])+'\n'
		f.write(st)
	f.close()
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
	try:
       		mas.append(get(line))
	except:
		print filename, line
 
save_file(mas)


