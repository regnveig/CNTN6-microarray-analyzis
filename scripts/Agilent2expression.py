import sys
import pandas
import matplotlib.pyplot as plt
sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/utils")
from miscellaneous import startPmDebug
startPmDebug()

import numpy as np
def readOne(filename):
    print ("Reading",filename)
    lines = open(filename).readlines()
    for l in lines[10:]:
        if len(l.split("\t")) != 62:
            print(l)
            raise
    global header
    header = lines[9].split("\t")
    data = np.genfromtxt(lines,skip_header=10,delimiter="\t",dtype=None,comments=None)
    data.dtype.names = header
    print ("Done")
    CNTN6 = list(data[data["GeneName"]=="CNTN6"][0])
    return CNTN6

design = pandas.read_table(sys.argv[1])
#design = design.head(2)
print (design)
def test(x):
    temp = np.genfromtxt(["1\t1\t1\tabc"],delimiter="\t",dtype=None)
    temp.dtype.names = ["a","b","c","d"]
    temp = list(temp[temp["a"]==1][0])
    return temp

res = pandas.DataFrame(design.apply(lambda x: readOne(x[0]),axis=1).values.tolist())
res.columns = header
print (res)

res = pandas.merge(left=design,right=res,left_index=True, right_index=True)
#res.to_txt("../results/python_CNTN6_expression.txt")
plt.errorbar(np.arange(len(res["Name"]))+1, y = res["gProcessedSignal"], yerr=res["gProcessedSigError"],fmt='o')
plt.gca().set_xlim(0,len(res["Name"])+2)
plt.gca().set_xticks(np.arange(len(res["Name"]))+1)
plt.gca().set_xticklabels(res["Name"],rotation=90)
plt.gca().set_xlabel("Cell line")
plt.gca().set_ylabel("gProcessedSignal")
plt.tight_layout()
plt.savefig("../results/python_CNTN6_expression.png")