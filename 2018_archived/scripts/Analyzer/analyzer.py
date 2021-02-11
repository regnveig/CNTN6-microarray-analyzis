import numpy as np
import os
import sys
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


sys.path.append("/mnt/storage/home/vsfishman/HiC/tutorial_Fishman/utils")
from miscellaneous import startPmDebug
startPmDebug()

design_file = sys.argv[1]
design = open(design_file).readlines()[1:]
design = dict([(line.split("\t")[0],line.split("\t")[1]) for line in design])
data = {}
dataMean = {}

header_line = "FEATURES	FeatureNum	Row	Col	accessions	chr_coord	SubTypeMask	SubTypeName	Start	Sequence	ProbeUID	ControlType	ProbeName	GeneName	SystematicName	Description	PositionX	PositionY	gSurrogateUsed	gIsFound	gProcessedSignal	gProcessedSigError	gNumPixOLHi	gNumPixOLLo	gNumPix	gMeanSignal	gMedianSignal	gPixSDev	gPixNormIQR	gBGNumPix	gBGMeanSignal	gBGMedianSignal	gBGPixSDev	gBGPixNormIQR	gNumSatPix	gIsSaturated	gIsFeatNonUnifOL	gIsBGNonUnifOL	gIsFeatPopnOL	gIsBGPopnOL	IsManualFlag	gBGSubSignal	gBGSubSigError	gIsPosAndSignif	gPValFeatEqBG	gNumBGUsed	gIsWellAboveBG	gBGUsed	gBGSDUsed	ErrorModel	gSpatialDetrendIsInFilteredSet	gSpatialDetrendSurfaceValue	SpotExtentX	SpotExtentY	gNetSignal	gMultDetrendSignal	gProcessedBackground	gProcessedBkngError	IsUsedBGAdjust	gInterpolatedNegCtrlSub	gIsInNegCtrlRange	gIsUsedInMD"
headers = header_line.split("\t")
index = headers.index("gProcessedSignal")
indexSigMean = headers.index("gMeanSignal")
indexSigBG =  headers.index("gBGMeanSignal")

for file in design:
    filepath = os.path.dirname(design_file)
    filepath = filepath + "/" + file
    print "Reading file ", filepath
    data[file] = []
    dataMean[file] = []
    for ind,l in enumerate(open(filepath).readlines()):
        if ind <= 10:
            continue
        l = l.split("\t")
        data[file].append(float(l[index]))
        dataMean[file].append(float(l[indexSigMean])-float(l[indexSigBG]))
    #temp = [len(line.split("\t")) for line in data[file]]
    #assert min(temp)==max(temp)
    #data[file] = np.genfromtxt(data[file],comments="*",delimiter="\t")
    #data[file] = np.genfromtxt(filepath,comments="*",skip_header=10,delimiter="\t")
    #data[file].dtype.names = headers
    #print data[file].dtype
    #raise

pca = PCA(n_components=2,whiten=True)
def do_pca(array,fout):
    pca.fit(array)
    reduced_data = pca.fit_transform(array)
    dist = (max(reduced_data[:, 0]) - min(reduced_data[:, 0])) / 50.
    for i in range(reduced_data.shape[0]):
        plt.text(reduced_data[i][0] + dist, reduced_data[i][1], design[sorted(design.keys())[i]])
        # print(reduced_data[i][0],reduced_data[i][1],design[sorted(design.keys())[i]])
    plt.scatter(reduced_data[:, 0], reduced_data[:, 1], color="blue")
    plt.xlabel("PC1 (" + str(pca.explained_variance_ratio_[0]) + ")")
    plt.ylabel("PC2 (" + str(pca.explained_variance_ratio_[1]) + ")")
    # plt.gca().set_xlim(min(reduced_data[:,0])-1,max(reduced_data[:,0])+1)
    # plt.gca().set_ylim(min(reduced_data[:,1])-1,max(reduced_data[:,1])+1)
    plt.savefig(fout+".png", dpi=300)
    plt.clf()


my_array = np.vstack(tuple([data[k] for k in sorted(data.keys())]))
#do_pca(my_array,"PCAwhiten")

#logArray = np.log(my_array)
#do_pca(my_array,"PCA_Loggedwhiten")

averages = np.average(my_array,axis=0)
assert len(averages)==my_array.shape[1]
t = 90
thr=np.percentile(averages,90)
indexes = averages>thr
high_exp = my_array[:,indexes]
print high_exp.shape
do_pca(high_exp,"PCA_highExpTop"+str(t))

t = 85
thr=np.percentile(averages,t)
indexes = averages<thr
low_exp = my_array[:,indexes]
do_pca(low_exp,"PCA_LowExpBot"+str(t))