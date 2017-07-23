import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

projectDir = '../work_solitary'
depth = (np.loadtxt('../work_solitary/depth_levee.txt'))*-1  # Change Depth to (-) under MWL and (+) over MWL.
[pts,numOfRows]  = depth.shape                   # num of Colums = number of points 

while True:
    try:
        numOfETA = int(input('Enter the Number of ETA Text Files: '))
        Lt = float(input('Enter the Total Horizontal Lenght (m): '))
        WaveMaker = float(input('Enter the Wavemaker x location (m): '))
        plotInt = float(input('Enter the plot interval: '))
        outputDataDir = projectDir
        break
    except ValueError:
        print("Oops! That was no valid number. Try again... \n")
        
# change to your corresponding directory which contains eta_00... files (text)
print("Output Directory is: ", outputDataDir)
postprocessDir = outputDataDir+'/output'

# data placeholder
eta = np.zeros((numOfRows,pts,numOfETA))
freeSurface = np.zeros((numOfRows,numOfETA))

def readETAData(num,numOfETA,postprocessDir):
    """Function that reads in the eta_00.. ASCII file and returns
    it in the 2D NumPy array of size [numOfRows,pts]."""

    assert num in range(1,numOfETA+1), "File Index:%d is not in range of station numbers." %(num,)
    if num < 10:
        fileIndex = str(num)
        fileName = postprocessDir+'/eta_000{0:s}'.format(fileIndex)

    elif num <100 and num > 9:
        fileIndex = str(num)
        fileName = postprocessDir+'/eta_00{0:s}'.format(fileIndex)

    elif num < 1000 and num>99: 
        fileIndex = str(num)
        fileName = postprocessDir+'/eta_0{0:s}'.format(fileIndex)

    else: 
        fileIndex = str(i)
        fileName = postprocessDir+'/eta_{0:s}'.format(fileIndex)

    fin = open(fileName, "r")
    print("READING IN: eta # %s" %(fileIndex))
    lines = fin.readlines()
    numpyFileArray = np.zeros((len(lines),pts))

    for lineIndex,line in enumerate(lines):
        rowDataList = line.strip().split()          # dropping all whitespace (leading & trailing)
        numpyData = np.array(rowDataList).astype('float')
        numpyFileArray[lineIndex,:] = numpyData[:]

    fin.close()
    return numpyFileArray

# Plot and save every ETA file:
for i in range(1, numOfETA+1):
    surf = readETAData(i,numOfETA,postprocessDir)   # Y axis = eta data
    x = np.linspace(0, Lt, pts)                     # X axis = Total Length (m)
    
    fig  = plt.figure(figsize=(18,4), dpi=600)
    ax = fig.add_subplot(1,1,1)
    plt.plot(x, surf[1,:], 'c-', linewidth = 0.2)
    plt.axis([-1.5,Lt,min(depth[1,:])-.05,5])

    plt.xlabel('Length (%4.2f m)' % (Lt), fontsize = 12, fontweight = 'bold')
    plt.ylabel('Height (m)', fontsize = 12, fontweight = 'bold')
    

    # Water Fill:
    plt.fill_between(x, depth[1,:], surf[1,:],
                     where = surf[1,:] > depth[1,:],
                     facecolor = 'cyan', interpolate =True)
    
    # Bottom Fill:
    plt.fill_between(x, min(depth[1,:])-.05, depth[1,:], 
                               where= depth[1,:] > (depth[1,:]-.05),facecolor = '0.35',
                               hatch = 'X')

    #Annotations:
    an = plt.annotate('\n~\n|', xy=(WaveMaker,.05), xycoords='data',fontsize = 20,
                                 ha='center', va='center')
    Time = plt.annotate('Time: %4.2f sec'%(i*plotInt-plotInt), xy=(Lt*.85,3),fontsize = 16,
                                 ha='center', va='center')

    if i < 10:
        fileIndex = '0'+str(i)
        fileName = postprocessDir+'/surf_00{0:s}'.format(fileIndex)
        plt.savefig(fileName, ext="png", bbox_inches='tight')
        plt.close()
    elif i < 100:
        fileIndex = str(i)
        fileName = postprocessDir+'/surf_00{0:s}'.format(fileIndex)
        plt.savefig(fileName, ext="png", bbox_inches='tight')
        plt.close()
    elif i < 1000: 
        fileIndex = str(i)
        fileName = postprocessDir+'/surf_0{0:s}'.format(fileIndex)
        plt.savefig(fileName, ext="png", bbox_inches='tight')
        plt.close()
    else: 
        fileIndex = str(i)
        fileName = postprocessDir+'/surf_{0:s}'.format(fileIndex)
        plt.savefig(fileName,ext="png", bbox_inches='tight')
        plt.close()
