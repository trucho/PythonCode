import os
import numpy as np
from pandas import read_csv
from skimage.measure import label, regionprops, regionprops_table

def loadPoints(viewer,dirPath,csvName,mColor='#555555',mSize=10,mSymbol="disc"):
    csvPath = dirPath + csvName + '.csv'
    if os.path.isfile(csvPath):
        if (read_csv(csvPath).empty==False):
            viewer.open(csvPath, name=csvName, plugin='builtins', symbol=mSymbol, size=mSize, face_color=mColor, edge_color='white')
        else:
            viewer.add_points(name=csvName, symbol=mSymbol, size=mSize, face_color=mColor, edge_color='white')
    else:
        viewer.add_points(name=csvName, symbol=mSymbol, size=mSize, face_color=mColor, edge_color='white')
from skimage.measure import label, regionprops, regionprops_table
def calculateCentroids(viewer,labelsName,outputName,mColor='#555555',mSize=10,mSymbol="disc"):
    props = regionprops_table(viewer.layers[labelsName].data,properties=('label','centroid'))
    centroidData = np.column_stack((props['centroid-0'],props['centroid-1']))
    viewer.add_points(data=centroidData, name = outputName, size = mSize, symbol = mSymbol, face_color= mColor, edge_color='white')
    return centroidData