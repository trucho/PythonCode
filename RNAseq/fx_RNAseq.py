import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from scipy.stats import zscore
from cmcrameri import cm #colormaps

def log2(m):
    if m>0:
        return np.log2(m)
    else:
        return 0

def log2series(m):
    tempseries=m
    for i in range(0,len(tempseries)):
        tempseries[i]=log2(tempseries[i])
    return tempseries

def log2matrix(m):
    np.vectorize(log2)

def approx(n):
    return round(n*100)/100


def formatBars(plotH):
    bars = plotH.get_children()
    for i in range(0,len(bars)):
        if i < 7:
            bars[i].set_facecolor('#747474')
        elif i > 6 & i < 12:
            bars[i].set_facecolor('#B540B7')
        elif i > 11 & i < 18:
            bars[i].set_facecolor('#4669F2')
        elif i > 17 & i < 25:
            bars[i].set_facecolor('#04CD22')
        elif i > 24 & i < 31:
            bars[i].set_facecolor('#CC2C2A')
        else:
            bars[i].set_facecolor('#000000')
            

def formatFigure_General(plottitle, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=18)
    fontLabels = font_manager.FontProperties(fname=font_path, size=22)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_title(plottitle, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    
    axH.set_xlabel(axH.get_xlabel(), fontproperties=fontLabels)
    axH.set_ylabel(axH.get_ylabel(), fontproperties=fontLabels)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    return fontLabels
            
def formatFigure_Opsins(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=18)
    fontLabels = font_manager.FontProperties(fname=font_path, size=22)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    # axH.set_xticks(n)
    # axH.set_xticklabels(gdf.iloc[0,h_start:h_end].index);
    axH.set_xticks([3.5,9.5,15.5,22.5,29.5])
    axH.set_xticklabels(['Rods','UV','S','M','L']);

    axH.set_yticks(np.arange(0, 3.5, step=0.5))
    # axH.set_yticklabels([str(i) + ' x 10$^6$' for i in np.arange(0, 2.5, step=0.5)]);
#     axH.set_ylabel('norm. counts x 10$^6$', fontproperties=fontLabels)
    axH.set_ylabel('fpkm', fontproperties=fontLabels)

    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)

#     formatBars(plotH)

def formatFigure_right(genename, figH, axH, plotH):
    formatFigure(genename, figH, axH, plotH)
    axH.set_ylabel('')

def defaultFonts():
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)
    return font_path, fontTicks, fontLabels, fontTitle

def formatFigure(genename, figH, axH, plotH):
    [font_path, fontTicks, fontLabels, fontTitle] = defaultFonts();
    axH.set_xticks([3.5,9.5,15.5,22.5,29.5])
    axH.set_xticklabels(['Rods','UV','S','M','L']);
    axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
#     axH.set_ylabel('norm. counts', fontproperties=fontLabels)
    axH.set_ylabel('fpkm', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)
    
def formatFigureDevCourse(genename, figH, axH, plotH):
    [font_path, fontTicks, fontLabels, fontTitle] = defaultFonts();
    axH.set_xticks([0,1,2,3,4,5])
    axH.set_xticklabels(['PRP','ePR','mPR','lPR','aPR']);
    axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
#     axH.set_ylabel('norm. counts', fontproperties=fontLabels)
    axH.set_ylabel('% expression', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)

def changeTitle(newtitle,axH,tColor):
    [font_path, fontTicks, fontLabels, fontTitle] = defaultFonts();
    axH.set_title(newtitle, fontproperties=fontTitle, color=tColor)


def formatFigure_list(genelist, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    xTicks = []
    xLabels = []
    for i in range(len(genelist)):
        xTicks = np.append(xTicks, [3.5+(10*i),9.5+(10*i),15.5+(10*i),22.5+(10*i),29.5+(10*i)])
        xLabels = np.append(xLabels, ['Rods','UV','S','M','L'])
    axH.set_xticks(xTicks)
    axH.set_xticklabels(xLabels);
    axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
#     axH.set_ylabel('norm. counts', fontproperties=fontLabels)
    axH.set_ylabel('fpkm', fontproperties=fontLabels)
    axH.set_title(genelist, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 18)
    axH.yaxis.offsetText.set_fontsize(18)

def formatFigure_nReads(plotname, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=18)
    fontLabels = font_manager.FontProperties(fname=font_path, size=22)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([3.5,9.5,15.5,22.5,29.5])
    axH.set_xticklabels(['Rods','UV','S','M','L']);
    axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
    axH.set_ylabel('mapped reads x $10^6$', fontproperties=fontLabels)
    axH.set_title(plotname, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 18)
    axH.yaxis.offsetText.set_fontsize(18)


def heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=36)
    fontLabels = font_manager.FontProperties(fname=font_path, size=22)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    if data.shape[0]==0:
        new_data = pd.DataFrame(np.zeros([1,data.shape[1]]), columns=data.columns)
        data = data.append(new_data)

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    # perceptually responsible colormaps are: inferno, viridis, plasma, magma, cividis
    im = ax.imshow(data, cmap = "bone", **kwargs)
#     im = ax.imshow(data, cmap = "inferno", **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', shrink=.9, pad=0.05, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, fontproperties=fontLabels, rotation=0,ha="right", va="center",rotation_mode="anchor")
    cbar.ax.tick_params(labelsize=22)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontproperties=fontLabels)
    ax.set_yticklabels(row_labels, fontproperties=fontLabels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_yticklabels(), rotation=30, ha="right", va="center",rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
#         spine.set_visible(False)
        spine.set_linewidth(.5)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

#     ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
#     ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)
#     #all-white grid
#     ax.grid(which="minor", color="w", linestyle='-', linewidth=1)


    # Custom grid according to photoreceptor subtype
    for h in np.arange(-.5,data.shape[0]+.5):
        ax.axhline(y = h, color = 'black', linewidth = 2, alpha = 1, solid_capstyle='butt')
    for v in np.arange(-.5,data.shape[1]+.5):
        ax.axvline(x = v, color = 'black', linewidth = 2, alpha = 1, solid_capstyle='butt')

    ax.axvline(x = -0.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = np.sum(groupsN)-.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')

    for i in np.arange(groupsN.shape[0]):
        ax.axvline(x = np.sum(groupsN[:i+1])-.5, color = 'white', linewidth = 3, alpha = 1, solid_capstyle='butt')
        ax.plot([np.sum(groupsN[:i])-.5,np.sum(groupsN[:i+1])-.5], [-.5,-.5], '-', lw=8, color = groupsColors[i], solid_capstyle='butt')
        ax.plot([np.sum(groupsN[:i])-.5,np.sum(groupsN[:i+1])-.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=8, color = groupsColors[i], solid_capstyle='butt')
        ax.text(((np.sum(groupsN[:i])+np.sum(groupsN[:i+1]))/2)-.5, -1.2, groupsLabels[i], color = groupsColors[i], horizontalalignment='center', fontproperties=fontTicks)

    return im, cbar

def heatmap(data, row_labels, col_labels, ax=None, cbarlabel="fpkm",
            cbar_kw={}, **kwargs):

    groupsN = np.array([6,5,6,7,6])
    groupsColors = np.array(['#747474','#B540B7','#4669F2','#04CD22','#CC2C2A'])
    groupsLabels = np.array(['Rods','UV','S','M','L'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel=cbarlabel)

    return im, cbar

def heatmap_general_z(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=18)
    fontLabels = font_manager.FontProperties(fname=font_path, size=12)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    if data.shape[0]==0:
        new_data = pd.DataFrame(np.zeros([1,data.shape[1]]), columns=data.columns)
        data = data.append(new_data)

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    # perceptually responsible colormaps are: inferno, viridis, plasma, magma, cividis
#     im = ax.imshow(data, cmap = "bone", **kwargs)
#     im = ax.imshow(data, cmap = "inferno", **kwargs)
    divnorm = matplotlib.colors.TwoSlopeNorm(vmin=data.min(), vcenter=0, vmax=data.max())
    im = ax.imshow(data, cmap = cm.berlin, norm=divnorm, **kwargs)
    
#     im = ax.imshow(data, cmap = cm.berlin, vmin=-5, vmax=5, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', shrink=.7, pad=0.02, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontproperties=fontLabels)
    ax.set_yticklabels(row_labels, fontproperties=fontLabels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
#         spine.set_visible(False)
        spine.set_linewidth(.5)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

#     ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
#     ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.tick_params(which="minor", bottom=False, left=False)
#     #all-white grid
#     ax.grid(which="minor", color="w", linestyle='-', linewidth=1)


    # Custom grid according to photoreceptor subtype
    for h in np.arange(-.5,data.shape[0]+.5):
        ax.axhline(y = h, color = 'black', linewidth = 2, alpha = 1, solid_capstyle='butt')
    for v in np.arange(-.5,data.shape[1]+.5):
        ax.axvline(x = v, color = 'black', linewidth = 2, alpha = 1, solid_capstyle='butt')

    ax.axvline(x = -.5, color = 'white', linewidth = 1, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = np.sum(groupsN)-.5, color = 'white', linewidth = 1, alpha = 1, solid_capstyle='butt')

    for i in np.arange(groupsN.shape[0]):
        ax.axvline(x = np.sum(groupsN[:i+1])-.5, color = 'white', linewidth = 1, alpha = 1, solid_capstyle='butt')
        ax.plot([np.sum(groupsN[:i])-.5,np.sum(groupsN[:i+1])-.5], [-.5,-.5], '-', lw=8, color = groupsColors[i], solid_capstyle='butt')
        ax.plot([np.sum(groupsN[:i])-.5,np.sum(groupsN[:i+1])-.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=8, color = groupsColors[i], solid_capstyle='butt')
        ax.text(((np.sum(groupsN[:i])+np.sum(groupsN[:i+1]))/2)-.5, -.8, groupsLabels[i], color = groupsColors[i], horizontalalignment='center', fontproperties=fontTicks)

    return im, cbar

def heatmap_z(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([6,5,6,7,6])
    groupsColors = np.array(['#747474','#B540B7','#4669F2','#04CD22','#CC2C2A'])
    groupsLabels = np.array(['Rods','UV','S','M','L'])
    
    data = zscore(data,axis=1, nan_policy='omit') # calculate z-score

    im,cbar = heatmap_general_z(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="z-score")

    return im, cbar

def heatmap_glia_showAllSamples(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([3,3,3,3,3,3, 3,3,3,3,3,3])
    groupsColors = np.array(['#c7ffec','#94ffdb','#57ffc7','#2effb9','#05ffac','#00d68f', '#e3e3e3','#c9c9c9','#ababab','#969696','#828282','#6b6b6b'])
    groupsLabels = np.array(['2','2.5','3','4','5','8', '2','2.5','3','4','5','8'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="ttm")

    return im, cbar

def heatmap_glia(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):
    if data.shape[0]==0:
        new_data = pd.DataFrame(np.zeros([1,data.shape[1]]), columns=data.columns)
        data = data.append(new_data)

    # Calculate averages
    meandata = pd.DataFrame(columns=['MG48','MG60','MG72','MG96','MG120','MG192','C48','C60','C72','C96','C120','C192'])
    meandata['MG48']=data.iloc[:,0:3].mean(axis='columns')
    meandata['MG60']=data.iloc[:,3:6].mean(axis='columns')
    meandata['MG72']=data.iloc[:,6:9].mean(axis='columns')
    meandata['MG96']=data.iloc[:,9:12].mean(axis='columns')
    meandata['MG120']=data.iloc[:,12:15].mean(axis='columns')
    meandata['MG192']=data.iloc[:,15:18].mean(axis='columns')
    meandata['C48']=data.iloc[:,18:21].mean(axis='columns')
    meandata['C60']=data.iloc[:,21:24].mean(axis='columns')
    meandata['C72']=data.iloc[:,24:27].mean(axis='columns')
    meandata['C96']=data.iloc[:,27:30].mean(axis='columns')
    meandata['C120']=data.iloc[:,30:33].mean(axis='columns')
    meandata['C192']=data.iloc[:,33:36].mean(axis='columns')


    groupsN = np.array([1,1,1,1,1,1, 1,1,1,1,1,1])
    groupsColors = np.array(['#c7ffec','#94ffdb','#57ffc7','#2effb9','#05ffac','#00d68f', '#e3e3e3','#c9c9c9','#ababab','#969696','#828282','#6b6b6b'])
    groupsLabels = np.array(['2','2.5','3','4','5','8', '2','2.5','3','4','5','8'])

    im,cbar = heatmap_general(meandata, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="ttm")

    return im, cbar

def formatFigure_glia(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=16)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([2,5,8,11,14,17,21.5,24.5,27.5,30.5,33.5,36.5])
    axH.set_xticklabels(['2','2.5','3','4','5','8','2','2.5','3','4','5','8']);
    axH.text(0.87, 0.060, "dpf", transform=figH.transFigure, fontproperties=fontTicks);
    axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
#     axH.set_ylabel('norm. counts', fontproperties=fontLabels)
    axH.set_ylabel('fpkm', fontproperties=fontLabels)
    axH.set_xlabel('Muller Glia', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 16)
    axH.yaxis.offsetText.set_fontsize(24)




def formatFigure_rods(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([2.5,7])
    axH.set_xticklabels(['Rods','Misc.']);
    axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
#     axH.set_ylabel('norm. counts', fontproperties=fontLabels)
    axH.set_ylabel('notsure', fontproperties=fontLabels)
#     axH.set_xlabel('Muller Glia', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)


def heatmap_rods(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([4,4])
    groupsColors = np.array(['#747474','#dac910'])
    groupsLabels = np.array(['rods','~rods'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="cpm")

    return im, cbar

def formatFigure_sqcones(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([3,6.5,10,13])
    axH.set_xticklabels(['S','S_hib','M','M_hib']);
    axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
#     axH.set_ylabel('norm. counts', fontproperties=fontLabels)
    axH.set_ylabel('counts', fontproperties=fontLabels)
#     axH.set_xlabel('Squirrel cones', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)

def heatmap_sqcones(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([5,2,3,3])
    groupsColors = np.array(['#4669F2','#548ced','#04CD22','#53e477'])
    groupsLabels = np.array(['S_awk','S_hib','M_awk','M_hib'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="counts")

    return im, cbar

def formatFigure_haircell(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([2.5,7])
    axH.set_xticklabels(['hairCell','Misc.']);
    axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
    axH.set_ylabel('rpkm', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)

def heatmap_haircell(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([3,3])
    groupsColors = np.array(['#e147c0','#dac910'])
    groupsLabels = np.array(['hairCells','~hairCells'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="rpkm")

    return im, cbar


def formatFigure_zfHoang2020(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([1,2,3,4,5,6,7,8,9,10,11])
    axH.set_xticklabels(['PRP','eslPR','mslPR','lslPR','adPR','lslR','R','UV','S','M','L']);
#     axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
    axH.set_ylabel('Mean Exp', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    axH.xaxis.set_tick_params(rotation=66)
    
#     axH.set_ylim([0,100])

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)


def heatmap_zfHoang2020(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([1,1,1,1,1,1,1,1,1,1,1])
    groupsColors = np.array(['#dfdac8','#dacd9a','#dcc360','#cca819','#ffd429','#a3a3a3','#7d7d7d','#B540B7','#4669F2','#04CD22','#CC2C2A'])
    groupsLabels = np.array(['PRP','eslPR','mslPR','lslPR','adPR','lslR','R','UV','S','M','L'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="avg.")

    return im, cbar


def formatFigure_zflarvaHoang2020(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([1,2,3,4,5,6])
    axH.set_xticklabels(['PRP','eslPR','mslPR','lslPR','adPR','lslR',]);
#     axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
    axH.set_ylabel('% Exp', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    
    axH.set_ylim([0,100])

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)

def heatmap_zflarvaHoang2020(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([1,1,1,1,1,1])
    groupsColors = np.array(["#dfdac8",'#dacd9a','#dcc360','#cca819','#ffd429','#474747'])
    groupsLabels = np.array(['PRP','eslPR','mslPR','lslPR','adPR','lslR',])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="avg.")

    return im, cbar

def formatFigure_zfadultHoang2020(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([1,2,3,4,5])
    axH.set_xticklabels(['R','UV','S','M','L']);
#     axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
    axH.set_ylabel('% Exp', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    
#     axH.set_ylim([0,100])

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)
    

def heatmap_zfadultHoang2020(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([1,1,1,1,1])
    groupsColors = np.array(['#747474','#B540B7','#4669F2','#04CD22','#CC2C2A'])
    groupsLabels = np.array(['R','UV','S','M','L'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="avg.")

    return im, cbar

def formatFigure_zfXu2020(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([1,2,3,4,5,6])
    axH.set_xticklabels(['mslR','mslUV','mslS','mslM','mslL','mslPR']);
#     axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
    axH.set_ylabel('Mean Exp', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    axH.xaxis.set_tick_params(rotation=66)
    
#     axH.set_ylim([0,100])

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)


def heatmap_zfXu2020(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([1,1,1,1,1,1])
    groupsColors = np.array(['#7d7d7d','#B540B7','#4669F2','#04CD22','#CC2C2A','#dcc360'])
    groupsLabels = np.array(['mslR','mslUV','mslS','mslM','mslL','mslPR'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="avg.")

    return im, cbar

def formatFigure_zfOgawa2021(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)

    axH.set_xticks([1,2,3,4,5,6,7,8])
    axH.set_xticklabels(['R','UV','S','M','L','M$_{4}$','B$_{on}$','B$_{off}$']);
#     axH.ticklabel_format(style='sci',axis='y',scilimits=(0,2))
#     axH.set_ylabel('% Exp', fontproperties=fontLabels)
    axH.set_ylabel('Avg. Exp', fontproperties=fontLabels)
    axH.set_title(genename, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    
#     axH.set_ylim([0,100])

    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.tick_params(axis = 'both', which = 'major', labelsize = 24)
    axH.yaxis.offsetText.set_fontsize(24)
    

def heatmap_zfOgawa2021(data, row_labels, col_labels, ax=None,
            cbar_kw={}, **kwargs):

    groupsN = np.array([1,1,1,1,1,1,1,1])
    groupsColors = np.array(['#7d7d7d','#B540B7','#4669F2','#04CD22','#CC2C2A','#cdcd04','#ccf2ff','#663d00'])
    groupsLabels = np.array(['R','UV','S','M','L','M$_{4}$','B$_{on}$','B$_{off}$'])

    im,cbar = heatmap_general(data, row_labels, col_labels, groupsN, groupsColors, groupsLabels, ax=ax, cbarlabel="avg.")

    return im, cbar