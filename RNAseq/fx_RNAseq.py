import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

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
    

def formatFigure(genename, figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=24)
    fontLabels = font_manager.FontProperties(fname=font_path, size=28)
    fontTitle = font_manager.FontProperties(fname=font_path, size=28)
    
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

    
def heatmap(data, row_labels, col_labels, ax=None,
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
    
    #photoreceptor Colors
    pC = {
        'r' : '#747474',
        'u' : '#B540B7',
        's' : '#4669F2',
        'm' : '#04CD22',
        'l' : '#CC2C2A',
        'plt' : '',
    }

    pC['plt']=[
        pC['r'],pC['r'],pC['r'],pC['r'],pC['r'],pC['r'],
        pC['u'],pC['u'],pC['u'],pC['u'],pC['u'],
        pC['s'],pC['s'],pC['s'],pC['s'],pC['s'],pC['s'],
        pC['m'],pC['m'],pC['m'],pC['m'],pC['m'],pC['m'],pC['m'],
        pC['l'],pC['l'],pC['l'],pC['l'],pC['l'],pC['l']
    ]

    if data.shape[0]==0:
        new_data = pd.DataFrame(np.zeros([1,data.shape[1]]), columns=data.columns)
        data = data.append(new_data)
        
    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

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
    ax.axvline(x = 5.5, color = 'white', linewidth = 1, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 10.5, color = 'white', linewidth = 1, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 16.5, color = 'white', linewidth = 1, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 23.5, color = 'white', linewidth = 1, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 29.45, color = 'white', linewidth = 1, alpha = 1, solid_capstyle='butt')
    
    plt.plot([-.5,5.5], [-.5,-.5], '-', lw=5, color = pC['r'], solid_capstyle='butt')
    plt.plot([5.5,10.5], [-.5,-.5], '-', lw=5, color = pC['u'], solid_capstyle='butt')
    plt.plot([10.5,16.5], [-.5,-.5], '-', lw=5, color = pC['s'], solid_capstyle='butt')
    plt.plot([16.5,23.5], [-.5,-.5], '-', lw=5, color = pC['m'], solid_capstyle='butt')
    plt.plot([23.5,29.5], [-.5,-.5], '-', lw=5, color = pC['l'], solid_capstyle='butt')
    
    plt.plot([-.5,5.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = pC['r'], solid_capstyle='butt')
    plt.plot([5.5,10.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = pC['u'], solid_capstyle='butt')
    plt.plot([10.5,16.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = pC['s'], solid_capstyle='butt')
    plt.plot([16.5,23.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = pC['m'], solid_capstyle='butt')
    plt.plot([23.5,29.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = pC['l'], solid_capstyle='butt')
    
    plt.text(2.5, -.8, 'Rods', color = pC['r'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(8, -.8, 'UV', color = pC['u'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(13.5, -.8, 'S', color = pC['s'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(20, -.8, 'M', color = pC['m'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(26.5, -.8, 'L', color = pC['l'], horizontalalignment='center', fontproperties=fontTicks)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[im.norm(data[i, j]) > threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def heatmap_glia_showAllSamples(data, row_labels, col_labels, ax=None,
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
    
    #glia Colors
    gC = {
        'm48' : '#c7ffec',
        'm60' : '#94ffdb',
        'm72' : '#57ffc7',
        'm96' : '#2effb9',
        'm120': '#05ffac',
        'm192': '#00d68f',

        'c48' : '#e3e3e3',
        'c60' : '#c9c9c9',
        'c72' : '#ababab',
        'c96' : '#969696',
        'c120': '#828282',
        'c192': '#6b6b6b',

        'plt' : '',
    }

    gC['plt']=[
        gC['m48'],gC['m48'],gC['m48'],
        gC['m60'],gC['m60'],gC['m60'],
        gC['m72'],gC['m72'],gC['m72'],
        gC['m96'],gC['m96'],gC['m96'],
        gC['m120'],gC['m120'],gC['m120'],
        gC['m192'],gC['m192'],gC['m192'],
        gC['c48'],gC['c48'],gC['c48'],
        gC['c60'],gC['c60'],gC['c60'],
        gC['c72'],gC['c72'],gC['c72'],
        gC['c96'],gC['c96'],gC['c96'],
        gC['c120'],gC['c120'],gC['c120'],
        gC['c192'],gC['c192'],gC['c192'],
    ]

    if data.shape[0]==0:
        new_data = pd.DataFrame(np.zeros([1,data.shape[1]]), columns=data.columns)
        data = data.append(new_data)
        
    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

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
        
    ax.axvline(x = -.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 17.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 35.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')

    
    plt.plot([-.5,2.5], [-.5,-.5], '-', lw=5, color = gC['m48'], solid_capstyle='butt')
    plt.plot([2.5,5.5], [-.5,-.5], '-', lw=5, color = gC['m60'], solid_capstyle='butt')
    plt.plot([5.5,8.5], [-.5,-.5], '-', lw=5, color = gC['m72'], solid_capstyle='butt')
    plt.plot([8.5,11.5], [-.5,-.5], '-', lw=5, color = gC['m96'], solid_capstyle='butt')
    plt.plot([11.5,14.5], [-.5,-.5], '-', lw=5, color = gC['m120'], solid_capstyle='butt')
    plt.plot([14.5,17.5], [-.5,-.5], '-', lw=5, color = gC['m192'], solid_capstyle='butt')
    
    plt.plot([-.5+18,2.5+18], [-.5,-.5], '-', lw=5, color = gC['c48'], solid_capstyle='butt')
    plt.plot([2.5+18,5.5+18], [-.5,-.5], '-', lw=5, color = gC['c60'], solid_capstyle='butt')
    plt.plot([5.5+18,8.5+18], [-.5,-.5], '-', lw=5, color = gC['c72'], solid_capstyle='butt')
    plt.plot([8.5+18,11.5+18], [-.5,-.5], '-', lw=5, color = gC['c96'], solid_capstyle='butt')
    plt.plot([11.5+18,14.5+18], [-.5,-.5], '-', lw=5, color = gC['c120'], solid_capstyle='butt')
    plt.plot([14.5+18,17.5+18], [-.5,-.5], '-', lw=5, color = gC['c192'], solid_capstyle='butt')
    
    plt.plot([-.5,2.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m48'], solid_capstyle='butt')
    plt.plot([2.5,5.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m60'], solid_capstyle='butt')
    plt.plot([5.5,8.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m72'], solid_capstyle='butt')
    plt.plot([8.5,11.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m96'], solid_capstyle='butt')
    plt.plot([11.5,14.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m120'], solid_capstyle='butt')
    plt.plot([14.5,17.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m192'], solid_capstyle='butt')

    plt.plot([-.5+18,2.5+18], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c48'], solid_capstyle='butt')
    plt.plot([2.5+18,5.5+18], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c60'], solid_capstyle='butt')
    plt.plot([5.5+18,8.5+18], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c72'], solid_capstyle='butt')
    plt.plot([8.5+18,11.5+18], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c96'], solid_capstyle='butt')
    plt.plot([11.5+18,14.5+18], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c120'], solid_capstyle='butt')
    plt.plot([14.5+18,17.5+18], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c192'], solid_capstyle='butt')
    
    plt.text(1, -.8, '2', color = gC['m48'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(4, -.8, '2.5', color = gC['m60'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(7, -.8, '3', color = gC['m72'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(10, -.8, '4', color = gC['m96'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(13, -.8, '5', color = gC['m120'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(16, -.8, '8', color = gC['m192'], horizontalalignment='center', fontproperties=fontTicks)

    plt.text(1+18, -.8, '2', color = gC['c48'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(4+18, -.8, '2.5', color = gC['c60'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(7+18, -.8, '3', color = gC['c72'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(10+18, -.8, '4', color = gC['c96'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(13+18, -.8, '5', color = gC['c120'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(16+18, -.8, '8', color = gC['c192'], horizontalalignment='center', fontproperties=fontTicks)
    

    return im, cbar

def heatmap_glia(data, row_labels, col_labels, ax=None,
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
    
    #glia Colors
    gC = {
        'm48' : '#c7ffec',
        'm60' : '#94ffdb',
        'm72' : '#57ffc7',
        'm96' : '#2effb9',
        'm120': '#05ffac',
        'm192': '#00d68f',

        'c48' : '#e3e3e3',
        'c60' : '#c9c9c9',
        'c72' : '#ababab',
        'c96' : '#969696',
        'c120': '#828282',
        'c192': '#6b6b6b',

        'plt' : '',
    }

    gC['plt']=[
        gC['m48'],gC['m48'],gC['m48'],
        gC['m60'],gC['m60'],gC['m60'],
        gC['m72'],gC['m72'],gC['m72'],
        gC['m96'],gC['m96'],gC['m96'],
        gC['m120'],gC['m120'],gC['m120'],
        gC['m192'],gC['m192'],gC['m192'],
        gC['c48'],gC['c48'],gC['c48'],
        gC['c60'],gC['c60'],gC['c60'],
        gC['c72'],gC['c72'],gC['c72'],
        gC['c96'],gC['c96'],gC['c96'],
        gC['c120'],gC['c120'],gC['c120'],
        gC['c192'],gC['c192'],gC['c192'],
    ]

    if data.shape[0]==0:
        new_data = pd.DataFrame(np.zeros([1,data.shape[1]]), columns=data.columns)
        data = data.append(new_data)
        
    if not ax:
        ax = plt.gca()

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
    # Plot the heatmap
    im = ax.imshow(meandata, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', shrink=.7, pad=0.02, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel)

    # We want to show all ticks...
    ax.set_xticks(np.arange(meandata.shape[1]))
    ax.set_yticks(np.arange(meandata.shape[0]))
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
    for h in np.arange(-.5,meandata.shape[0]+.5):
        ax.axhline(y = h, color = 'black', linewidth = 2, alpha = 1, solid_capstyle='butt')
    for v in np.arange(-.5,meandata.shape[1]+.5):
        ax.axvline(x = v, color = 'black', linewidth = 2, alpha = 1, solid_capstyle='butt')
        
    ax.axvline(x = -.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 5.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 11.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')

   
    plt.plot([-.5,.5], [-.5,-.5], '-', lw=5, color = gC['m48'], solid_capstyle='butt')
    plt.plot([.5,1.5], [-.5,-.5], '-', lw=5, color = gC['m60'], solid_capstyle='butt')
    plt.plot([1.5,2.5], [-.5,-.5], '-', lw=5, color = gC['m72'], solid_capstyle='butt')
    plt.plot([2.5,3.5], [-.5,-.5], '-', lw=5, color = gC['m96'], solid_capstyle='butt')
    plt.plot([3.5,4.5], [-.5,-.5], '-', lw=5, color = gC['m120'], solid_capstyle='butt')
    plt.plot([4.5,5.5], [-.5,-.5], '-', lw=5, color = gC['m192'], solid_capstyle='butt')

    plt.plot([-.5+6,.5+6], [-.5,-.5], '-', lw=5, color = gC['c48'], solid_capstyle='butt')
    plt.plot([.5+6,1.5+6], [-.5,-.5], '-', lw=5, color = gC['c60'], solid_capstyle='butt')
    plt.plot([1.5+6,2.5+6], [-.5,-.5], '-', lw=5, color = gC['c72'], solid_capstyle='butt')
    plt.plot([2.5+6,3.5+6], [-.5,-.5], '-', lw=5, color = gC['c96'], solid_capstyle='butt')
    plt.plot([3.5+6,4.5+6], [-.5,-.5], '-', lw=5, color = gC['c120'], solid_capstyle='butt')
    plt.plot([4.5+6,5.5+6], [-.5,-.5], '-', lw=5, color = gC['c192'], solid_capstyle='butt')

    plt.plot([-.5,.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m48'], solid_capstyle='butt')
    plt.plot([.5,1.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m60'], solid_capstyle='butt')
    plt.plot([1.5,2.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m72'], solid_capstyle='butt')
    plt.plot([2.5,3.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m96'], solid_capstyle='butt')
    plt.plot([3.5,4.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m120'], solid_capstyle='butt')
    plt.plot([4.5,5.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['m192'], solid_capstyle='butt')

    plt.plot([-.5+6,.5+6], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c48'], solid_capstyle='butt')
    plt.plot([.5+6,1.5+6], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c60'], solid_capstyle='butt')
    plt.plot([1.5+6,2.5+6], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c72'], solid_capstyle='butt')
    plt.plot([2.5+6,3.5+6], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c96'], solid_capstyle='butt')
    plt.plot([3.5+6,4.5+6], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c120'], solid_capstyle='butt')
    plt.plot([4.5+6,4.5+6], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = gC['c192'], solid_capstyle='butt')

    plt.text(0, -.8, '2', color = gC['m48'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(1, -.8, '2.5', color = gC['m60'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(2, -.8, '3', color = gC['m72'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(3, -.8, '4', color = gC['m96'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(4, -.8, '5', color = gC['m120'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(5, -.8, '8', color = gC['m192'], horizontalalignment='center', fontproperties=fontTicks)

    plt.text(0+6, -.8, '2', color = gC['c48'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(1+6, -.8, '2.5', color = gC['c60'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(2+6, -.8, '3', color = gC['c72'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(3+6, -.8, '4', color = gC['c96'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(4+6, -.8, '5', color = gC['c120'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(5+6, -.8, '8', color = gC['c192'], horizontalalignment='center', fontproperties=fontTicks)
    

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

    #rod Colors
    rC = {
    'r' : '#747474',
    'c': '#dac910',
    'plt' : '',
    }

    rC['plt']=[
    rC['r'],rC['r'],rC['r'],rC['r'],
    rC['c'],rC['c'],rC['c'],rC['c'],
    ]

    if data.shape[0]==0:
        new_data = pd.DataFrame(np.zeros([1,data.shape[1]]), columns=data.columns)
        data = data.append(new_data)

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

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
        
    ax.axvline(x = -.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 3.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')
    ax.axvline(x = 7.5, color = 'white', linewidth = 2, alpha = 1, solid_capstyle='butt')

    n = 4
    plt.plot([-.5,3.5], [-.5,-.5], '-', lw=5, color = rC['r'], solid_capstyle='butt')    
    plt.plot([-.5+n,3.5+n], [-.5,-.5], '-', lw=5, color = rC['c'], solid_capstyle='butt')
    
    plt.plot([-.5,3.5], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = rC['r'], solid_capstyle='butt')
    plt.plot([-.5+n,3.5+n], [data.shape[0]-.5,data.shape[0]-.5], '-', lw=5, color = rC['c'], solid_capstyle='butt')
    
    plt.text(1.5, -.8, 'rods', color = rC['r'], horizontalalignment='center', fontproperties=fontTicks)
    plt.text(1.5+n, -.8, '!rods', color = rC['c'], horizontalalignment='center', fontproperties=fontTicks)
    

    return im, cbar