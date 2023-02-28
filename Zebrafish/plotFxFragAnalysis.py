from cmcrameri import cm #colormaps
import numpy as np
from numpy import floor
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

batlow = cm.vik(range(255))
batlow = batlow[0:255:int(floor(255/8)),:]

def formatPeaks(figH, axH, plotname):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=6)
    fontLabels = font_manager.FontProperties(fname=font_path, size=9)
    fontTitle = font_manager.FontProperties(fname=font_path, size=9)
    
    axH.set_xlabel('size (bp)', fontproperties=fontLabels)
    axH.set_ylabel('fluo (au)', fontproperties=fontLabels)

    axH.set_title(plotname, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    axH.yaxis.set_label_coords(-.14,.5)
    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    
def plotABIcsv(dirPath, fileName, gdf, gene):
    fH, axH = plt.subplots(figsize=[4,2])
    fH.subplots_adjust(left=0.16, bottom=0.2, wspace=0.1)

    lim1, lim2 = getLims(gene)
    ilim1 = np.searchsorted(gdf['size'],lim1)
    ilim2 = np.searchsorted(gdf['size'],lim2)
    pH = axH.plot(gdf['size'][ilim1:ilim2],gdf['fluo'][ilim1:ilim2], color='#4669F2', zorder = 1)
    formatPeaks(fH, axH,fileName)
    return fH, axH
    
def formatPeaks_withladder(figH, axH1, axH2, plotname):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=14)
    fontLabels = font_manager.FontProperties(fname=font_path, size=18)
    fontTitle = font_manager.FontProperties(fname=font_path, size=18)
    
    axH2.set_xlabel('size (bp)', fontproperties=fontLabels)
    axH1.set_ylabel('fluo (au)', fontproperties=fontLabels)
    axH2.set_ylabel('fluo (au)', fontproperties=fontLabels)

    axH1.set_title(plotname, fontproperties=fontTitle)
    axH1.get_xaxis().set_visible(False)
    axH1.spines['bottom'].set_visible(False)
    axH1.spines['top'].set_visible(False)
    axH1.spines['right'].set_visible(False)
    axH2.spines['top'].set_visible(False)
    axH2.spines['right'].set_visible(False)

    for label in (axH1.get_xticklabels() + axH1.get_yticklabels()+ axH2.get_xticklabels() + axH2.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    
    figH.align_ylabels()

def plotABIcsv_withladder(dirPath, fileName, gdf, gene):
    fH = plt.figure(constrained_layout=True, figsize=[8,4])
    gs = fH.add_gridspec(12, 1)
    ax1 = fH.add_subplot(gs[0:3, 0])
    ax2 = fH.add_subplot(gs[4:, 0], sharex = ax1)
    
    lim1, lim2 = getLims(gene)
    ilim1 = np.searchsorted(gdf['size'],lim1)
    ilim2 = np.searchsorted(gdf['size'],lim2)
    pH = ax1.plot(gdf['size'][ilim1:ilim2],gdf['ladder'][ilim1:ilim2], color='#CC2C2A', zorder = 1)
    pH = ax2.plot(gdf['size'][ilim1:ilim2],gdf['fluo'][ilim1:ilim2], color='#4669F2', zorder = 1)
    formatPeaks(fH, ax1,ax2,fileName)
    return fH, [ax1,ax2]

def getLims(gene):
    if gene == 'eml1': lim1 = 150; lim2 = 300;
    elif gene == 'sema7a': lim1 = 250; lim2 = 450;
    elif gene == 'gnat2': lim1 = 200; lim2 = 350;
    elif gene == 'syt5a': lim1 = 400; lim2 = 550;
    elif gene == 'efna1b': lim1 = 400; lim2 = 550;
    elif gene == 'tbx2a': lim1 = 400; lim2 = 550; # lim1 = 350; lim2 = 650;
    elif gene == 'tbx2aFiiRii': lim1 = 400; lim2 = 550;
#     elif gene == 'tbx2a': lim1 = 450; lim2 = 600;
#     elif gene == 'tbx2b': lim1 = 250; lim2 = 400;
    elif gene == 'tbx2b': lim1 = 280; lim2 = 380;
#     elif gene == 'ntng2b': lim1 = 100; lim2 = 300;
    # elif gene == 'foxq2': lim1 = 350; lim2 = 500; # for F0
    elif gene == 'foxq2': lim1 = 410; lim2 = 440; # for Laura's -4bp mutants
    elif gene == 'ntng2b': lim1 = 200; lim2 = 300;
    elif gene == 'syt5b': lim1 = 250; lim2 = 400;
    elif gene == 'nr2f1b': lim1 = 300; lim2 = 500;
    elif gene == 'xbp1': lim1 = 150; lim2 = 400;
    elif gene == 'lrrfip1a': lim1 = 200; lim2 = 450;
    elif gene == 'skor1a': lim1 = 250; lim2 = 450;
    elif gene == 'sall1a': lim1 = 100; lim2 = 250;
    elif gene == 'tgif': lim1 = 200; lim2 = 350;
    elif gene == 'lhx1a': lim1 = 300; lim2 = 450;
    elif gene == 'tefa': lim1 = 250; lim2 = 400;
    elif gene == 'empty': lim1 = 100; lim2 = 600;
    else: lim1 = 100; lim2 = 600;
    return lim1, lim2


def plotABIfig(dirPath, fileName, gdf, gene):
    fH, axH = plt.subplots(figsize=[11,8])
    fH.subplots_adjust(left=0.16, bottom=0.2, wspace=0.1)

    lim1, lim2 = getLims(gene)
    ilim1 = np.searchsorted(gdf['size'],lim1)
    ilim2 = np.searchsorted(gdf['size'],lim2)
    pH = axH.plot(gdf['size'][ilim1:ilim2],gdf['fluo'][ilim1:ilim2]/100, color='#4669F2', zorder = 1, linewidth=4)
#     formatPeaksFig(fH, axH,filename)
    formatPeaksFig(fH, axH,'')
    return fH, axH

def formatPeaksFig(figH, axH, plotname):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=50)
    fontLabels = font_manager.FontProperties(fname=font_path, size=50)
    fontTitle = font_manager.FontProperties(fname=font_path, size=24)
    
    axH.set_xlabel('size (bp)', fontproperties=fontLabels)
    axH.set_ylabel('x 10$^2$ fluo (au)', fontproperties=fontLabels)
    
    plt.locator_params(axis='x', nbins=6)
    plt.locator_params(axis='y', nbins=2)
    axH.set_title(plotname, fontproperties=fontTitle)
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    axH.spines['bottom'].set_linewidth(4)
    axH.spines['left'].set_linewidth(4)
    axH.yaxis.set_label_coords(-.12,.48)
    axH.xaxis.set_tick_params(width=3, length=10)
    axH.yaxis.set_tick_params(width=3, length=10)
    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)

print('defs ready!')