from cmcrameri import cm #colormaps
from numpy import floor
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

batlow = cm.vik(range(255))
batlow = batlow[0:255:int(floor(255/8)),:]

def applyPlotStyle(plotStyleString):
    if plotStyleString=='Dark':
        # dark background
        params = {"ytick.color" : "w",
                  "xtick.color" : "w",
                  "axes.labelcolor" : "w",
                  "axes.edgecolor" : "w",
                 "axes.linewidth" : 3,
                 "xtick.major.width" : 3,
                 "ytick.major.width" : 3,
                 "xtick.major.size" : 8,
                 "ytick.major.size" : 8,
                 "text.color" : "w"}
        plt.rcParams.update(params)
        plt.style.use('dark_background')
    elif plotStyleString=='Light':
        # white background
        params = {"ytick.color" : "k",
                  "xtick.color" : "k",
                  "axes.labelcolor" : "k",
                  "axes.edgecolor" : "k",
                 "axes.linewidth" : 3,
                 "xtick.major.width" : 3,
                 "ytick.major.width" : 3,
                 "xtick.major.size" : 8,
                 "ytick.major.size" : 8,
                 "text.color" : "k"}
    plt.rcParams.update(params)
    font_prop = font_manager.FontProperties(fname='/System/Library/Fonts/Avenir.ttc')
    matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=batlow)
    print('Plotting style is ' + plotStyleString)

#gene Colors
zfC = {
    'R' : '#7d7d7d',
    'U' : '#B73AB9',
    'S' : '#4364F6',
    'M' : '#59CB3B',
    'L' : '#CE2A22',
    'rho'  : '#7d7d7d',
    'sws1' : '#B73AB9',
    'sws2' : '#4364F6',
    'mws1' : '#59CB3B',
    'mws2' : '#59CB3B',
    'mws3' : '#59CB3B',
    'mws4' : '#59CB3B',
    'lws1' : '#CE2A22',
    'lws2' : '#CE2A22',
    'actb2': '#926645',
    'tbx2a': '#c92675',
    'tbx2b': '#7526c9',
    'six7' : '#d6ab00',
}

zfG = {
    'wt' : '#000000',
    'tbx2a' : '#ab266b',
    'tbx2b' : '#421f8e',
    'foxq2' : '#001dd6',
    'nr2e3' : '#7d7d7d',
    'lrrfip1a' : '#00CC6A',
    'xbp1' : '#B73AB9',
    'sall1a' : '#B89504',
}

zfGm = {
    'wt' : 'o',
    'tbx2a' : 'P',
    'tbx2b' : 'X',
    'foxq2' : '^',
    'nr2e3' : '+',
    'lrrfip1a' : 'D',
}

prLabel = {
    'R'  : 'Rods',
    'U' : 'UV',
    'S' : 'S',
    'M' : 'M',
    'L' : 'L',
}

def formatFigureMain(figH, axH, plotH):
    font_path = '/System/Library/Fonts/Avenir.ttc'
    fontTicks = font_manager.FontProperties(fname=font_path, size=30)
    fontLabels = font_manager.FontProperties(fname=font_path, size=36)
    fontTitle = font_manager.FontProperties(fname=font_path, size=32)
    axH.set_xscale('linear')
    axH.spines['top'].set_visible(False)
    axH.spines['right'].set_visible(False)
    
    for label in (axH.get_xticklabels() + axH.get_yticklabels()):
        label.set_fontproperties(fontTicks)
    axH.set_xlabel(axH.get_xlabel(), fontproperties = fontTicks)
    axH.set_ylabel(axH.get_ylabel(), fontproperties = fontTicks)
    return fontLabels

def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import colorsys
    try:
        c = matplotlib.colors.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*matplotlib.colors.to_rgb(c))
    return matplotlib.colors.rgb2hex(colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2]))

def estimateJitter(dataArray):
    """ creates random jitter scaled by local density of points"""
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(dataArray)
    density = kde(dataArray)
    jitter = np.random.randn(len(dataArray))*density
    return jitter

def formatFigure(figH, axH, plotH):
    fontLabels = formatFigureMain(figH, axH, plotH)
#     axH.set_xlabel('wt vs. cr', fontproperties=fontLabels)
    axH.set_ylabel('cells per 64 x 64 $\mu$m$^2$', fontproperties=fontLabels)
    axH.xaxis.set_tick_params(rotation=45)

def formatFigureRvU(figH, axH, plotH):
    fontLabels = formatFigureMain(figH, axH, plotH)
    axH.set_xlabel('Rods per 64 x 64 $\mu$m$^2$', fontproperties=fontLabels)
    axH.set_ylabel('UV cones per 64 x 64 $\mu$m$^2$', fontproperties=fontLabels)
    axH.xaxis.set_tick_params(rotation=45)
    
def formatFigureMvS(figH, axH, plotH):
    fontLabels = formatFigureMain(figH, axH, plotH)
    axH.set_xlabel('M cones per 64 x 64 $\mu$m$^2$', fontproperties=fontLabels)
    axH.set_ylabel('S cones per 64 x 64 $\mu$m$^2$', fontproperties=fontLabels)
    axH.xaxis.set_tick_params(rotation=45)