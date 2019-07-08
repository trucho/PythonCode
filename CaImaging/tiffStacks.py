import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os

from ScanImageTiffReader import ScanImageTiffReader as sitr
from tifffile import TiffFile, imwrite, xml2dict
from skimage.viewer import ImageViewer
import skimage
import cmocean


from suite2p.run_s2p import run_s2p


class sciscanTiff:
    """Python loader for tiff stack obtained using galvo2P SciScan"""
    """Created May_2019 (Angueyra)"""
    """Rescales 12-bit images from SciScan into 16-bit images, runs registration using suite 2P and saves a new registered, rescaled stack"""
    def __init__(self, filepath, filename):
        ### directories and paths ###
        self.basedir = '/Users/angueyraaristjm/Documents/LiImaging/TwoPhoton/'
        self.filepath = filepath + '/'
        self.filename = filename 
        self.savename_rs = 'rs_' + filename # path for rescaled tiff
        self.savename_rf = 'rf_' + filename # path for reference tiff
        self.savename_rg = 'rg_' + filename # path for registered tiff
        self.savepath = self.basedir + self.filepath + "analysis/"
        self.loadpath = self.basedir + self.filepath + self.filename + ".tif"
        self.rs_file = self.savepath + self.savename_rs + '.tif'
        self.rf_file = self.savepath + self.savename_rf + '.tif'
        self.rg_file = self.savepath + self.savename_rg + '.tif'
        self.s2p_regpath = self.savepath + 'suite2p/plane0/reg_tif/'
        ### Load tif file
        self.tif = sitr(self.loadpath)
        ### Run default methods
        self.metadata = self.getMetadata()
        
    def runRegistration(self):
        self.saveRescaledTiff();
        print('\nsuite2P start:\n')
        self.s2p_Run();
        print('\nsuite2P end:\n')
        self.s2p_saveRegStack();


    # ANALYSIS
    def s2p_Run(self):
        # create ops dictionary
        ops, db = self.s2p_ops();
        #run suite2P
        opsEnd = run_s2p(ops=ops,db=db);
        # save reference Image
        refImg = opsEnd.item()['refImg']
        refImg[refImg<0] = 0;
        refImg = refImg * 2; #rescaling to 16-bit
        refImg = np.ndarray.astype(refImg,'uint16')
        skimage.external.tifffile.imsave(self.rf_file, refImg, compress=0, metadata=None)
        return opsEnd
    
    def s2p_saveRegStack(self):
        rg_fnames = np.sort([f for f in os.listdir(self.s2p_regpath) if f.endswith('.tif')])
        if np.size(rg_fnames)>0:
            onceFlag = False
            for rT in rg_fnames:
                rgTif = sitr(self.s2p_regpath + rT)
                if onceFlag:
                    rgData = np.concatenate((rgData, rgTif.data()),axis=0)
                else:
                    onceFlag = True
                    rgData = rgTif.data();
            rgData[rgData<0] = 0;
            rgData = rgData * 2; #rescaling to 16-bit
            rgData = np.ndarray.astype(rgData,'uint16')
            skimage.external.tifffile.imsave(self.rg_file, rgData, compress=0, description=self.tif.description(0), metadata=None)
            print('Saved registered stack: ' + self.rg_file)
        elif np.size(rg_fnames)==0:
            print('Could not find .tif in ' + self.s2p_regpath)
    
    # rescale images from 11 bits to 16 bits
    def rescaleSciScanData(self):
        rsData = self.getData();
        rsData[rsData<np.power(2,15)] = np.power(2,15) # remove spurious negative values
        rsData = rsData - np.power(2,15) # subtract baseline
        rsData = rsData * 32 # multiply by 2^16 - 2^11
        return rsData
    
    def saveRescaledTiff(self):
        if not os.path.exists(self.savepath):
            os.mkdir(self.savepath)
        if not os.path.exists(self.rs_file):
            rsData = self.rescaleSciScanData();
            skimage.external.tifffile.imsave(self.rs_file, rsData, compress=0, description=self.tif.description(0), metadata=None)
            print('Saved rescaled stack:' + self.rs_file)
        else:
            print('Stack already rescaled:' + self.rs_file)
    
    # CONVENIENCE
    def getMetadata(self):
        imgDescription = self.getImageDescription()
        metadata = {
            'Lt' : imgDescription['OME']['Image']['Pixels']['SizeT'],
            'Lx' : imgDescription['OME']['Image']['Pixels']['SizeX'],
            'Ly' : imgDescription['OME']['Image']['Pixels']['SizeY'],
            'nChannels' :  imgDescription['OME']['Image']['Pixels']['SizeC'],
            'realX' : np.multiply(imgDescription['OME']['StructuredAnnotations']['XMLAnnotation'][3]['Value']['ImagePhysicalDimensions']['PhysicalSizeX'],1e6), # in um
            'realY' : np.multiply(imgDescription['OME']['StructuredAnnotations']['XMLAnnotation'][3]['Value']['ImagePhysicalDimensions']['PhysicalSizeY'],1e6), # in um
            'realT' : imgDescription['OME']['StructuredAnnotations']['XMLAnnotation'][3]['Value']['ImagePhysicalDimensions']['PhysicalSizeT'], # in s
            'samplingRate' : 0,
            'cellBodyDiameterX': 0,
            'cellBodyDiameterY': 0,
        }
        metadata['samplingRate'] = np.divide(metadata['Lt'],metadata['realT'])
        metadata['cellBodyDiameterX'] = np.divide(metadata['Lx'],metadata['realX']) * 5
        metadata['cellBodyDiameterY'] = np.divide(metadata['Ly'],metadata['realY']) * 5
        return metadata
    
    def getImageDescription(self):
        imgDescription = xml2dict(self.tif.description(0))
        return imgDescription
    
    def getData(self):
        img = self.tif.data();
        return img
    
    
    # Suite2P specific
    def s2p_ops(self):
        ops = {
            'fast_disk': [], # used to store temporary binary file, defaults to save_path0 (set as a string NOT a list)
            'save_path0': [], # stores results, defaults to first item in data_path
            'delete_bin': False, # whether to delete binary file after processing
            # main settings
            'nplanes' : 1, # each tiff has these many planes in sequence
            'nchannels' : 1, # each tiff has these many channels per plane
            'functional_chan' : 1, # this channel is used to extract functional ROIs (1-based)
            'diameter': [self.metadata['cellBodyDiameterX'], self.metadata['cellBodyDiameterY']], # this is the main parameter for cell detection, 2-dimensional if Y and X are different (e.g. [6 12])
            'tau':  2.0, # this is the main parameter for deconvolution (GCaMP6f = 0.7; GCaMP6m = 1.25; GCaMP6s = 2.0)
            'fs': self.metadata['samplingRate'],  # sampling rate (total across planes) (128 x 128 = 12.2; 256 x 128 = 6.1)
            # output settings
            'save_mat': False, # whether to save output as matlab files
            'combined': True, # combine multiple planes into a single result /single canvas for GUI
            # parallel settings
            'num_workers': 0, # 0 to select num_cores, -1 to disable parallelism, N to enforce value
            'num_workers_roi': -1, # 0 to select number of planes, -1 to disable parallelism, N to enforce value
            # registration settings
            'do_registration': True, # whether to register data
            'nimg_init': 200, # subsampled frames for finding reference image
            'batch_size': 200, # number of frames per batch
            'maxregshift': 0.1, # max allowed registration shift, as a fraction of frame max(width and height)
            'align_by_chan' : 1, # when multi-channel, you can align by non-functional channel (1-based)
            'reg_tif': True, # whether to save registered tiffs
            'subpixel' : 10, # precision of subpixel registration (1/subpixel steps)
            'nonrigid': False, # wheter to perform non-rigid registration
            # cell detection settings
            'roidetect': False, #whether of not to run cell-detection algorithm
            'connected': False, # whether or not to keep ROIs fully connected (set to 0 for dendrites)
            'navg_frames_svd': 5000, # max number of binned frames for the SVD
            'nsvd_for_roi': 1000, # max number of SVD components to keep for ROI detection
            'max_iterations': 20, # maximum number of iterations to do cell detection
            'ratio_neuropil': 6., # ratio between neuropil basis size and cell radius
            'ratio_neuropil_to_cell': 3, # minimum ratio between neuropil radius and cell radius
            'tile_factor': 1., # use finer (>1) or coarser (<1) tiles for neuropil estimation during cell detection
            'threshold_scaling': 1., # adjust the automatically determined threshold by this scalar multiplier
            'max_overlap': 0.75, # cells with more overlap than this get removed during triage, before refinement
            'inner_neuropil_radius': 2, # number of pixels to keep between ROI and neuropil donut
            'outer_neuropil_radius': np.inf, # maximum neuropil radius
            'min_neuropil_pixels': 350, # minimum number of pixels in the neuropil
            # deconvolution settings
            'baseline': 'maximin', # baselining mode
            'win_baseline': 60., # window for maximin
            'sig_baseline': 10., # smoothing constant for gaussian filter
            'prctile_baseline': 8.,# optional (whether to use a percentile baseline)
            'neucoeff': .7,  # neuropil coefficient
        }
        db = {
            'h5py': [], # a single h5 file path
            'h5py_key': 'data',
            'look_one_level_down': False, # whether to look in ALL subfolders when searching for tiffs
            'data_path': [self.savepath], # a list of folders with tiffs 
                                                             # (or folder of folders with tiffs if look_one_level_down is True, or subfolders is not empty)
            'subfolders': [], # choose subfolders of 'data_path' to look in (optional)
            'fast_disk': [], # string which specifies where the binary file will be stored (should be an SSD)
        }
        return ops, db
    
        # def s2p_rgOps(self):
        #     ops = {
        #         'reg_file': self.rs_file,
        #         'Ly': self.metadata['Ly'],
        #         'Lx': self.metadata['Lx'],
        #         'nplanes' : self.metadata['Lt'],
        #         'nchannels' : 1, # each tiff has these many channels per plane
        #         'align_by_chan' : 1, # when multi-channel, you can align by non-functional channel (1-based)
        #         'batch_size': 200, # number of frames per batch
        #         'maxregshift': 0.1, # max allowed registration shift, as a fraction of frame max(width and height)
        #         'nimg_init': 200,
        #         'subpixel' : 10, # precision of subpixel registration (1/subpixel steps)
        #         'smooth_sigma': 1.0, # ~1 good for 2P recordings, recommend >5 for 1P recordings
        #         'th_badframes': 1.0, # determines frames to exclude when determining cropping - set it smaller to exclude more frames
        #         'pad_fft': False,
        #         'do_phasecorr': True,
        #         'nonrigid': False,
        #         'num_workers': 0, # 0 to select num_cores, -1 to disable parallelism, N to enforce value
        #         'non-rigid': False,
        #
        #
        #         'fast_disk': [], # used to store temporary binary file, defaults to save_path0 (set as a string NOT a list)
        #         'save_path0':  self.savepath + 'suite2p/', # stores results, defaults to first item in data_path
        #         'delete_bin': False, # whether to delete binary file after processing
        #     }
        #     return ops
    

    # DISPLAY
    
    def histRescaled(self):
        rsData = self.rescaleSciScanData();
        plotHistogram(rsData);
    
    def histRescaled_single(self,t=0):
        rsData = self.rescaleSciScanData();
        if t>=0 & t<np.shape(rsData)[0]:
            rsData = rsData[t,:,:];
        else:
            t=0
            print('Showing first plane; t is out-of-bounds')
        plotHistogram(rsData);
            

            
"""Common methods"""
def getMean(tiffData):
    meanData = np.mean(tiffData,axis=0)
    return meanData

def getTSeries(maskedData):
    tData = np.mean(np.mean(maskedData,axis=1), axis=1)
    return tData

def norm16bit(tiffData):
    nData = tiffData
    nData = nData - np.min(nData)
    nData = (nData)/np.max(nData)
    nData = nData * np.power(2,16)
    return nData
    
def getSD(tiffData):
    sdData = np.std(tiffData,axis=0)
    return sdData

def getMIP(tiffData):
    mipData = np.max(tiffData,axis=0)
    return mipData

def calcDS(peaks,angles):
    dsR, dsT = vectorSum(peaks,angles);
    return dsR, dsT

def cart2pol(x,y):
    rho = np.sqrt(np.power(x,2)+np.power(y,2))
    theta = np.arctan2(y,x)
    return(rho,theta)

def pol2cart(rho, theta):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return(x, y)

def vectorSum(rho, theta):
    x,y = pol2cart(rho,theta)
    return(cart2pol(np.sum(x),np.sum(y)))

def calc_corrcoeff(tiffData, stim):
    nt = np.shape(tiffData)[0];
    ny = np.shape(tiffData)[1];
    nx = np.shape(tiffData)[2];
    if nt != np.size(stim):
        print('input dimensions do no match: tiffData_t = %f vs. stim = %f',nt, np.size(stim))
    else:
        ccData = np.zeros(ny,nx)
        for y in range(ny):
            for x in range(nx):
                ccData[y,x] = np.max(signal.correlate(tiffData[:,y,x],stim, mode = 'full', method ='fft'));
        return ccData

def autoThresh(tiffData):
    threshold = skimage.filters.threshold_mean(getMean(tiffData))
    thMask = getMean(tiffData) >= threshold
    return thMask        

def plotMean(tiffData):
    fig, ax0, ax1, ax2 = plotHistogram(getMean(tiffData));
    return fig, ax0, ax1, ax2

def plotSD(tiffData):
    fig, ax0, ax1, ax2 = plotHistogram(getSD(tiffData));
    return fig, ax0, ax1, ax2

def plotMIP(tiffData):
    fig, ax0, ax1, ax2 = plotHistogram(getMIP(tiffData));
    return fig, ax0, ax1, ax2
       
def plotHistogram(tiffData):
    if np.ndim(tiffData)==2:
        imgData = tiffData
    elif np.ndim(tiffData)==3:
        imgData = getMean(tiffData)
    # [hist, bin_edges] = np.histogram(tiffData[tiffData>0], bins='fd', density=True) #ignore black pixels
    # bins = (bin_edges[:-1] + bin_edges[1:]) / 2
    [hist, bins] = skimage.exposure.histogram(tiffData[tiffData>0], normalize=True) #ignore black pixels
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(18, 9), gridspec_kw={'width_ratios': [1, .1, 1]})

    plt.sca(ax0);
    plt.imshow(imgData, cmap=plt.cm.gray, interpolation='nearest')
    plt.xticks([]); plt.yticks([]); plt.colorbar();
    
    plt.sca(ax1);
    plt.plot([0,0], [0,np.size(tiffData[tiffData==0])/np.size(tiffData)], lw=2)
    plt.xticks([0])
    ax1.spines['top'].set_visible(False); ax1.spines['right'].set_visible(False)
    
    plt.sca(ax2);
    plt.plot(bins, hist, lw=2)
    ax2.spines['top'].set_visible(False); ax2.spines['right'].set_visible(False)
    return fig, ax0, ax1, ax2

def phaseColormap():
    return cmocean.cm.phase
    