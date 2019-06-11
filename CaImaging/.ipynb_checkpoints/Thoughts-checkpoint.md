# Calcium Image Analysis thoughts:
#### Angueyra (June 2019)

1. Use ScanImageTiffReader to load SciScan tiff stack into numpy array:  
	```python
	from ScanImageTiffReader import ScanImageTiffReader
	tif = ScanImageTiffReader('test16bit.tif')
	im = tif.data(beg=0, end=800)
	s = tif.shape()
	tif.close()
	s
	```  
	
1. Use scikit to rescale stack (https://scikit-image.org/docs/dev/user_guide):  
	```python
	 from skimage import exposure
	 image = exposure.rescale_intensity(img10bit, in_range=(0, 2**10 - 1))
	 ```  
	 or  
	 ```python
	 image = exposure.rescale_intensity(img10bit, in_range='uint10') 
	 ```  
	 can also accept out-range, or assumes that rescaling should be done to image dtype  
1. Save as new tiff stack
```python
???
```
    * ImageDescription required some elements (https://docs.openmicroscopy.org/ome-model/5.6.3/specifications/index.html):
    ```xml
    <?xml version="1.0"?>
    <OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xmlns:str="http://exslt.org/strings"
         xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06                              http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">

      <Instrument ID="Instrument:0">
        <Detector ID="Detector:0:0" Model="COOLSNAP_HQ / ICX285" Type="CCD"/>
        <Objective ID="Objective:10002" Immersion="Oil" LensNA="1.4" Manufacturer="Olympus" NominalMagnification="100"/>
      </Instrument>
      <Image ID="Image:1" Name="example_R3D_D3D.dv">
        <AcquisitionDate>2005-01-28T13:50:08</AcquisitionDate>
        <Description>An example OME compliant file, 
                based on a wide-field microscope image</Description>
        <ObjectiveSettings ID="Objective:10002" Medium="Oil" RefractiveIndex="1.52" CorrectionCollar="7"/>
        <Pixels DimensionOrder="XYCZT" ID="Pixels:1" PhysicalSizeX="0.06631" PhysicalSizeY="0.06631" PhysicalSizeZ="0.2" SizeC="3" SizeT="1" SizeX="480" SizeY="480" SizeZ="5" Type="int16">
          <Channel EmissionWavelength="457" ExcitationWavelength="360" ID="Channel:1:0" NDFilter="0.5" Name="DAPI" SamplesPerPixel="1" Fluor="DAPI" IlluminationType="Epifluorescence" ContrastMethod="Fluorescence" AcquisitionMode="WideField" Color="65535">
            <DetectorSettings Binning="1x1" Gain="0.5" ID="Detector:0:0" ReadOutRate="10.0"/>
          </Channel>
          <Channel EmissionWavelength="528" ExcitationWavelength="490" ID="Channel:1:1" NDFilter="0.0" Name="FITC" SamplesPerPixel="1" Fluor="GFP" IlluminationType="Epifluorescence" ContrastMethod="Fluorescence" AcquisitionMode="WideField" Color="16711935">
            <DetectorSettings Binning="1x1" Gain="0.5" ID="Detector:0:0" ReadOutRate="10.0"/>
          </Channel>
          <Channel EmissionWavelength="617" ExcitationWavelength="555" ID="Channel:1:2" NDFilter="0.0" Name="RD-TR-PE" SamplesPerPixel="1" Fluor="TRITC" IlluminationType="Epifluorescence" ContrastMethod="Fluorescence" AcquisitionMode="WideField" Color="-16776961">
            <DetectorSettings Binning="1x1" Gain="0.5" ID="Detector:0:0" ReadOutRate="10.0"/>
          </Channel>
          <BinData BigEndian="false" Length="0"/>
          <Plane DeltaT="0.0" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.496" TheC="0" TheT="0" TheZ="0"/>
          <Plane DeltaT="0.294" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.696" TheC="0" TheT="0" TheZ="1"/>
          <Plane DeltaT="0.587" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.896" TheC="0" TheT="0" TheZ="2"/>
          <Plane DeltaT="0.881" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-22.096" TheC="0" TheT="0" TheZ="3"/>
          <Plane DeltaT="1.174" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-22.296" TheC="0" TheT="0" TheZ="4"/>
          <Plane DeltaT="9.625" ExposureTime="0.3" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.496" TheC="1" TheT="0" TheZ="0"/>
          <Plane DeltaT="10.12" ExposureTime="0.3" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.696" TheC="1" TheT="0" TheZ="1"/>
          <Plane DeltaT="10.613" ExposureTime="0.3" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.896" TheC="1" TheT="0" TheZ="2"/>
          <Plane DeltaT="11.106" ExposureTime="0.3" PositionX="3316.37" PositionY="-646.46" PositionZ="-22.096" TheC="1" TheT="0" TheZ="3"/>
          <Plane DeltaT="11.599" ExposureTime="0.3" PositionX="3316.37" PositionY="-646.46" PositionZ="-22.296" TheC="1" TheT="0" TheZ="4"/>
          <Plane DeltaT="25.447" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.496" TheC="2" TheT="0" TheZ="0"/>
          <Plane DeltaT="25.739" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.696" TheC="2" TheT="0" TheZ="1"/>
          <Plane DeltaT="26.033" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-21.896" TheC="2" TheT="0" TheZ="2"/>
          <Plane DeltaT="26.326" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-22.096" TheC="2" TheT="0" TheZ="3"/>
          <Plane DeltaT="26.619" ExposureTime="0.1" PositionX="3316.37" PositionY="-646.46" PositionZ="-22.296" TheC="2" TheT="0" TheZ="4"/>
        </Pixels>
      </Image>        
    </OME>
    ```
1. Use suite2p for image registration (https://github.com/MouseLand/suite2p/wiki/Registration) 
	* First find target reference image  
	```python
	from suite2p import register
	refImg = register.pick_init(ops)
	```  
	where ops requires: nonrigid=False, num_workers, and maxregshift
	* Then run rigid registration using reference:  
	```python
	maskMul,maskOffset,cfRefImg = register.prepare_masks(refImg)
refAndMasks = [maskMul,maskOffset,cfRefImg]
aligned_data, yshift, xshift, corrXY, yxnr = register.phasecorr(data, refAndMasks, ops)	
	```  
	where ops needs the reg_file, Ly = data.shape[1], Lx= data.shape[2], and nimg_init parameters  
	* Replace tif data with aligned_data and resave tiff
	* Save rest of things in new file (since there is some random sampling when calculating refImg, probably don't want to run registration repeatedly)
1. Could run suite2p here without registration
	* Could also try to just follow suite2p functions to save registration parameters.
1. calculate DSI /PD on a pixel by pixel basis:
  * findPeaks? watersheding?
  * required inputs are prepts, stmpts, nDirections (extract directly from hdf5 file in python? go thtough Ovation/Matlab and make csv? match based on timestamps?)
  * this should be parallelized
1. display results as image
	* make/find a circular color map to display PD
	* display DSI as grayscale  
	```python
	from skimage import data
	from skimage.viewer import ImageViewer

	image = data.coins()
	viewer = ImageViewer(image)
	viewer.show()
	```  
	* images could be combined, especially if colormap is isoluminant. Maybe be easier with just 2 colors (red vs. blue?).
1. try image segmentation algorithms in this image
1. For Masks, skimage has easy invert:
	```python
	from skimage import util
	img = data.camera()
	inverted_img = util.invert(img)
	```  
	
---

---
Good resources:
https://scikit-image.org/docs/dev/user_guide
https://scikit-image.org/docs/dev/user_guide/tutorial_segmentation.html
https://scikit-image.org/docs/dev/user_guide/tutorial_parallelization.html
https://scikit-image.org/docs/dev/user_guide/viewer.html
https://scikit-image.org/docs/dev/auto_examples/index.html#examples-gallery   

---