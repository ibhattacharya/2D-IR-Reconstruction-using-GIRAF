Matlab scripts for Reconstruction of 2 dimensional Infrared (2D IR) Spectroscopic Data 
======================================================================================

(c) *Ipshita Bhattacharya*, *Jonathan J. Humston*, *Christopher M. Chetaum*, *Mathews Jacon*, Nov 2017.

[Optics Letters]

**University of Iowa**

Code Structure
--------------
### Demo scripts
**main.m**
Reconstruction script to demonstrate GIRAF reconstruction from under-sampled time domian 2D-IR data. 

### Data 
**2D_IRSpec_GIRAF_data.mat**:
2D IR data from a sample of 50 mM sodium cyanate in methanol held in a sample cell with a 50 μm path
length.

**undersaMplingMask_US5.mat**
Example undersampling mask of factor 5

### Functions: 
**CreateAxesAll.m**: 
	Creates axes using calibration file.  
	
**LinearGriddedData.m**: 
	The w3 axis is slightly nonlinear. This code is to be used to interpolate the data onto a linearized grid.  
	
**GIRAF_preprocessing.m**: 
	Sets up GIRAF pre-processing steps, e.g. reindexing of data.  
	
**cosine_filtering.m**: 
	Processes raw data with cosine windowing.
	
**Hybrid_DownsamplingMASK.m**: 
	Create random undersampling mask.
	
**defAAt.m**: 
	Create forward and backward Fourier undersamling operator.
	
**runGIRAF.m**: 
	Main code executing GIRAF algorithm.
	
**RMSE.m**: 
	Calculates root mean square error of reconstructed data.
	
**left_fun.m**: 
	Solves the least square problem step of GIRAF reconstruction.
	
**cls_analysis.m**: 
	Calculates center line slope of the reconstructed data.
