January 28, 2017
NuclearMorphology.m

Code is designed to quantify nuclear and cellular shape factors. It performs nuclear and cellular segmentation from confocal .lsm images, and uses calcien AM segementation as a live/dead mask. 

Code from Dr. Amy Rowat's Lab
UCLA Department of Integrative Biology and Physiology
Los Angeles, CA

Image format:
File name format: '#.lsm'
LSM files have data in two color channels corresponding to:
  - R channel: Hoechst stain (nuclear)
  - G channel: Calcein AM (cellular/live stain)
    LSM files also contain images captured at 4 different locations in the
    wells.

Output:
  - creates output data folder in image parent directory
  - writes segmented image files
  - saves .mat file containing output data cell:
	- Row 1: cellular data
	- Row 2: nuclear data
	- Column Indx: corresponds to each image
	- Each cell contains a matrix of data:
		- Row Indx: each object
		- Column 1: Area
		- Column 2: Perimeter

Code originally by Kendra Nyberg (December 2016)
Code updated by Eliza King Lassman (January 2017)
  - Added bwselect feature to manually verify segmentation