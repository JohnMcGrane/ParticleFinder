"""*******************************************************************************
 *	
 *  Author: John M. Fox
 *  Institution: University of Minnesota Duluth
 *  College: Swenson College of Science and Engineering
 *  Department: Chemistry & Biochemistry
 *  Copyright: John Fox, 2021
 *
 *  Description: An API for the FilterView class. This class is designed to take 
 *  a directory of CSV files exported from Thermo Fisher Scientific's OMNIC 
 *  spectroscopy software. CSV files of every spectra in an OMNIC chemigram can
 *  be exported by selecting the "Split Map" option from the Atlus drop down menu.
 *  The exported file type can then be selected as CSV and then OMNIC will export
 *  all spectra in the chosen chemigram as CSV files to the selected destination 
 *  folder. This program performs a rough analysis on those csv spectra to determine 
 *  which likely originate from microplastic particles. This program assumes 
 *  that spectra that display enhanced reflectance (reduced absorbance) at 2850 
 *  wavenumbers, corresponding to the sp3 CH stretching region of the IR 
 *  spectrum, are likely to be microplastics. Spectra that meet the "enhanced
 *  reflectance" criterium are determined by subtracting the reflectance 
 *  intensity at 2850 wavenumbers from the reflectance intensity at 2750 
 *  wavenumbers (a proxy for the spectral baseline) for all spectra in the 
 *  given directory. Spectra that are outliers in the resulting Gaussian
 *  distribution of reflectance differences are assuemd to be spectra that were
 *  obtained from microplastic particles on the filter surface. The 
 *  functionality of the FilterView class and API is detailed below. For those 
 *  who are interested, the provided code framework could be expanded upon to 
 *  accomodate a more sophisticated form of spectral matching. A good place to 
 *  start would be in the diffList function and the pixelList function, both 
 *  located in the utilities file. A spectral matching algorithm was not 
 *  included in this version because spectral quality was not sufficient to
 *  perform spectral matching, therefore the more rudimentary approach based
 *  on reflectance differences was carried out, as discussed above. Of note, it 
 *  is acknowledged that both function implementations and program runtimes at 
 *  the moment have not been optimized. This program can accomodate all spectral
 *  ranges and assumes that the chemigramstep size matches the aperture x & y 
 *  dimensions, which should be the same (i.e. a square aperture).
 *
 ******************************************************************************"""

from utilities import *

class FilterView:

	def __init__(self,width,height,directory, um = 50, large = False, particles = []):
		"""Instantiate a FilterView object. 

		Arguments: 
		width 			-- the number of spectra along the x-axis of the field of view
		width type		-- int
		height 			-- the number of spectra along the y-axis of the field of view
		height type		-- int
		directory		-- the name of the directory containing all CSV files in the field of view 
		directory type	-- str

    	Keyword arguments:
    	um 		   		-- the stepsize, in microns, between consecutive spectra in the field of view (default = 50)
    	um type			-- int
    	large		 	-- if dealing with particles larger than 200 microns this should be set to True (default = False)
    	large type		-- bool
    	particles 		-- list of xy filter pixel coordinates found by manual spectral analysis (default = [])
    	particles type	-- list
    	"""
		self.width = width
		self.height = height
		self.name = directory
		self.stepsize = um
		self.large = large
		self.filelist = createFileList(self.name)
		self.manual = particles
		self.differencelist = diffList(self.filelist,self.name)
		self.pixellist = pixelList(self.filelist,self.differencelist,self.cutoffNumber())

	def cutoffNumber(self):
		"""Return a cutoff based on the distribution of reflectance differences for all spectra in the field of view.

    	Arguments:
    	None
    	"""
		return cutoffMaker(self.large,self.differencelist)

	def particleCount(self):
		"""Return a computationally determined microplastic particle count for the field of view.

    	Arguments:
    	None
    	"""
		particles = particleCount(self.pixellist,[self.width,self.height])[0]
		return particles

	def particleSize(self):
		"""Return a the average particle diameter assuming that all particles are circular.

    	Arguments:
    	None
    	"""
		groupList = particleCount(self.pixellist,[self.width,self.height])[1]
		sizeList = particleCount(self.pixellist,[self.width,self.height])[2]
		return sizer(groupList,sizeList,self.stepsize)

	def histogram(self):
		"""Creates a histogram of reflectance differences. Reflectance at 2750 cm-1 - reflectance at 2850 cm-1. 

    	Arguments:
    	None
    	"""
		if self.large == True:
			paredlist = toGauss(self.differencelist)
		elif self.large == False:
			paredlist = noGauss(self.differencelist)
		createHist(paredlist,self.differencelist,self.cutoffNumber())

	def showView(self):
		"""Display a screenshot of the filter if a screenshot is located in directory with the CSV files.

    	Arguments:
    	None
    	"""
		showScreenShot(self.name)

	def showChemigram(self):
		"""Display a visualization of the plastic particles identified by the computational analysis.
		
    	Arguments:
    	None
    	"""
		showChemigram(self.pixellist,[self.width,self.height])

	def visualCompare(self,startx,starty):
		"""Display a visualization comparing manually identified particles and computer identified particles.
		
    	Arguments:
    	startx		-- x coordinate of the leftmost pixel in the field of view
    	startx type -- int
    	starty		-- y coordinate of the lowest pixel in the field of view
    	starty type -- int
    	"""
		humanlist = self.manualPixels()
		likely = self.candidatePixels(startx,starty)
		humanpixels, comppixels = transformer(humanlist,likely,startx,starty,self.width,self.stepsize)
		printout(humanpixels,comppixels,[self.width,self.height])

	def textCompare(self,startx,starty):
		"""Display a text-based comparison of the manually identified particles and computer identified particles.
		Counts for both approaches are displayed as well as a count for particles that were identified by both techniques.
		
    	Arguments:
    	startx		-- x coordinate of the leftmost pixel in the field of view
    	startx type -- int
    	starty		-- y coordinate of the lowest pixel in the field of view
    	starty type -- int
    	"""
		humanlist = self.manualPixels()
		likely = self.candidatePixels(startx,starty)
		summary(humanlist,likely,self.particleCount(),self.manual)

	def candidatePixels(self,startx,starty):
		"""Display a list of xy pixel coordinates that the computer identified as plastics. 
		
    	Arguments:
    	startx		-- x coordinate of the leftmost pixel in the field of view
    	startx type -- int
    	starty		-- y coordinate of the lowest pixel in the field of view
    	starty type -- int
    	"""
		likely = candidates(self.pixellist,self.filelist,startx,starty,[self.width,self.height],self.stepsize)
		return likely

	def manualPixels(self):
		"""Display a list of xy pixel coordinates that the researcher identified as plastics. 
		
    	Arguments:
    	None

    	Note: Requires the user to pass in a list of xy coordinates using the "particles" keyword argument
    	"""
		manuallist = listcleaner(self.manual)
		return manuallist


if __name__ == '__main__':

	# The following script is designed to be run from the command line. 
	# User particles can be copied and pasted in to the command line to 
	# enable comparison of the computer program and the manual analysis. 
	file = input("Enter directory name: ")
	xdim = int(input("Enter x-dimension: "))
	ydim = int(input("Enter y-dimension: "))
	ask = input("Compare to manual analysis: (y/n)")
	if ask == 'y':
		print("Enter manually identified pixel coordinates: ")
		lines = []
		while True:
			line = input()
			if line:
				lines.append(line)
			else:
				break
		foundlist = '\n'.join(lines)
		startx = int(input("First x coordinate: "))
		starty = int(input("First y coordinate: "))
	else:
		foundlist = []
		startx = 0
		starty = 0


	view1 = FilterView(xdim,ydim,file,particles = foundlist)
	view1.showView()
	view1.showChemigram()
	view1.textCompare(startx,starty)
	print(f'Particle Size           : {view1.particleSize()}')
	print(f'Cutoff                  : {view1.cutoffNumber()}')
	print(f'Particle Count          : {view1.particleCount()}')
	view1.visualCompare(startx,starty)
	view1.histogram()
	print(view1.candidatePixels(startx,starty))
	print(view1.manualPixels())


















