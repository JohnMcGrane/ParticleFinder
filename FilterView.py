"""*******************************************************************************
 *	
 *  Author: John M. Fox
 *  Institution: University of Minnesota Duluth
 *  College: Swenson College of Science and Engineering
 *  Department: Chemistry & Biochemistry
 *  Copyright: John Fox, 2021
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
    	um 		   	-- the stepsize, in microns, between consecutive spectra in the field of view (default = 50)
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


















