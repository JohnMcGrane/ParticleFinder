import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from scipy.stats import norm
import pandas as pd
import os
import re
from scipy import ndimage
from scipy.stats.mstats import mquantiles
from PIL import Image 
import matplotlib.image as mpimg


def createFileList(name):
    """A function that takes a directory name and creates 
    a list of all the files in that directory"""
    filelist = []
    ct=0
    path = os.path.abspath(os.getcwd())
    for filename in os.listdir(f"{path}/{name}/"):
        if filename.endswith(".CSV"): 
            ct+=1
            filelist.append(os.path.basename(filename)) # Add .csv files to a list
    else:
        pass
    
    df = pd.DataFrame(filelist) # Add filenames to a pandas dataframe
    df.rename(columns = {df.columns[0]:"filename"},inplace=True) # Rename the column with the files to "filename"
    df = df.sort_values("filename",ignore_index = True, ascending = True) # Sort the files in the list
    return df

def cutoffMaker(booly,dlist):
    """Function that returns a cutoff value either immediately for small 
    particles or after using toGauss function for large particles"""
    if booly == False:
        thirdquartile = mquantiles(dlist)[2]
        iqr = stats.iqr(dlist)
        return thirdquartile + 1.5*iqr # Cutoff based on typical statistical outlier
    elif booly == True:
        paredlist = toGauss(dlist) # toGauss removes outliers then checks to remove further outliers
        thirdquartile = mquantiles(paredlist)[2] 
        iqr = stats.iqr(paredlist)
        return thirdquartile + 1.5*iqr
    else:
        print("Keyword argument \"largeParticles\" must be a Boolean")

def diffList(flist,name):
    """Function takes list of filenames and a directory name. For each .csv file
    the function calculates the difference between the reflectance at 2750 cm-1
    and 2850 cm-1 and returns a list of these differences. Based on .csv files with 
    any spectral range and any resolution"""
    numlistin=[]
    for j in range(flist.size):
        df1 = pd.read_csv(f"/Users/johnfox/Desktop/Microplastics/Code/{name}/{flist['filename'][j]}")
        df1.rename(columns = {df1.columns[0]:"wavenumber",df1.columns[1]:"reflectance"},inplace=True) # Rename columns
        if j == 0:
            i2850 = ((df1['wavenumber']-2850).abs().argsort()[:1])[0]
            i2750 = ((df1['wavenumber']-2750).abs().argsort()[:1])[0]
        
        num = df1.iloc[i2750][1] - df1.iloc[i2850][1]
        numlistin.append(num)
    return numlistin


def pixelList(flist,dlist,coff):
    """Function to return a list of pixels in the chemigram that have the 
    characteristic reflectance differences of microplastic particles. Takes 
    the list of .csv files, the list of reflectance differences, and a cutoff
    value."""
    pixellistin = []
    for j in range(flist.size):
        if dlist[j] > coff:
            pixellistin.append(j) 
    return pixellistin


def candidates(plist,flist,x,y,dimz,stepsize):
    """Function to return a list xy-coordinates that indicate which
    pixels are likely microplastics based on reflectance differences.
    This function takes the list of pixels with a positive reflectance difference,
    the whole list of .csv files, the x and y coordinates of the bottom left hand
    corner of the chemigram, and a tuple that indicates the (x,y)size in spectra 
    of the chemigram"""
    startx, starty = x,y
    state = np.zeros((dimz[0],dimz[1]),float) # Create a 41 x 41 matrix
    likely = []

    for istep in plist: # Iterate through the positive pixels and convert those into (x,y) coordinates
        ix1 = np.remainder((istep),dimz[0])
        iy1 = int(np.floor((istep)/(dimz[0])))
        likely.append((startx+ix1*stepsize,starty+iy1*stepsize))
    return likely



def listcleaner(newlist):
    """Function that returns list of tuples - (x,y) coordinates - found by 
    a researcher. This function takes an unformatted string of points in 
    order xy xy xy xy etc..."""
    if bool(newlist) == False:
        return []
    newlist = re.split('\s',newlist)
    foundlist = []
    if len(newlist) != 0:
        for k in range(0,len(newlist),2):
            foundlist.append((int(newlist[k]),int(newlist[k+1])))
    else:
        pass
    return foundlist


def updateState(plist,dimz,update = False, statein = 0,incrementby = 2):
    """Function that creates or updates a given state matrix by 2 unless 
    otherwise specified. This takes a list of pixels to update and the dimensions 
    of the chemigram. This can either update or create based on the update 
    keyword argument. The statein kwarg will be 0 unless specified. The incrementby
    kwarg will be 2 unless otherwise specified. Returns a state with the values of 
    positive pixels increased by two"""
    if update == False:
        state = np.zeros((dimz[0],dimz[1]),float)
    else:
        state = statein
    for istep in plist: # Iterate through the positive pixels and converts to xy coordinates
        ix1 = np.remainder((istep),dimz[0])
        iy1 = int(np.floor((istep)/(dimz[0])))
        state[ix1,iy1] = state[ix1,iy1] + incrementby  # Update positive pixels to a different color
    return state

def particleCount(plist,dimz):
    """Function that takes a list of positive pixels and the dimensions of the chemigram
    and returns a count of the particles."""
    state = updateState(plist,dimz)
    t = 1
    for istep in plist: # Iterate through the positive and gets their matrix coordinates
        ix1 = np.remainder((istep),dimz[0])
        iy1 = int(np.floor((istep)/(dimz[0])))
        visiter(ix1,iy1,state,t,dimz) # Recursively visits all connected 
        # particles and updates them to the same unique value
        t+=2

    particles = 0 
    grouplist = []
    
    sizelist = []
    negative = 2
    
    for istep in plist: # 
        ix1 = np.remainder((istep),dimz[0])
        iy1 = int(np.floor((istep)/(dimz[0])))
        sizelist.append(state[ix1,iy1]) # THIS IS FOR SIZING
        if state[ix1,iy1] == 2:
            particles+=1
            grouplist.append(negative)
            negative -=3
        elif state[ix1,iy1] != 0 and state[ix1,iy1] not in grouplist:
            particles+=1
            grouplist.append(state[ix1,iy1])


    return particles, grouplist,sizelist
    # particles is a count of the number of positive particles
    # grouplist is a list of the unique values that each group has. Its length is the number of particles. Has no repeat values.
    # sizelist is a list of the values of all positive pixels. Unless all positive pixels are not attached to one another, 
    # this list will be larger than grouplist. For example, if there is a particle that is made up of 10 pixels, the unique value
    # of those pixels will be added to sizelist ten times. grouplist and sizelist are passed to the "sizer" function later on.


def visiter(x,y,current,num,dimz): 
    """Function that updates the value of all matrix coordinates of attached pixels to the same unique value. 
    This is a recursive function. All diagonal pixels are assumed to be attached and originate from the same particle.
    Function takes x&y matrix coordinates, current matrix state, the incrementing value, and the dimensions of the matrix
    as arguments."""
    if x != (dimz[0]-1) and current[x+1,y] == 2:
        current[x+1,y] = num
        visiter(x+1,y,current,num,dimz)

    if x != (dimz[0]-1) and y != (dimz[1]-1) and current[x+1,y+1] == 2:
        current[x+1,y+1] = num
        visiter(x+1,y+1,current,num,dimz)

    if x != 0 and current[x-1,y] == 2:
        current[x-1,y] = num
        visiter(x-1,y,current,num,dimz)

    if x != 0 and y!= (dimz[1]-1) and current[x-1,y+1] == 2:
        current[x-1,y+1] = num
        visiter(x-1,y+1,current,num,dimz)

    if y != 0 and current[x,y-1] == 2:
        current[x,y-1] = num
        visiter(x,y-1,current,num,dimz)

    if x != 0 and y!= 0 and current[x-1,y-1] == 2:
        current[x-1,y-1] = num
        visiter(x-1,y-1,current,num,dimz)

    if y != (dimz[1]-1) and current[x,y+1] == 2:
        current[x,y+1] = num
        visiter(x,y+1,current,num,dimz)

    if x != (dimz[0]-1) and y!= 0 and current[x+1,y-1] == 2:
        current[x+1,y-1] = num
        visiter(x+1,y-1,current,num,dimz)
    else: 
        pass


def summary(hlist,likly,parts,truth):
    """Function that takes a list of researcher identified pixels and computer identified pixels and calculates
    which pixels are the same. This function mainly serves as a text summary of the computational analysis. Also 
    takes the number of particles as calculated by the particleCount function. If only a computational analysis, 
    no comparison to researcher data is displayed."""
    numright = 0
    for i in hlist:
        if i in likly:
            numright+=1

    if bool(truth) == False:
        print(f'Computer Found Particles: {parts}')
    else:
        print(f'User Found Particles    : {len(hlist)}')
        print(f'Computer Found Particles: {parts}')
        print(f'Intersection Number     : {numright}')
    


def transformer(hlist,likly,startx,starty,xdim,stepsize):  
    """This function takes a list of xy pixel coordinates determined by the researcher and a list of xy pixel
    coordinates determined by the computer and converts them into xy matrix coordinates. The starting x and y 
    coordinates of the field of view are require arguments as well as the x dimension of the chemigram."""
    human_pixels = []
    for i,j in hlist: 
        x = int(((i-startx)/stepsize))
        y = int(((j-starty)/stepsize)*xdim)
        human_pixels.append(x+y)
        
    comp_pixels = []
    for i,j in likly:  
        x = int(((i-startx)/stepsize))
        y = int(((j-starty)/stepsize)*xdim)    
        comp_pixels.append(x+y)  
        
    return human_pixels, comp_pixels


def sizer(grouplist,sizelist,stepsize):
    """Takes the sizelist and grouplist from particleCount function and determines the average diamter of each
    particle assuming that each particle is circular."""
    sizes=[]
    for i in grouplist:
        sizes.append(sizelist.count(i))
    return (np.sqrt((np.mean(sizes)*(stepsize*stepsize))/np.pi))*2


def printout(hpixels,cpixels,dimz): 
    """Function that creates a graphic comparing the manually identified particles and the computationally
    identified particles. Takes a list of pixel coordinates from the researcher and computer, as created by
    the "transformer" function"""
    state1 = updateState(hpixels,dimz,incrementby = 1)
    state2 = updateState(cpixels,dimz,update = True,statein = state1,incrementby = 2)

    img = state2
    plt.figure(figsize=(7,7))
    rotated_img = ndimage.rotate(img, 90)
    plt.imshow(rotated_img)
    plt.axis('off')
    plt.show()


def showChemigram(plist,dimz):
    """Function that creates a chemigram graphic showing the comptuationally identified particles"""
    img = updateState(plist,dimz)
    plt.figure(figsize=(7,7))
    rotated_img = ndimage.rotate(img, 90)
    plt.imshow(rotated_img)
    plt.axis('off')
    plt.show()



def showScreenShot(name):
    """If there is a screenshot of the filter surface in the given directory, that image is shown. 
    Takes directory name."""
    for filename in os.listdir(f"/Users/johnfox/Desktop/Microplastics/Code/{name}/"):
        if filename.endswith(".JPG"):
            newname = f"/Users/johnfox/Desktop/Microplastics/Code/{name}/" + filename
            im = Image.open(newname)
            width, height = im.size
              
            # Setting the points for cropped image
            # left = 1150
            # top = 0
            # right = width-100
            # bottom = height/2*1.2

            # im1 = im.crop((left, top, right, bottom))

            im1 = im
            plt.imshow(im1)
            plt.show()
            # display(im1)
            return
    print("No screenshot available in current directory")

def createHist(plist,dlist,coff):
    """Function that takes a list of the pixels that are positively identified as plastic by the computer,
    a list of reflection differences, and a cutoff. Prints out a histogram of the differences in reflection 
    between 2750 and 2850 wavenumbers. The cutoff is also visually shown as a vertical blueline. A red normal
    distribution is overlaid so that the assumption of a normal distribution can be qualitatively checked."""
    bins = np.linspace(start = min(plist), stop = max(plist), num=int(abs(max(plist) - min(plist)))*4)
    # sns.boxplot(x=dlist)
    plt.figure(figsize=(9,5))
    plt.hist(dlist,np.arange(-10,50,1))
    plt.vlines(coff,ymin = 0, ymax = 210,color = "blue",linestyle = "dotted")
    plt.ylim(0,210)
    plt.xlabel("Difference",size = 20)
    plt.ylabel("Occurrences",size = 20)
    plt.plot(bins,len(plist)*norm.pdf(bins,loc=np.mean(plist), scale=np.std(plist)),color = "red",
             linestyle = "dashed")
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.show()
    # sns.boxplot(x=plist)
    
def toGauss(ls):
    """Takes a list of reflection differences with outliers and returns a list of reflection without outliers.
    This function gets rid of outliers and then recursively calls itself until no outliers remain and then returns
    the list that is now "pared" or "pruned" of all outliers"""
    cutoff = mquantiles(ls)[2] + 1.5*stats.iqr(ls)
    newnumlist = []
    for i in ls:
        if i < cutoff:
            newnumlist.append(i)
    if len(newnumlist) == len(ls):
        return ls
    else:
        return toGauss(newnumlist)


def noGauss(ls):
    """"Similar to the toGauss function but only the first set of outliers are removed. No recursive calls to itself 
    and therefore fewer particles are removed as outliers which means that fewer particles are labelled as microplastics"""
    cutoff = mquantiles(ls)[2] + 1.5*stats.iqr(ls)
    newnumlist = []
    for i in ls:
        if i < cutoff:
            newnumlist.append(i)
    return newnumlist











