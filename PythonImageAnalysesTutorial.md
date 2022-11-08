# Tutorial for image analyses using Python #
Tutorial from Dr. Sreenivas Bhattiprolu: [github](https://github.com/bnsreenu/python_for_microscopists); [Youtube](https://www.youtube.com/watch?v=_xIcIrdIqlA&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr)\

Work in specific conda environment
```bash
conda create --name Image
conda activate Image
python
```
Use `python COMAND.py` to see output directly. Use *CTRL+Z* to exit\
to save image you can write clic on pic and save in jpg



## Bases ##
Video04-XXX\

Use python script to mark border/froniter of image object
```python
import cv2
from skimage.filters import sobel #script to detect frontier

img = cv2.imread("images/test.jpg", 0) # Import image, 0 indicate grey level ie. 2D; 1 is color
img2 = sobel(img)

cv2.imshow("pic", img) # Show image
cv2.imshow("edge", img2) # Show image
print(img.shape) # print the shape of the image

cv2.waitKey(0)
cv2.destroyAllWindows()
```
To run python script, tape `python temp_demo.py`\
**troubleshooting solution:** Module missing:
- `No module named cv2`; so install it `conda install opencv`
- `No module named 'skimage'`; so install it `conda install scikit-image`


Basic operation, directly in `python`, **`math` library can be used for basic operation**
```python
a = 5
b = 3

if a <= 5:
  print("Yeah thats right bro")
else:
  print("no dude!") # Print "Yeah thats right as a is inferior or less than 5
 ```
 Test a basic python script math operator, 
```pyton
# Math package
import math

mass = input("Please enter the mass:") # in kg
c = 10000000 # in m per s
a = 9

Force = mass * a
print("Force=", Force)
```
--> It ask for the input value and output the Force!\
## Test *if else elif statements* ##
```python
## Raw
num = int(input("Please enter an integer: ")) # convert string (default) into integer
if num >= 0:
  print("This number is positive")
else:
  print("This number is negative")

## More elaborate
num = int(input("Please enter an integer: ")) # convert string (default) into integer
if num >= 0:
  print(num, "is a positive number")
else:
  print(num, "is a negative number")
  
## Even more
num = int(input("Please enter an integer: ")) # convert string (default) into integer
if num >= 0:
  print("%d is a positive number" %(num))  # %d is used to define number, "d" indicate that is an integer; "f" can be used for float(=continuous) %f4.2f = show only 2 decimals
else:
  print("%d is a negative number" %(num)) 

## More more
num = int(input("Please enter an integer: ")) # convert string (default) into integer
if num >= 0:  # If over zero print
  print("%d is a positive number" %(num))
elif num==0: # AND if == 0 print
  print("%d is zero" %(num))
else:
  print("%d is a negative number" %(num))
  
## More more more
if num >= 0: # If that is true, go inside and execute line, if not, go to the else at the SAME LEVEL
  if num==0:
    print("Zero")
  else:
    print("%d is a positive number" %(num))
else: # else here refer to the first if (position~column is important!!)
  print("%d is a negative number" %(num))
```
## lists tuples and dictionaries ##
list is a list of number (image is also a list of number)
```python
## List

a = [1,2,3,4] # this a list
a = list(range(5))
type(a) # output that is a list
a[0] # output the first entry of the list (a[-1] = 5)
a.pop() # Remove the last entry of the a list
print(a[2:4])) # Here INCLUDE:EXCLUDE; thus print(a[2:3])) will show only 1 number

## Tuples (is a list unmutable, cannot be changed)
 a=(1,2,3)
 type(a) # output that is a tupple
 b = list(a) # convert tupple to list

## Dictionnary
fruit = {'banana':'yellow', 'apple':'red', 'pear':'green'}
print(fruit['apple']) # print output value, so red
fruit['grape'] = 'black' # Add grape to our dictionnary
del fruit['apple'] # delete apple from our dictionnary
```


## Numpy: Do math on a list = "play with color (=do math) within an image (=list)" ##
[Lecture 11](https://www.youtube.com/watch?v=4uFs1qouPEI&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr&index=12)\
how multiply value within a list (as `a*3` output our list concatenated 3 times)
**Using loop:**
```python
a = [1,3,5,1,4] # create a list
for i in a: # for each entry within the a list
  print(i*2) # multiply per 2
```
But result is not a list... To have list as output:
```python
a = [1,3,5,1,4] # create a list
square = [] # create an empty list
for i in a: # for each entry within the a list
  square.append(i*2) # multiply per 2 the "square" variable and populate it with i*2
  
print(square)
```
Simplify **loop** to **comprehenssion**
```python
square = [i**2 for i in a]
print(square)
# this also give the same:
print([i**2 for i in a])

## Now if we want to give the square of only even number using ifelse
square = [i**2 for i in a if i%2 == 0] # if you divide per two and the remainder is zero == "if i is an even number"
print(square) # here only show the even number
```
Using **numpy**
```python
import numpy as np # import numpy library and call it using np

a = [1,2,3,4,5,6]
b = np.array(a) # convert the a list into a numpy array

print(b * 2) # will multiply every element in a by 2 ! And not concatenated the list two times
```
Here is an example to create a matrix with a desired number; could be use to create a digital filter for an image (filter out pixel value over 5):
```python
a = np.full((3,3), 5) # create a 3 column 3 rows matrix with value of 5
```
Using **Slicing** to subset an array fron an array = part of an image/pixel
```python
a = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]]) # generate a matrix 3 * 4 arrays (4 columns 3 rows)
b = a[:2] # keep only the first 2 rows; same as --> a[0:2] 
b = a[:, :2] # ":" print all rows, ":2" from only the two columns
b = a[:,2] # print all rows but only the third columns
np.array([a[0,0], a[1,1]]) # print number from 1st column first row and second column second row

```
--> slice need to be read as `a[ROW,COLUMN]`\
To do sum of two matrices:
```python
x = np.array([[1,2], [3,4]], dtype=np.float64) # adding data type is optional, here float64 refer to a numeric character with decimal
y = np.array([[5,6], [7,8]], dtype=np.float64)
print(x+y) # print sum of the two matrices
print(np.add(x,y)) # same
np.sum(x, axis = 0) # comput the sum of each columns; axis = 1 for each rows
print(x.T) # transpose (column to rows and rows to column); usefull to invert color!
```
To do **broadcasting** = adding an array to a multidimensional array
```python
a = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12]]) 
b = np.array([[1,1,1,1]])
print(a+b) # will add one to each rows of the a multidimensional array
```
Now let's **read/modify an image with numpy**\
```python
## import libraries
import numpy as np
from imageio import imread, imsave

## import image
img = imread('images/Meristem_profile.jpg')

print(img) # print array of the image
print(img.dtype, img.shape) # print some image information; the last number is 3 as refer to red green blue

## image manipulation
img_tinted = img * [0., 1., 0.] # keep only the green; multiply per zero the other color
imsave('images/tinted_Meristem_profile.jpg', img_tinted)
```
--> We can play with color intensity of an image (keep only green; decrease other color...)

## for and while Loops ##

```python
while (count<10): 
  print("The count is: ", count) # So here it will print that the count is "0"
  count = count + 1 # here add +1 at the count
```
Conversion C to F using **while**
```python
# F = 9/5C + 32

C = -40
while C<= 40:
  F = ((9./5)*C) + 32
  print(C,F)
  C = C + 5
```
play with **break** in **loop** or **while**; to stop a loop or a while.
```
for letter in 'I love python':
  if letter == 'p':
    break # break the code if the letter is p
  print("Current letter is: ", letter)
```
```
a = 10
while a>0:
  print("Current value is: ", a)
  a = a-1
  if a==5:
    break
```
```
for num in range(10,20): # range start at 10 and exclude 20
  if num%2 ==0: # when you divide your number per two, if your reminder is 0 = even number
    print("%d is an even number" %(num))
  else:
    print("%d is an odd number" %(num))
```

## Python functions ##
write function for specific task and run it (e.g. to rotate; save; invert color...)
```python
## Create some functions
def my_function(): # here define function name as "my_function()"
  print("Hello I am a function")
def my_function(name): # here define function name as "my_function()"
  print("Hello your name is: ", name)
## Launch the function
my_function() # call the function
my_function("Thomas") # call the function here with an argument

## Add a default argument (parameter) to a function
def my_country(country = "France"):
  print("My country is: ", country)

## Function with a loop
### Function that print all parameters within a list
def my_function(my_list): 
  for x in my_list:
    print(x) 
    
fruits = ['apple', 'pear', 'banana']

my_function(fruits)

## Function using return
def my_function(x):
  return 5*x
  
my_function(3) # output 15
```
Create a function to change F into C
```python
# Simple version
def F(C):
  F_value = ((9.0/5)*C) + 32
  print('%d degrees Centigrade is equal to %d Fahreneit' %(C, F_value)) # Add place holder, first one is centigrade then farneeight

F(-40) # output C to F

# Version with a list
Centi = [-40, -20, -10, 0, 10, 20, 30, 40] # Create the list

def F(C):
  for i in C: # C is the list; i is the value in the loop
  
    F_value = ((9.0/5)*i) + 32
    print('%d degrees Centigrade is equal to %d Fahreneit' %(i, F_value))

F(Centi) # output math for the list
```
Create a function to modify image
```python
## import libraries
from scipy import misc, ndimage # For basic image modifications
import imageio # To save and read image

def rotated(img): # take an image as an argument
  rotated_img = ndimage.rotate(img, 45) # take img and rotate by 45 deg
  imageio.imsave('images/rotated.jpg', rotated_img) # save the img

def blurred(img): 
  blurred_img = ndimage.gaussian_filter(img, sigma=5) # take img and blur it
  imageio.imsave('images/blurred.jpg', blurred_img) # save the img

def denoised(img): 
  denoised_img = ndimage.median_filter(img, 3) # take img and blur it
  imageio.imsave('images/denoised.jpg', denoised_img) # save the img
  
img1 = imageio.imread('images/Meristem_profile.jpg')

rotated(img1)
blurred(img1)
denoised(img1)
```
## Python classes ##
*attributes* = parameter into the function
Count each cells and define the distance from other cell

**class** work almost like a **function**

```python
class CellDemo: # Define the class "Cell"
  pass
  
cell_1 = CellDemo() # instance 1 (= object 1)
cell_2 = CellDemo() # instance 2

print(cell_1)
print(cell_2)

# specify position to cell
cell_1.name = "cell_one"
cell_1.x = 0
cell_1.y = 0

cell_2.name = "cell_two"
cell_2.x = 10
cell_2.y = 5
```
Now let's uses classes to simplify the code
```python
from math import sqrt

# Define the class "Cell"
class Cell: 
  def __init__(self, name, x, y): # def for defining a function; initialize a function with the parameters name x and y
    self.name = name   # almost equal to: cell_1.name = "cell_one"
    self.x = x
    self.y = y
    
# function to calculate distance within 2 cells
  def cell_distance(self, other_cell):
    distance = sqrt((self.x - other_cell.x)**2 + (self.y - other_cell.y)**2) # basic math for distance between coordinates
    return distance 
 
# now create instances (= objects)
cell_1 = Cell('cell_one', 0, 0)
cell_2 = Cell('cell_two', 10, 5)
cell_3 = Cell('cell_two', 20, 9)


#
distance = cell_1.cell_distance(cell_2) # cell1 vs distance to cell_2 which is called using the cell_distance() function changing 'other_cell' to 'cell_2'
print("The cells are %f units apart." % distance)
```
--> Instances could be defined automatically if we used cell/object detection before

## scikit-image ##
Edge detection using skimage.filters
```python
from skimage import io
from matplotlib import pyplot as plt # to install: conda install -c conda-forge matplotlib
import imageio # To save and read image

img = imageio.imread('images/protoplast_1.jpg', as_gray = True) # read image and scale to grey

from skimage.filters import roberts, sobel, scharr, prewitt

egde_roberts = roberts(img)
edge_sobel = sobel(img)
edge_scharr = scharr(img)
edge_prewitt = prewitt(img)

fig, axes = plt.subplots(nrows = 2, ncols = 2, sharex = True, sharey = True, figsize = (8,8)) # used to make multi plots

ax = axes.ravel()

ax[0].imshow(img, cmap=plt.cm.gray)
ax[0].set_title('Original image')

ax[1].imshow(edge_sobel, cmap=plt.cm.gray)
ax[1].set_title('sobel image')

ax[2].imshow(edge_scharr, cmap=plt.cm.gray)
ax[2].set_title('scharr image')

ax[3].imshow(edge_prewitt, cmap=plt.cm.gray)
ax[3].set_title('prewitt image')

for a in ax:
  a.axis('off')
 
plt.tight_layout()
plt.show()
```
**Kenny** Edge detection using skimage.features = output is a binary image
```python
from skimage import io
from matplotlib import pyplot as plt # to install: conda install -c conda-forge matplotlib
import imageio # To save and read image

img = imageio.imread('images/protoplast_1.jpg', as_gray = True) # read image and scale to grey

from skimage.feature import canny # canny do noise detection, gradient detection...
edge_canny = canny(img, sigma=3)
plt.imshow(edge_canny)
plt.show()
```
Lets calculate the speed by which cell colonize an area ("scratch assay")
Import scratch images from tutorial [repo](https://github.com/bnsreenu/python_for_microscopists/tree/master/images/scratch_assay)
```python
# Import packages
import matplotlib.pyplot as plt
from skimage import io, restoration
from skimage.filters.rank import entropy # entropy detect the mess; not clean area; in an image, 
from skimage.morphology import disk
import numpy as np

# Import image
img = io.imread("images/scratch_assay/Scratch0.jpg")

# Apply entropy filter
entr_img = entropy(img, disk(3)) # disk is the size of the disk that detect clean area
plt.imshow(entr_img, cmap = 'gray')
plt.show() # We see the scratch (clean area) is black and side with cells is white

# Add treshold on entropy image to detect the middle/scratch/clean area

from skimage.filters import try_all_threshold
fig, ax = try_all_threshold(entr_img, figsize=(10,8), verbose=False) # Calculate some different tresholdfrom the image
plt.show() # Otsu looks better

# Apply the treshold to the entropy image

from skimage.filters import threshold_otsu
thresh = threshold_otsu(entr_img) # Give the threshold value of Osu filter

binary = entr_img <= thresh # Modify binary/entropy image as all dark pixel are False  and white are True
plt.imshow(binary, cmap = 'gray')
plt.show() # We have a clean black and white image

# Then calculate the number of white pixel = area that needs to be colonize
print("The percent bright pixels is: ", (np.sum(binary==1)/((np.sum(binary==1))+(np.sum(binary==0))))) # Here we count the number of True / total values = Percent of bright pixels
```
Simplify the scratch assay code with loop:
```python
# Import packages
import matplotlib.pyplot as plt
from skimage import io, restoration
from skimage.filters.rank import entropy # entropy detect the mess; not clean area; in an image, 
from skimage.morphology import disk
import numpy as np
from skimage.filters import try_all_threshold
from skimage.filters import threshold_otsu
import glob
import os # To make sure files are in the correct order

# Prerequistes
time = 0 # Image 0 = time is 0, image 1; time is 1, etc...
time_list = [] # Create an empty list to be populated
area_list = [] # Create an empty list to be populated
path = "images/scratch_assay/*.*" # Read all files within this folder
list_of_files = sorted(filter(os.path.isfile, glob.glob(path)))

for file in list_of_files: 
  img = io.imread(file) # Import image
  entr_img = entropy(img, disk(10)) # Apply entropy
  thresh = threshold_otsu(entr_img) # Apply threshold
  binary = entr_img <= thresh 
  scratch_area = np.sum(binary == True)
#  print(time, scratch_area)
  time_list.append(time) 
  area_list.append(scratch_area)
  time+=1 # Time equal to time + 1

plt.plot(time_list, area_list, 'bo')
plt.show() # Show the time as function of pixels numbers

# To add linear regression

from scipy.stats import linregress
slope, intercept, r_value, p_value, std_err = linregress(time_list, area_list)
print("y= ", slope, "x", " + ", intercept)
print("R\N{SUPERSCRIPT TWO} = ", r_value**2)
```
Real examples with Shalini data:
```python
# Import packages
import matplotlib.pyplot as plt
from skimage import io, restoration
from skimage.filters.rank import entropy # entropy detect the mess; not clean area; in an image, 
from skimage.morphology import disk
import numpy as np
from skimage.filters import try_all_threshold
from skimage.filters import threshold_otsu
import glob
import os # To make sure files are in the correct order

# Modify images to keep only the green color
from imageio import imread, imsave

# Lets do it ugly way one by one but a loop could be performed

FRP_0_green = io.imread('images/Shalini/FRP_0.jpg') * [0., 1., 0.] # keep only the green; multiply per zero the other color
imsave('images/Shalini/FRP_0_green.jpg', FRP_0_green)
FRP_1_green = io.imread('images/Shalini/FRP_1.jpg') * [0., 1., 0.] # keep only the green; multiply per zero the other color
imsave('images/Shalini/FRP_1_green.jpg', FRP_1_green)
FRP_2_green = io.imread('images/Shalini/FRP_2.jpg') * [0., 1., 0.] # keep only the green; multiply per zero the other color
imsave('images/Shalini/FRP_2_green.jpg', FRP_2_green)
FRP_3_green = io.imread('images/Shalini/FRP_3.jpg') * [0., 1., 0.] # keep only the green; multiply per zero the other color
imsave('images/Shalini/FRP_3_green.jpg', FRP_3_green)
FRP_4_green = io.imread('images/Shalini/FRP_4.jpg') * [0., 1., 0.] # keep only the green; multiply per zero the other color
imsave('images/Shalini/FRP_4_green.jpg', FRP_4_green)
FRP_5_green = io.imread('images/Shalini/FRP_5.jpg') * [0., 1., 0.] # keep only the green; multiply per zero the other color
imsave('images/Shalini/FRP_5_green.jpg', FRP_5_green)
FRP_6_green = io.imread('images/Shalini/FRP_6.jpg') * [0., 1., 0.] # keep only the green; multiply per zero the other color
imsave('images/Shalini/FRP_6_green.jpg', FRP_6_green)


# Prerequistes
time = 0 # Image 0 = time is 0, image 1; time is 1, etc...
time_list = [] # Create an empty list to be populated
area_list = [] # Create an empty list to be populated
path = "images/Shalini_green/*.*" # Read all files within this folder
list_of_files = sorted(filter(os.path.isfile, glob.glob(path)))

#############################################
##### Find the best entropy parametrs #####
#############################################
img = io.imread("images/Shalini_green/FRP_1_green.jpg", as_gray = True)


# Apply entropy filter
entr_img = entropy(img, disk(1)) # disk is the size of the disk that detect clean area
plt.imshow(entr_img, cmap = 'gray')
plt.show() # We see the scratch (clean area) is black and side with cells is white

# Add treshold on entropy image to detect the middle/scratch/clean area

from skimage.filters import try_all_threshold
fig, ax = try_all_threshold(entr_img, figsize=(10,8), verbose=False) # Calculate some different tresholdfrom the image
plt.show() # Otsu looks better

#############################################

# Launch the program on all pictures

for file in list_of_files: 
  img = io.imread(file, as_gray = True) # Import image
  entr_img = entropy(img, disk(1)) # Apply entropy
  thresh = threshold_otsu(entr_img) # Apply threshold
  binary = entr_img <= thresh 
  scratch_area = np.sum(binary == False)
#  print(time, scratch_area)
  time_list.append(time) 
  area_list.append(scratch_area)
  time+=1 # Time equal to time + 1

plt.plot(time_list, area_list, 'bo')
plt.show() # Show the time as function of pixels numbers
```

 ## Denoising microscope images ##
Gaussian kurnel, the sum of all pixel within the image is 1. 
```python
from skimage import io # import image
from scipy import ndimage as nd # denoise image
from matplotlib import pyplot as plt # make subplots
from skimage.restoration import denoise_nl_means, estimate_sigma

# Gaussian denoising
img = io.imread("images/denoised.jpg")
guassian_img = nd.gaussian_filter(img, sigma = 3)

plt.imshow(guassian_img)
plt.show() 

# Median denoising
img = io.imread("images/denoised.jpg")
median_img = nd.median_filter(img, size = 3)

plt.imshow(median_img)
plt.show() 


# nlm denoising
sigma_est = np.mean(estimate_sigma(img, multichannel=True))
patch_kw = dict(patch_size=5, patch_distance=3, multichannel=True)

nlm = denoise_nl_means(img, h= 1.15 * sigma_est, fast_mode=True, **patch_kw) # Here **patch_kw  mean we are unpacking a dictionnary called "patch"; it is the same as taping what is inside dict

# plt.imsave("images/denoised_gaussian.jpg", guassian_img)
```
--> Can be used to denoise noisy, blurry images with not clear edges

## Histogram based image segmentation ##
```python
from skimage.restoration import denoise_nl_means, estimate_sigma
from skimage import data, img_as_float, img_as_ubyte # convert images to float and ubyte
import numpy as np
from matplotlib import pyplot as plt

# import img as float
img = img_as_float(io.imread("images/protoplast_1.jpg"))

# denoise img
patch_kw = dict(patch_size=5, patch_distance=3, multichannel=True) # Create the patch_kw dictionary
denoise = denoise_nl_means(img, h= 1.15 * sigma_est, fast_mode=True, **patch_kw) 
plt.imshow(denoise)
plt.show() 

##  --> Denoising is not needed for our example here!

# make histogram segmentation
img_ubyte = img_as_ubyte(img) # convert img to ubyte = uint8 (from 0 to 255 range)
plt.hist(img_ubyte.flat, bins = 100, range=(0,255)) # 
plt.show() # Here lets manually segment our picture by individuals peaks

segm1 = (img_ubyte <= 50)
segm2 = (img_ubyte > 50) & (img_ubyte <= 100)
segm3 = (img_ubyte > 100)

all_segments = np.zeros((img_ubyte.shape[0], img_ubyte.shape[1], 3)) # Create a blamk image as the same size as our image; the 3 correspond to the three color channel

all_segments[segm1] = (1,0,0) # all segment 1 pixel color in red
all_segments[segm2] = (0,1,0)
all_segments[segm3] = (0,0,1)

plt.imshow(all_segments)
plt.show()

# clean the image --> That help removing tiny patch that spread in non-desired region
from scipy import ndimage as nd

segm1_opened = nd.binary_opening(segm1, np.ones((3,3)))
segm1_closed = nd.binary_closing(segm1_opened, np.ones((3,3)))

segm2_opened = nd.binary_opening(segm2, np.ones((3,3)))
segm2_closed = nd.binary_closing(segm2_opened, np.ones((3,3)))

segm3_opened = nd.binary_opening(segm3, np.ones((3,3)))
segm3_closed = nd.binary_closing(segm3_opened, np.ones((3,3)))

all_segments_cleaned = np.zeros((img_ubyte.shape[0], img_ubyte.shape[1], 3))

all_segments_cleaned[segm1_closed] = (1,0,0) # all segment 1 pixel color in red
all_segments_cleaned[segm2_closed] = (0,1,0)
all_segments_cleaned[segm3_closed] = (0,0,1)

plt.imshow(all_segments_cleaned)
plt.show()
```
Let's try histogram segmentation with a high quality microscopy image
```python
XXX Or follow tutorial as there is probably better segmentation method
```

## Random Walker segmentation ##











CHUI AL :  https://www.youtube.com/watch?v=jcUx-TQpcM8&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr&index=22



## Tutorial I looked quickly: ##
17 - Reading images in Python: https://www.youtube.com/watch?v=52pMFnkDU-4&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr&index=18
Various way to import image
18 - Image processing using pillow in Python https://www.youtube.com/watch?v=s_hDL2fGvow&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr&index=19
Pillow is a package for basic image processing: croping, resizing, changing color. Not for Machine Learning stuff such as segmentation, object recognition (opencv, sicket.image, ...)
19 - image processing using scipy in Python https://www.youtube.com/watch?v=s_hDL2fGvow&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr&index=20
scipy library part of numpy stack. Could be used to plot different images in one plot. To filter some pixel values from an image




https://www.youtube.com/watch?v=sLc0O-8RpW0&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr&index=15
# Usefull command: #
- srun --x11 --nodelist=node03 --mem=20g --pty bash -l
- type(VARIABLE) = say which type is it (integer, string,..)

```






