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






