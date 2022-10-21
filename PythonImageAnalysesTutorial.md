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


## Numpy: Do math on a list = play with color (=do math) within an image (=list) ##
Lecture 11; how multiply value within a list (as `a*3` output concatened list 3 times)
https://www.youtube.com/watch?v=4uFs1qouPEI&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr&index=12
```python
XXX
```



```

# Usefull command: #
- type(VARIABLE) = say which type is it (integer, string,..)









