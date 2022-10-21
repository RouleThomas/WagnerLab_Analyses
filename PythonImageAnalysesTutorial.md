# Tutorial for image analyses using Python #
Tutorial from Dr. Sreenivas Bhattiprolu: [github](https://github.com/bnsreenu/python_for_microscopists); [Youtube](https://www.youtube.com/watch?v=_xIcIrdIqlA&list=PLZsOBAyNTZwYHBIlu_PUO19M7aHMgwBJr)\

Work in specific conda environment
```bash
conda create --name Image
conda activate Image
python
```
Use `python COMAND.py` to see output directly

## Bases ##
Video04-XXX\

Use python script to mark border/froniter of image object
```python
import cv2
from skimage.filters import sobel #script to detect frontier

img = cv2.imread("images/test.jpeg", 0) # Import image, 0 indicate black/white ie. 2D
img2 = sobel(img) # Transform img
print(img.shape)

cv2.waitKey(0)
cv2.destroyAllWindows()
```
To run python script, tape `python NAMESCRIPT.py`\
**troubleshooting:** Module missing:
- `No module named cv2`; so install it `conda install opencv`
- ` No module named 'skimage'`; so install it `conda install scikit-image`
- 

maybe modify script to save image
