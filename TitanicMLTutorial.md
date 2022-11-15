Tutorial to learn Machine Learning in Python can be found [here](https://www.kaggle.com/code/alexisbcook/titanic-tutorial/notebook)\
Let's use Image conda environment for the tutorial:
```bash
srun --x11 --nodelist=node03 --mem=20g --pty bash -l
conda activate Image
```
Download and import manually dataset to HPC and unzip with `unzip titanic.zip`\
Using the patterns you find in train.csv, you have to predict whether the other 418 passengers on board (in test.csv) survived ( 1 = survived).
**Read data path**
```python
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)

# Input data files are available in the "../input/" directory.

import os
for dirname, _, filenames in os.walk('input/'):
    for filename in filenames:
        print(os.path.join(dirname, filename))
```
--> list all files under the input directory\
**import data**
```python
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import os

train_data = pd.read_csv("input/train.csv")
test_data = pd.read_csv("input/test.csv")

train_data.head() # print head(5)

# Percent of woman that survive
women = train_data.loc[train_data.Sex == 'female']["Survived"] # keep only Sex and survive column and keep the woman
print("% of women who survived:", sum(women)/len(women))

# Percent of men that survive
men = train_data.loc[train_data.Sex == 'male']["Survived"]
rate_men = sum(men)/len(men)
print("% of men who survived:", rate_men)
```
More female survive, now let's try to take into account ALL columns into our data to better predict survival\
**Random Forest Model**: 
```python
from sklearn.ensemble import RandomForestClassifier

y = train_data["Survived"]

features = ["Pclass", "Sex", "SibSp", "Parch"] # feature to take into account
X = pd.get_dummies(train_data[features]) # converts categorical data, here sex into dummy= 0 or 1
X_test = pd.get_dummies(test_data[features])

model = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=1) # n_estimator =nb of trees
model.fit(X, y)
predictions = model.predict(X_test)

output = pd.DataFrame({'PassengerId': test_data.PassengerId, 'Survived': predictions})
output.to_csv('submission.csv', index=False)
print("Your submission was successfully saved!")
```
--> Create a file that predict whether passenger survived taking features into account...\
Another way tutorial from [here](https://keytodatascience.com/machine-learning-titanic-disaster-problem/)
older rich women and children were the most likely to survive and poor middle-aged men were the least likely to survive.
```python
from sklearn.model_selection import train_test_split

train_data.describe() # give summary statistics

train_x = train_data[["Pclass", "Sex"]] # select only class and sex, not age as there is NA value and Sklearn cannot deal with that, need permutation
train_y = train_data[["Survived"]] # this is the output target variable

# Clean data: Making Male/ Female to integer numbers
train_x["Sex"].replace("male", 1, inplace = True)
train_x["Sex"].replace("female", 0, inplace = True)
train_x.head()

# create test file from the train file (take 30% of the train file to use it as test file)
tr_x, cv_x, tr_y, cv_y   = train_test_split(train_x, train_y, test_size = 0.30)
```
--> r_x & tr_y are the training input and output and cv_x & cv_y are cross-validation input and output.\
Now machine learning with RandomForestClassifier:\
```python
# Call the Machine Learning Algorithm
rf = RandomForestClassifier()

# Fitting and training the above called algorithm
rf.fit(tr_x, tr_y) # fit() train the algorithm; take the input tr_x and learn the expected output tr_y

#check how well the model perform on cross-validation data
Accuracy_RandomForest = rf.score(cv_x, cv_y)
print("Accuracy = {}%".format(Accuracy_RandomForest * 100))
```
Now machine learning with regression linear
```python
from sklearn.linear_model import LogisticRegression


# Call the Machine Learning Algorithm
lgr = LogisticRegression()
lgr.fit(tr_x, tr_y)

# Fitting and training the above called algorithm
rf.fit(tr_x, tr_y) # fit() train the algorithm; take the input tr_x and learn the expected output tr_y

#check how well the model perform on cross-validation data
Accuracy_LogisticRegression = lgr.score(cv_x, cv_y)
print("Accuracy = {}%".format(Accuracy_LogisticRegression * 100))

###### now do prediction on the test file
test_x = test_data[["Pclass", "Sex"]] # select only class and sex, not age as there is NA value and Sklearn cannot deal with that, need permutation

# Clean data: Making Male/ Female to integer numbers
test_x["Sex"].replace("male", 1, inplace = True)
test_x["Sex"].replace("female", 0, inplace = True)
test_x.head()

# prediction
prd = rf.predict(test_x) # random forest prediction
op = test_data[['PassengerId']]

op['Survived'] = prd
op.head()

op.to_csv("Submission1.csv",index=False)

```






