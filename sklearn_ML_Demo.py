# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 13:14:53 2018

@author: ste_j
"""
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
#from sklearn.ensemble import AdaBoostClassifier
from sklearn.svm import SVC
#from sklearn.neural_network import MLPClassifier

from sklearn.metrics import confusion_matrix, roc_curve, roc_auc_score, matthews_corrcoef, accuracy_score

from random import shuffle

import pandas
import numpy as np

data = pandas.read_csv('C:\Users\ste_j\Desktop\UCL-Test\Github\Cryo-Lead-ML\Individual_RF_mod.csv')

#data_shuffle = data.sample(frac=1).reset_index(drop=True)

data_x = data[['STD','Skew','kurtosis','PP','Sigma']]
data_y = data [['Label']]

X_train, X_test, Y_train, Y_test = train_test_split(data_x, data_y, test_size=0.33, random_state=42)

# Random Forests
rf_clf = RandomForestClassifier(n_estimators=1000)

rf_clf.fit(X_train, Y_train)

prediction = rf_clf.predict(X_test)
prediction_list = prediction.tolist()

cnf_matrix = confusion_matrix(Y_test, prediction)
accuracy = accuracy_score(Y_test, prediction)

# Support Vector Machines
svc_clf = SVC()
svc_clf.fit(X_train, Y_train)

prediction = svc_clf.predict(X_test)
prediction_list = prediction.tolist()

cnf_matrix = confusion_matrix(Y_test, prediction)
accuracy = accuracy_score(Y_test, prediction)
