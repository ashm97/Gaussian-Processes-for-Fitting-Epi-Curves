# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 08:41:59 2020

@author: Ash
"""


import numpy as np
from fitderiv import fitderiv
import matplotlib.pyplot as plt

# import data
d= np.genfromtxt('./Data/y_true_20_04_2020.csv', delimiter= ',')
d= np.delete(d, 0, 0) # remove header
t= d[:, 0]
od= d[:, 2] # take just mean y_true column








# Loop GP over expanding window 
for i in range(1,(len(t) + 1)):
    #print(i)
    t_prime = np.resize(t,i)
    od_prime = np.resize(od,i)
    
    plot_filename = "plot{0}.png".format(i)
    dat_filename = "fitDat{0}.csv".format(i)
    
    q_exp= fitderiv(t_prime, od_prime, bd= {0: [-10, 10],1: [-10, 10],2: [-10, 10]}, cvfn= 'sqexp', stats= True, esterrs= True)
    
    #Save Figure
    plt.figure()
    plt.subplot(2,1,1)
    q_exp.plotfit('f')
    plt.ylabel('Log(Y_true)')
    plt.subplot(2,1,2)
    #q_exp.plotfit('df')
    plt.plot(np.exp(q_exp.df))
    plt.xlabel('Time')
    plt.ylabel('Growth Rate')
    plt.savefig(plot_filename)
    
    #save Data
    q_exp.export(dat_filename)
    
    
    
    
    
    
    