# -*- coding: utf-8 -*-
"""
Project: SMD-MSR_Modelica
Author: Visura Pathirana
Advisor: Dr. Ondrej Chvala 

The script performs following,
1. Read depletion data from a text file
2. Interpolate depletion data 
3. create a folder for depletion point desired 
4. open MSRE.mo modelica model file in simBase and write MSRE_deplPoint.mo with neutronics data
5. copy model files from SimBase 
6. write a modelica .mos file to load model and SMD_MSR library and to save results with deplPont.csv  
7. Execute the model by running .mos script with omc command 
"""

# import packages
import re
import numpy as np
import os
from shutil import copyfile
from scipy import interpolate

simBasePath = os.getcwd() 

#Import depletion data from .txt file 
#First line of the text file is coloumn headings and coloumns are divided by tab delimeter
deplTime,beta1,beta2,beta3,beta4,beta5,beta6,ngt,alphaFuel,alphaGrap = np.loadtxt('kin_dyn_edit.txt', skiprows=1, delimiter='\t',unpack=True)

#Interpolate data 
beta1Intpl = interpolate.InterpolatedUnivariateSpline(deplTime,beta1)
beta2Intpl = interpolate.InterpolatedUnivariateSpline(deplTime,beta2)
beta3Intpl = interpolate.InterpolatedUnivariateSpline(deplTime,beta3)
beta4Intpl = interpolate.InterpolatedUnivariateSpline(deplTime,beta4)
beta5Intpl = interpolate.InterpolatedUnivariateSpline(deplTime,beta5)
beta6Intpl = interpolate.InterpolatedUnivariateSpline(deplTime,beta6)
ngtIntpl = interpolate.InterpolatedUnivariateSpline(deplTime,ngt)
alphaFuelIntpl = interpolate.InterpolatedUnivariateSpline(deplTime,alphaFuel)
alphaGrapIntpl = interpolate.InterpolatedUnivariateSpline(deplTime,alphaGrap)

#Depletion range
depletion_range = np.linspace(deplTime[0],deplTime[len(deplTime)-1],num=1)

for depletionPoint in depletion_range:
    
    workPath = "../depl"+'{:.0f}'.format(depletionPoint) 
    os.mkdir(workPath)

    betaG1 = beta1Intpl(depletionPoint)
    betaG2 = beta2Intpl(depletionPoint)
    betaG3 = beta3Intpl(depletionPoint)
    betaG4 = beta4Intpl(depletionPoint)
    betaG5 = beta5Intpl(depletionPoint)
    betaG6 = beta6Intpl(depletionPoint)
    betaG = [betaG1,betaG2,betaG3,betaG4,betaG5,betaG6]
    beta =   'beta = {' +  ','.join(map(str, betaG)) +'}'
    LAM = 'Lam = ' + str(ngtIntpl(depletionPoint))
    a_F = 'a_F = ' + str(alphaFuelIntpl(depletionPoint))
    a_G = 'a_G = ' + str(alphaGrapIntpl(depletionPoint))
    
    omega = 'omega = ' + str(1)
    sin_mag = 'sin_mag = ' + str(0)
    
    #Write depl script
    reading_file = open("MSRE.mo", "r")

    modelFileText = ""
    
    for line in reading_file:
        searchLine = line.strip()
        newLine1 = re.sub('beta = {1,2,3,4,5,6}',beta,searchLine)
        newLine2 = re.sub('Lam = 1',LAM,newLine1)
        newLine3 = re.sub('a_F = 1',a_F,newLine2)
        newLine4 = re.sub('a_G = 1',a_G,newLine3)
        newLine5 = re.sub('omega = 1',omega,newLine4)
        newLine6 = re.sub('sin_mag = 1',sin_mag,newLine5)
            
        
        modelFileText = modelFileText + newLine6 +"\n"
        
    reading_file.close()
   
    modelFileName = '{}'.format(workPath)+"/MSRE_depl"+'{:.0f}'.format(depletionPoint)+".mo"
    modelFile  = open(os.path.join(workPath, modelFileName), mode='w+')
    modelFile.write(modelFileText)
    modelFile.close()
    
    runDeplSim = '//SMD-MSR_Modelica \n\
    loadFile("SMD_MSR_Modelica.mo"); \n\
    loadFile("MSRE_depl' + '{:.0f}'.format(depletionPoint)+'.mo"); \n\
    simulate(MSRE,startTime=0,stopTime=10000,numberOfIntervals=100000,tolerance=1E-6,method=dassl,outputFormat = "csv",fileNamePrefix ="depl'+'{:.0f}'.format(depletionPoint)+'"); \n'
    
    #Save depl script in work dir
    runDeplName = "{}/runModelica.mos".format(workPath)
    runDepl = open(os.path.join(workPath, runDeplName), mode='w+')
    runDepl.write(runDeplSim)
    runDepl.close()
    
    #Copy SMD_Modelica Library 
    modelicaLibrary = "{}/SMD_MSR_Modelica.mo".format(workPath)
    copyfile("SMD_MSR_Modelica.mo", modelicaLibrary) 
        
    bashCommand1 = "chmod +x {}".format(workPath) 
    os.system(bashCommand1)
    
    bashCommand2 = "cd {} && omc runModelica.mos".format(workPath)
    os.system(bashCommand2)

    bashCommand3 = "cd {}".format(simBasePath)
    os.system(bashCommand3)