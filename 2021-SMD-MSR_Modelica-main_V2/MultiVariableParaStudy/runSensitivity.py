# -*- coding: utf-8 -*-


# import packages
import re
import numpy as np
import os
import subprocess
import random
from shutil import copyfile
from scipy import interpolate

simBasePath = os.getcwd() 

numberOfRuns = 1000
rndRangeLow = 30
rndRangeHigh = 200

run = []
CpCoefFuelVal = []
HTCcoefCoreVal = []
HTCcoefPHXVal = []

for runNumber in range(1, numberOfRuns+1):

    CpCoefFuel = (random.randint(rndRangeLow , rndRangeHigh))/100
    HTCcoefCore = (random.randint(rndRangeLow , rndRangeHigh))/100
    HTCcoefPHX = (random.randint(rndRangeLow , rndRangeHigh))/100
    
    CpFuel = 'parameter Real CpCoefFuel = ' + str(CpCoefFuel)
    HTCcore = 'parameter Real HTCcoefCore = ' + str(HTCcoefCore)
    HTCphx = 'parameter Real HTCcoefPHX = ' + str(HTCcoefPHX)
    
    #Make dir 
    workPath = "../run"+'{:.0f}'.format(runNumber) 
    os.mkdir(workPath)
    
    reading_file = open("MSRE.mo", "r")

    modelFileText = ""
    
    for line in reading_file:
        searchLine = line.strip()
        newLine1 = re.sub('parameter Real CpCoefFuel = 1',CpFuel,searchLine)
        newLine2 = re.sub('parameter Real HTCcoefCore = 1',HTCcore,newLine1)
        newLine3 = re.sub('parameter Real HTCcoefPHX = 1',HTCphx,newLine2)
            
        modelFileText = modelFileText + newLine3 +"\n"
        
    reading_file.close()
    
    modelFileName = '{}'.format(workPath)+"/MSRE.mo"
    modelFile  = open(os.path.join(workPath, modelFileName), mode='w+')
    modelFile.write(modelFileText)
    modelFile.close()
    
    #Write depl script
    runSim = '//SMD-MSR_Modelica \n\
    loadFile("SMD_MSR_Modelica.mo"); \n\
    loadFile("MSRE.mo"); \n\
    simulate(MSRE,startTime=0,stopTime=10000,numberOfIntervals=10000,tolerance=1E-6,method=dassl,outputFormat = "csv",fileNamePrefix ="run'+'{:.0f}'.format(runNumber)+'"); \n'
    
    #Save depl script in work dir
    runSimFileName = "{}/runMSRE.mos".format(workPath)
    deplFile = open(os.path.join(workPath, runSimFileName), mode='w+')
    deplFile.write(runSim)
    deplFile.close()
    
    #Copy SMD_Modelica Library 
    modelicaLibrary = "{}/SMD_MSR_Modelica.mo".format(workPath)
    copyfile("SMD_MSR_Modelica.mo", modelicaLibrary) 
    
    
    #Give permission to create an executable
    bashCommand1 = "chmod +x {}".format(workPath) 
    os.system(bashCommand1)
    
    #Execute omc with bash time as a python subprocess at workPath  
    simOutput,timeOutput = subprocess.Popen(['time', '-p', 'omc', 'runMSRE.mos'],cwd=workPath,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
    run.append(runNumber)
    CpCoefFuelVal.append(CpCoefFuel)
    HTCcoefCoreVal.append(HTCcoefCore)
    HTCcoefPHXVal.append(HTCcoefPHX)

#Write timeResults as a .m file in proper format to plot using matlab
with open("senResultsX.m", "w") as text_file:
    print("Run = {}".format(run)+";", file=text_file)
    print("CpCoefFuel = {}".format(CpCoefFuelVal)+";", file=text_file)
    print("HTCcoefCore = {}".format(HTCcoefCoreVal)+";", file=text_file)
    print("HTCcoefPHX = {}".format(HTCcoefPHXVal)+";", file=text_file)
    
