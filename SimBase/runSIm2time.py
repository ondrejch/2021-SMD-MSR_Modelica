# -*- coding: utf-8 -*-
"""
Project: SMD-MSR_Modelica
Author: Visura Pathirana
Advisor: Dr. Ondrej Chvala 

The script performs following for a given numberOfRuns
1. create a folder 
2. create .mos file to load model and simulate 
3. run the simulation with timing using bash time and omc 
4. extract times using regex and split str list to make a number list using map
5. print time results to a .m file for easy plotting using matlab 

"""

# import packages
import re
import numpy as np
import os
import subprocess
from shutil import copyfile
from scipy import interpolate

simBasePath = os.getcwd() 

numberOfRuns = 5

realTime = []
userTime = []
sysTime = []

for runNumber in range(1, numberOfRuns+1):

    #Make dir 
    workPath = "../run"+'{:.0f}'.format(runNumber) 
    os.mkdir(workPath)
    
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
    
    #Copy SMD_Modelica Model File
    modelFile = "{}/MSRE.mo".format(workPath)
    copyfile("MSRE.mo", modelFile)
    
    #Give permission to create an executable
    bashCommand1 = "chmod +x {}".format(workPath) 
    os.system(bashCommand1)
    
    #Execute omc with bash time as a python subprocess at workPath  
    simOutput,timeOutput = subprocess.Popen(['time', '-p', 'omc', 'runMSRE.mos'],cwd=workPath,stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    
    #Write simOutput to a text file for inspection
    simOutputFileName = "{}/simOutput.txt".format(workPath)
    simOutputFile = open(os.path.join(workPath, simOutputFileName), mode='w+')
    simOutputFile.write(simOutput.decode("utf-8"))
    simOutputFile.close()
    
    #Write timeOutput to a text file for inspcection
    timeOutputFileName = "{}/timeOutput.txt".format(workPath)
    timeOutputFile = open(os.path.join(workPath, timeOutputFileName), mode='w+')
    timeOutputFile.write(timeOutput.decode("utf-8"))
    timeOutputFile.close()
  
    #Extract time values by stripping str and create a list using append
    runTimeVal = re.findall('\d*\.\d+', timeOutput.decode("utf-8"), re.IGNORECASE)
    realTime.append(runTimeVal[0])
    userTime.append(runTimeVal[1])
    sysTime.append(runTimeVal[2])

#Strip list using map to make char list to numeric list    
realTime = list(map(float, realTime))
userTime = list(map(float, userTime))
sysTime = list(map(float, sysTime))

#Write timeResults as a .m file in proper format to plot using matlab
with open("timeResults.m", "w") as text_file:
    print("realTime = {}".format(realTime)+";", file=text_file)
    print("userTime = {}".format(userTime)+";", file=text_file)
    print("sysTime = {}".format(sysTime)+";", file=text_file)
