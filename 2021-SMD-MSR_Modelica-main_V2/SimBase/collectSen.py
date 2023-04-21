# -*- coding: utf-8 -*-

# import packages
import re
import numpy as np
import numpy as np
import pandas as pd 
import os
import subprocess
import random
from shutil import copyfile
from scipy import interpolate

simBasePath = os.getcwd() 

numberOfRuns = 1000

run = []

maxPower = []
minPower = []
maxDerPower = []
maxDerFuelTemp1 = []
maxDerFuelTemp2 = []
maxDerGrapTemp = []
maxInletTemp = []
maxOutletTemp = []
	
SSPower = []
SSDerPower = []
SSDerFuelTemp1 = []
SSDerFuelTemp2 = []
SSDerGrapTemp = []
SSInletTemp = []
SSOutletTemp = []
	
for runNumber in range(1, numberOfRuns+1):

	workPath = "../run"+'{:.0f}'.format(runNumber) 
	dataFileName = '/run'+'{:.0f}'.format(runNumber)+'_res.csv'
        
	simData = pd.read_csv('{}'.format(workPath) + dataFileName)
        
	#simData headers contain conflicting characters that need removing
        
	simDataHeader = str(list(simData.columns)) 
        
	chars_to_remove = ['[' , ']' , '.' , '(' , ')' , '_', '\'']
	rx = '[' + re.escape(''.join(chars_to_remove)) + ']'
	newSimDataHeader = re.sub(rx,'', simDataHeader)
        
	newSimDataHeader = newSimDataHeader.replace(" ", "")
	newSimDataHeader = newSimDataHeader.split(',')

	simData.columns = newSimDataHeader
        
	#Data analysis 
        
	time = simData['time']
	power = simData.FuelChannelNomPower
        
	indexSS = simData[simData['time']==1990].index.item()        
	index = simData[simData['time']==2000].index.item()
        
	steadyState = simData.iloc[indexSS]
	selectData = simData.iloc[index:]
	
	run.append(runNumber)
	
	maxPower.append(max(selectData.FuelChannelNomPower))
	minPower.append(min(selectData.FuelChannelNomPower))
	maxDerPower.append(max(selectData.dermPKEnpopulationn))
	maxDerFuelTemp1.append(max(selectData.derFuelChannelfuelNode1T))
	maxDerFuelTemp2.append(max(selectData.derFuelChannelfuelNode2T))
	maxDerGrapTemp.append(max(selectData.derFuelChannelgrapNodeT))
	maxInletTemp.append(max(selectData.FuelChanneltempInT))
	maxOutletTemp.append(max(selectData.FuelChanneltempOutT))
	
	SSPower.append(steadyState.FuelChannelNomPower)
	SSDerPower.append(steadyState.dermPKEnpopulationn)
	SSDerFuelTemp1.append(steadyState.derFuelChannelfuelNode1T)
	SSDerFuelTemp2.append(steadyState.derFuelChannelfuelNode2T)
	SSDerGrapTemp.append(steadyState.derFuelChannelgrapNodeT)
	SSInletTemp.append(steadyState.FuelChanneltempInT)
	SSOutletTemp.append(steadyState.FuelChanneltempOutT)

with open("senResultsY.m", "w") as text_file:
	print("run = {}".format(run)+";", file=text_file)
	print("maxPower = {}".format(maxPower)+";", file=text_file)
	print("minPower = {}".format(minPower)+";", file=text_file)
	print("maxDerPower = {}".format(maxDerPower)+";", file=text_file)
	print("maxDerFuelTemp1 = {}".format(maxDerFuelTemp1)+";", file=text_file)
	print("maxDerFuelTemp2 = {}".format(maxDerFuelTemp2)+";", file=text_file)
	print("maxDerGrapTemp = {}".format(maxDerGrapTemp)+";", file=text_file)
	print("maxInletTemp = {}".format(maxInletTemp)+";", file=text_file)
	print("maxOutletTemp = {}".format(maxOutletTemp)+";", file=text_file)
	print("SSPower = {}".format(SSPower)+";", file=text_file)
	print("SSDerPower = {}".format(SSDerPower)+";", file=text_file)
	print("SSDerFuelTemp1 = {}".format(SSDerFuelTemp1)+";", file=text_file)
	print("SSDerFuelTemp2 = {}".format(SSDerFuelTemp2)+";", file=text_file)
	print("SSDerGrapTemp = {}".format(SSDerGrapTemp)+";", file=text_file)
	print("SSInletTemp = {}".format(SSInletTemp)+";", file=text_file)
	print("SSOutletTemp = {}".format(SSOutletTemp)+";", file=text_file)
	
