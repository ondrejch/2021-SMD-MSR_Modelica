%%% Visura Pathirana
%%% MSR Physical Parameter Sensitivity 
%%% Plot correlation matrix

%Dataset 1 - Results during steady state
%Dataset 2 - Results during 100pcm step insertion

%% Save figures in the size of
figLocX=0;
figLocY=0;
figWidth=1100;
figHeight=1050;

%% Load data, organize and add labels
run('senResultsY.m');
run('senResultsX.m');

X1 = [SSPower',SSInletTemp',SSGrapTemp'...
     ,SSOutletTemp'];
 
inputVarLabel1 = {'SSnomPower','SSInletTemp','SSGrapTemp'...
     ,'SSOutletTemp'};


X2 = [maxPower',maxInletTemp',maxGrapTemp',maxOutletTemp'];

inputVarLabel2 = {'MaxNomPower','MaxInletTemp','MaxGrapTemp','MaxOutletTemp'};


Y = [CpCoefFuel',HTCcoefCore',HTCcoefPHX']; 
outputVarLabel = {'CpCoefFuel','HTCcoefCore','HTCcoefPHX'};

dataAll1 = [X1,Y];
[numData1,numInputVar1] = size(X1);

dataAll2 = [X2,Y];
[numData2,numInputVar2] = size(X2);


%% Correlation Matrix
cc1 = corrcoef([X1 Y]);

figure(1)
ccplot([X1 Y],[],[inputVarLabel1  outputVarLabel],cc1);
set(gcf,'position',[figLocX,figLocY,figWidth,figHeight])

cc2 = corrcoef([X2 Y]);

figure(2)
ccplot([X2 Y],[],[inputVarLabel2  outputVarLabel],cc2);
set(gcf,'position',[figLocX,figLocY,figWidth,figHeight])
