%%% Visura Pathirana
%%% NE597 Final Project 
%%% MSR Physical Parameter Sensitivity 
%%% Linear and PCR

close all
clear
clc

%% Model select 
%reactorCond=1 for steady state
%reactorCond=2 for transient

%para=1 for Cp
%para=2 for HTC core
%para=3 for HTC PHX

reactorCond = 1;
para = 2;

%% Load data, organize and add labels
run('senResultsY.m');
run('senResultsX.m');

if reactorCond == 1
    X = [SSPower',SSInletTemp',SSFuelTemp1',SSFuelTemp2',SSGrapTemp'...
     ,SSOutletTemp'];
 
    inputVarLabel = {'SSPower','SSInletTemp','SSFuelTemp1','SSFuelTemp2','SSGrapTemp'...
     ,'SSOutletTemp'};

elseif reactorCond == 2
    X = [maxPower',minPower',maxInletTemp',maxFuelTemp1',maxFuelTemp2',...
    maxGrapTemp',maxOutletTemp',maxDerPower',...
    maxDerFuelTemp1',maxDerFuelTemp2',maxDerGrapTemp'];

    inputVarLabel = {'maxPower','minPower','maxInletTemp','maxFuelTemp1' ...
     'maxFuelTemp2','maxGrapTemp','maxOutletTemp',...
     'maxDerPower','maxDerFuelTemp1','maxDerFuelTemp2','maxDerGrapTemp'};
end

if para == 1
    Y = CpCoefFuel'; 
    outputVarLabel = 'CpCoefFuel';

elseif para == 2

    Y = HTCcoefCore'; 
    outputVarLabel = 'HTCcoefCore';

elseif para == 3

    Y = HTCcoefPHX'; 
    outputVarLabel = 'HTCcoefPHX';
end    

%% Save figures in the size of
figLocX=0;
figLocY=0;
figWidth=1500;
figHeight=1050;

%% Dividing dataset to three groups for training, validation and testing

data = [X,Y];
[numData,numInputVar] = size(X);

trainIndx = 2:2:numData;
valIndx = 1:4:numData;
testIndx = 3:4:numData;

trainData = data(trainIndx,:);
valData = data(valIndx,:);
testData = data(testIndx,:);

  
%% Test if minimum and maximum data points are in train dataset
%Train dataset must include maximum and minimum valus
if sum(min(trainData)>min(data))>0

    indx = find(min(trainData)>min(data));
    indxMin = zeros(1,numel(indx));

    for i = 1:numel(indx)
        indxMin(i) = find(data(:,indx(i)) == min(data(:,indx(i))),1);
    end
    clear indx
    trainData = [trainData;data(unique(indxMin),:)];  
end

if sum(max(data)>max(trainData))>0

    indx = find(max(data)>max(trainData));
    indxMax = zeros(1,numel(indx));
    for i = 1:numel(indx)
        indxMax(i) = find(data(:,indx(i)) == max(data(:,indx(i))),1);
    end
    clear indx
    trainData = [trainData;data(unique(indxMax),:)];  
end

%% Rearrange dataset by seperating input data and output data
x_train = trainData(:,1:end-1);
y_train = trainData(:,end);
x_val = valData(:,1:end-1);
y_val = valData(:,end);
x_test = testData(:,1:end-1);
y_test = testData(:,end);

%% Linear Regression

RMSE_single = NaN(1,numInputVar);
R2_single = NaN(1,numInputVar);

for i = 1:numInputVar
    xt = [x_train(:,i) ones(size(x_train,1),1)];
    xval = [x_val(:,i) ones(size(x_val,1),1)];
    B = regress(y_train,xt); % simple linear regression
    y_pred = xval*B;        % predictions for validation data
    RMSE_single(i) = sqrt(mean((y_val-y_pred).^2));
    R2_single(i) = 1 - sum((y_val-y_pred).^2)/sum((y_val-mean(y_val)).^2);
end

% plot results %
figure(1);
bar(RMSE_single)
xlabel('Input Variable')
ylabel('RMSE')
set(gca,'XTickLabel',inputVarLabel,'XTickLabelRotation',45)
set(gcf,'position',[figLocX,figLocY,figWidth,figHeight])

RMSE_Linear = min(RMSE_single(:))
varSingle = find(RMSE_single == min(RMSE_single(:)));


figure(2)
plot(x_train(:,varSingle),y_train,'.')
xlabel(inputVarLabel{varSingle}) 
ylabel(outputVarLabel)
set(gcf,'position',[figLocX,figLocY,figWidth,figHeight])

figure(3)
plot(min([y_val; y_pred]):max([y_val;y_pred]),min([y_val; y_pred]):max([y_val;y_pred]),'k--')
hold on
plot(y_val,y_pred,'.')
xlabel(['Actual ' outputVarLabel])
ylabel(['Predicted ' outputVarLabel])
set(gcf,'position',[figLocX,figLocY,figWidth,figHeight])

%% PCR

[xs,x_m,x_std] = zscore1(x_train);
[loadings,~,latent,~,explained] = pca(xs);
cum_explained = cumsum(explained); 
xs_val = zscore1(x_val,x_m,x_std);

xPC = [xs*loadings ones(size(x_train,1),1)]; % add column of ones for non-zero y-intercept %

condNum = cond(xPC'*xPC)

b = regress(y_train,xPC);

xPC_val = [xs_val*loadings ones(size(x_val,1),1)];
y_pred = xPC_val*b;

RMSE_PCR = sqrt(mean((y_val - y_pred).^2))

% Plot actual vs predicted to look for trends %
figure(4); 
plot(min([y_val; y_pred]):max([y_val;y_pred]),min([y_val; y_pred]):max([y_val;y_pred]),'k--')
hold on
plot(y_val,y_pred,'.')
xlabel(['Actual ' outputVarLabel])
ylabel(['Predicted ' outputVarLabel])
set(gcf,'position',[figLocX,figLocY,figWidth,figHeight])


