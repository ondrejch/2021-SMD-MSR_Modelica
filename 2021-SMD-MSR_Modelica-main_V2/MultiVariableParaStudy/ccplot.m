function varargout = ccplot(data,Titles,Vnames,Nums)
%Correlation Coeffecent Plot Maker
%Makes color plot of the Absolute Values of the Correlation Coeffecent Matrix
%       -input of a Cell of multiple Matrices results in multiple plots
%
%    CCoeff = ccplot(Data,Titles,Vnames,Nums)
%INPUT::
%data     - Data Matrix
%OUTPUT::
%CCoeff - Correlation Coeffecent Matrix
%
%Optional Inputs::
%Titles - Cell Aray of Plot Titles
%Vnames - Cell Aray of Variable Names
%Nums   - Logical: Display Numeric Values (True/False)

%Author Michael E. Sharp
%December 2007
%
% Updated 13 September 2018 to force color map to run from 0 to 1
% Jamie Coble
% University of Tennessee, Knoxville



% if gcf==1; subplot(111);end
%figure('name','CorrCoeff Matrix')

    
     
    
    
    [o v] = size(data);
    if (nargin < 3 | isempty(Vnames)) & v < 10

        for i = 1:v
        Vnames(i) = {['V' num2str(i)]};
        end
    end
    
    
    if nargin < 4;
        if v <=10
            Nums = true;
        else
            Nums = false;
        end
    end
    
    
    [CCoeff PP] = corrcoef(data);

    h=bar3(abs(CCoeff));
    for i = 1:length(h)
        zdata = ones(6*length(h),4);
        k = 1;
        for j = 0:6:(6*length(h)-6)
            zdata(j+1:j+6,:) = abs(CCoeff(k,i));
            if Nums
            text(k,i,1,num2str(CCoeff(k,i),'%3.2g'),...
                'horizontalalignment','center')
            end
            k = k+1;
        end
        set(h(i),'Cdata',zdata)
    end
    view(2)
    colormap jet
    caxis([0 1])
    colorbar
    try
    set(gca,'XTickLabel',Vnames)
    set(gca,'YTickLabel',Vnames)
    catch
        xlabel('Variable Indices')
        ylabel('Variable Indices')
    end
    try
    title(Titles)
    catch
    title(['ABS Correlation Coefficient Matrix'])
    end
axis tight
axis square


%Output
if nargout >=1
varargout{1} = CCoeff;
varargout{2} = PP;
end