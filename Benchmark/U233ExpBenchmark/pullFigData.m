fig = gcf;
axObjs = fig.Children;
dataObjs = findobj(fig,'-property','YData');
dataObjs = findobj(fig,'-property','XData');
y = dataObjs(1).YData;
x = dataObjs(1).XData;

save('U233EXP8MW.mat','x','y')

