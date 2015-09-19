%function batchpltscat(bt,NGROUP);
% batchpltscat(bt,NGROUP);

bt='batch';
NGROUP=5;

%load up that data files structure
getfilesstr;

%some constants
comb = [1,2;1,3;1,4;2,3;2,4;3,4];

% cluster colors
NCOLOR=1;
COLORS=['rgmkcy'];

% save the cluster data
NCLUST=1;
clustdat=[];

%open up the plot
CLUSTPLT=figure;
set(CLUSTPLT,'Position',[5,162,783,779]);


YP=5;
COLORBOX=uicontrol('Style','edit','Position',[100,YP,20,20],...
        'String',COLORS(1));
MAKENEWCLUST=uicontrol('Style','pushbutton',...
        'Position',[130,YP,70,20],...
        'String','Cluster','Callback','NewCluster2');
CURRCLUSTERBOX=uicontrol('Style','edit','Position',[210,YP,30,20],...
        'String','1');
NEXTBTN=uicontrol('Style','pushbutton',...
        'Position',[250,YP,70,20],...
        'String','Next','Callback','NextFile');
GFITBTN=uicontrol('Style','pushbutton',...
        'Position',[340,YP,70,20],...
        'String','Guass Fit','Callback','GFit');
REPLOTBTN=uicontrol('Style','pushbutton',...
        'Position',[420,YP,70,20],...
        'String','Replot','Callback','Replot');
UNDOLASTCLUST=uicontrol('Style','pushbutton',...
        'Position',[500,YP,70,20],...
        'String','Undo Clust','Callback','UndoClust');

AMPPLOTS=zeros([6,1]);
for ii = 1:6
    AMPPLOTS(ii)=subplot(3,2,ii);
end

curfileind=1;
set(CLUSTPLT,'Toolbar','figure');

NextFile;
