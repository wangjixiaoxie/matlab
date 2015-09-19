%figure1 master script.

%first 
%option 1, plot time course for upshift and downshift,only
%runs up until asymptote...



%option 2, plot six panels
%of example reversion with raw data and histograms

%option 3, plot all runs, with controls and baseline
%on scatter plot.

%option 4, plot all offsets as percent reversion.

function [avl_struct]=inactiv_revfig1(bs,avlsind)


%figstoplot
figstoplot=[2]
avlsind=[8 5]

avl_struct=mk_avl_struct(sumbs,avlsind);


% %option 1, 
% call, plot_tmcourse2

%option 2, 
% call some variant of inactivexbk61w42fn.m
%function below, plotexamples.

if(ismember(2,figstoplot))
   
    [axrawpts,axhst]=plotexamples(avl_struct);
end

function [avl_struct]=mk_avl_struct(bs,avlsind)
    for ii=1:length(avlsind)
        crind=avlsind(ii);
        crbs=bs(crind);
        pth=crbs.path;
        bt=crbs.matfilename;
        cmd=['cd ' pth]
        eval(cmd);
        cmd=['load ' bt]
        eval(cmd);
        avl_struct(ii)=avls;
    end


function [axrawpts,axhst]=plotexamples(avl_struct)
    figure;
    rowfig=5;
    colfig=2;
    axrawind=[2 6 10]
    axhstind=[1 3 5 7 9]
    avlsind=[1 1 2]
    ptsps.marksize=14
    ptsps.rawvl='raw'
    ptsps.ntvl=1;
    indtoplot=[5 6 7]
    col{1}=[0 0 0]
    col{2}=[0 0 0]
    col{3}=[0.4 0.4 1]
    col{4}=[0.4 0.4 1]
    ptsps.plotextra=1
    ptsps.col=col;
    
    for ii=1:length(axrawind)
        subplotind=axrawind(ii);
        ptsps.indtoplot=indtoplot(ii);
        ptsps.ax=subplot(rowfig,colfig,subplotind);
        avls=avls_struct{ii}
        axhst(ii)=inactiv_rawpoints(avls,graphvals,ptsps)
    end