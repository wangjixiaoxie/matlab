%inactivexbk61w42.m.m
%possibility to make four distinctive figures here.

%figure 1 is spectrogram.
%with pitch measurement marks.

%figure 2 are histograms

%figure 3 is time course of inactivation run,
%comparison of average pitch measurements.

%raw data panel is made with calls to inactivplotoneday.m
%top pael pulled off illustrator file.
%to make acsf day 1, target note.
%ps is short for plotstruct.

figstoplot=[2]
clear plotboxes
%example song specs.
%ac1204
%This is a single motif.
% exsong.path=avls.pvls{15}
% % exsong.fn='bk63w42_041208_1243.229.cbin';
% % exsong.bnds=[6.38 6.88]

exsong.path=avls.pvls{15}
exsong.fn='bk63w42_041208_1649.512.cbin';
exsong.bnds=[6.35 8.7]
exsong.boxtoplot=[4:6]
exsong.scalfreqs=[6000 8000]

%This is two motifs, one with notched white noise, one without.

plotboxes(1).ntnm=1;

plotboxes(1).ofnm=58;
plotboxes(1).lbl=47
for ii=1:length(plotboxes)
    ct=mod((ii-1),3)+1;
    plotboxes(ii).fbins=[2150 2650]
    plotboxes(ii).ntnm=ct;
    plotboxes(ii).tshft=.03-512/32000
    plotboxes(ii).NFFT=.08
end

ptsps.marksize=14
ptsps.rawvl='raw'
ptsps.ntvl=1;
ptsps.indtoplot=3
ptsps.plotextra=1;
ptsps.sfact=3;
col{1}=[0 0 0]
col{2}=[0 0 0]
col{3}=[0.4 0.4 1]
col{4}=[0.4 0.4 1]

ptsps.col=col;
normfreqs=1;


ps.nts=[1 2]
ps.muind=[3 4 12]
ps.pvltoplot=[13 14 16 17 32 33]

%exsong.bnds
%exsong.plotlns

%figure 1 is the example song
if(find(figstoplot==1))
    ax=figure
    exsong.ax=ax;
%     axnewvl(1:6)=subplot(3,3,4:9);
    [ax]=plotcbin(exsong,plotboxes,normfreqs )
end

%avn plot
if(find(figstoplot==3))
    figure;
    [ax]=plotavn(ps, avls)
end

%figure 2 are histograms
if(find(figstoplot==2))
    figure;
    
    [ax]=plothistinactivtw(avls,ps )
end

if(find(figstoplot==4))
    figure;
    ptsps.ax=subplot(1,1,1);
    inactiv_rawpoints(avls,graphvals,ptsps,ptsps.sfact)
end
    
 

