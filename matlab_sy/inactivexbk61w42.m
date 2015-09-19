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
%test
figstoplot=[2]

%example song specs.
%ac1204
%This is a single motif.
% exsong.path=avls.pvls{15}
% % exsong.fn='bk63w42_041208_1243.229.cbin';
% % exsong.bnds=[6.38 6.88]
% 
% exsong.path=avls.pvls{16}
% exsong.fn='bk63w42_051208_1706.122.cbin';
% exsong.bnds=[5.7 7.3]
% exsong.boxtoplot=[4:6]
% exsong.scalfreqs=[6000 8000]

%This is two motifs, one with notched white noise, one without.

plotboxes(1).ntnm=1;
plotboxes(2).ntnm=2;
plotboxes(3).ntnm=3;
plotboxes(4).ntnm=1;
plotboxes(5).ntnm=2;
plotboxes(6).ntnm=3;

plotboxes(1).ofnm=27;
plotboxes(2).ofnm=28;
plotboxes(3).ofnm=29;

plotboxes(4).ofnm=34;
plotboxes(5).ofnm=35;
plotboxes(6).ofnm=36;
for ii=1:6
    ct=mod((ii-1),3)+1;
    plotboxes(ii).fbins=avls.fbins{ct}
    plotboxes(ii).ntnm=ct;
    plotboxes(ii).tshft=avls.tshft{ct}
    plotboxes(ii).NFFT=avls.NFFT(ct)
end

ptsps.marksize=14
ptsps.rawvl='raw'
ptsps.ntvl=1;
ptsps.indtoplot=5
col{1}=[0 0 0]
col{2}=[0 0 0]
col{3}=[0.4 0.4 1]
col{4}=[0.4 0.4 1]

ptsps.col=col;
normfreqs=1;


ps.nts=[1 2]
ps.muind=[1 4 5]
ps.pvltoplot=[13 14 16 17 32 33]

%exsong.bnds
%exsong.plotlns

%figure 1 is the example song
if(find(figstoplot==1))
    figure
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
    
    [ax]=plothistinactivtw3(avls,ps )
end

if(find(figstoplot==4))
    figure;
    ptsps.ax=subplot(1,1,1);
    ptsps.plotextra=1
    inactiv_rawpoints(avls,graphvals,ptsps)
end
    
 

