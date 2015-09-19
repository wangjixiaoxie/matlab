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

figstoplot=[4]
clear plotboxes
%example song specs.
%ac1204
%This is a single motif.
% exsong.path=avls.pvls{15}
% % exsong.fn='bk63w42_041208_1243.229.cbin';
% % exsong.bnds=[6.38 6.88]
figure
clear ax
ax(1)=subplot(311)
ax(2)=subplot(312)
ax(3)=subplot(313);

exsong(1).ax=ax(1);
exsong(1).path='/oriole6/bk61w42/ac1204/'
exsong(1).fn='bk63w42_041208_1649.512.cbin';
exsong(1).bnds=[7.90 8.75]
exsong(1).boxtoplot=[1]
exsong(1).scalfreqs=[6000 8000]


exsong(2).ax=ax(2);
exsong(2).path='/oriole6/bk61w42/ac1204/'
exsong(2).fn='bk63w42_041208_1649.512.cbin';
exsong(2).bnds=[7.90 8.75]
exsong(2).boxtoplot=[1]
exsong(2).scalfreqs=[6000 8000]


exsong(3).ax=ax(3);
exsong(3).path='/oriole6/bk61w42/ac1204/'
exsong(3).fn='bk63w42_041208_1649.512.cbin';
exsong(3).bnds=[7.90 8.75]
exsong(3).boxtoplot=[1]
exsong(3).scalfreqs=[6000 8000]


%This is two motifs, one with notched white noise, one without.
plotboxes(1).ofnm=58;
plotboxes(1).lbl=47
plotboxes(1).sfact=3;

plotboxes(2).ofnm=59;
plotboxes(2).lbl=48;
plotboxes(2).sfact=1;
% plotboxes(6).ofnm=36;
for ii=1:2
    ct=mod((ii-1),2)+1;
    plotboxes(ii).fbins=avls.fbins{ct}
    plotboxes(ii).ntnm=ct;
    plotboxes(ii).tshft=avls.tshft{ct}
    plotboxes(ii).NFFT=avls.NFFT(ct)
end

ptsps.marksize=4
ptsps.rawvl='raw'
ptsps.ntvl=1;
ptsps.indtoplot=[3 ]
sfact=1;
col{1}=[0 0 0]
col{2}=[0 0 0]
col{3}=[0.4 0.4 1]
col{4}=[0.4 0.4 1]

ptsps.col=col;
ptsps.plotextra=1
normfreqs=1;


ps.nts=[1 2]
ps.muind=[1 2 3]
ps.pvltoplot=[1 2 3]
ps.plotextra=0

%exsong.bnds
%exsong.plotlns

%figure 1 is the example song
if(find(figstoplot==1))
   
%     axnewvl(1:6)=subplot(3,3,4:9);
for ii=1:length(exsong)
    
    plotcbin(exsong(ii),plotboxes,normfreqs )
    axes(ax(1))
    plotbox([.35 .50 1500 3000],'c')
    plotbox([.57 .63 500 2500],'c')
    axes(ax(2))
    axis([.35 .50 1500 3000])
    axes(ax(3))
    axis([.59 .63 500 2500])
end
end


%avn plot
if(find(figstoplot==3))
    figure;
    [ax]=plotavn(ps, avls)
    
end

%figure 2 are histograms
if(find(figstoplot==2))
    figure;
    
    [ax]=plothistinactivtw5(avls,ps )
end
sfact=3000
if(find(figstoplot==4))
    figure;
    ptsps.ax=subplot(1,1,1);
    inactiv_rawpoints(avls,graphvals,ptsps,sfact)
   
end
    
 

