%This script is used for postanalysis plotting of individual files...
h=figure
%desired
load stimf.mat

COLORS=['rkmc']
bt='batch';
stimnum=1;
%number of files for plotting cluster
numfiles=1;
song_chan=1;
tet_chans=2:5;
clustnum=[1 2 3 4];
comb = [1,2;1,3;1,4;2,3;2,4;3,4];
TH=-2000;
fs=32000;
numstimfiles=length(stimf(stimnum).outind);
%act_file=floor(rand(1)*numstimfiles);
act_file=252;
rawwav='/cobain4/twarren4/g100o55/stim/1055.wav';
if(~exist('rawsong'))
        [rawsong,fs_sound]=wavread(rawwav);
end

%take an input file, w/cbin endininputfile
getfilesstr;
%read in and plot the raw data, with actual voltages on the y axis.
[data,fs,spkind,spkamp]=tetanal(files(act_file).fn,TH,song_chan,tet_chans);
%divide value by 2^15, to get value in volts, and then multiply by 1e06 to
%get actual value in microvolts.
sfct=(1e3/2^15)
spkmax=max(max(sfct*abs(data(:,tet_chans))));


ax(1)=subplot(5,1,1)
[sm,sp,t,f]=evsmooth(rawsong,44100,0.01);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
axis off;
box off;

plotcount=1;
for ii=1:4
    stimleng=length(data)/fs;
    ax(ii+1)=subplot(5,1,ii+1)
    plot(0:1/fs:(length(data)/fs-1/fs),sfct*data(:,ii+1))
    hold on;
    axis([0 stimleng -1.2*spkmax 1.2*spkmax])
    out=strcat('Ch',num2str(ii));
    if(ii==4)
        out=strcat(out,' uvolts')
    end
    YLabel(out,'Fontsize',12,'Color','k');    
    
    if (plotcount==4)
            XLabel('Time (s)','Fontsize',16);
            
            
    end
    box off;
    set(gca,'YTick',[ceil(-1.2*spkmax) 0  floor(1.2*spkmax)])
        if(plotcount~=4)
            set(gca,'XTick',[])
            set(gca,'YTick',[])
            set(gca,'xcolor','w')
            set(gca,'ycolor','w')
        end
    plotcount=plotcount+1;
end
%open the .spk file

spkamps=[];
spkts=[];
for(ii=1:numfiles)
    eval(['load -mat ',files(act_file+ii-1).fn,'.spk']);

        spkamps=[spkamps ;(spkamp)];
        spkts=[spkts ;spkt];
end    


%open the .clust file
eval(['load -mat ',files(act_file).fn,'.clust']);
%make the plot of this cluster, and others if necessary.

figure
%subplotind=[4 5 9 10 14 15]                    
for ii=1:size(comb,1)
    %subind=subplotind(ii);
    subplot(3,2,ii);cla;
    plot(-sfct*spkamps(:,comb(ii,1)),-sfct*spkamps(:,comb(ii,2)),...
                                        '.','MarkerSize',1);
    title([num2str(comb(ii,1)),' vs.' ,num2str(comb(ii,2))]);
end

for ii = 1:length(clustdat)
    tmp=clustdat(ii).comb;
    ind=find((comb(:,1)==tmp(1))&(comb(:,2)==tmp(2)));
    subplot(3,2,ind);
    hold on;
    plot(clustdat(ii).bnds(:,1),clustdat(ii).bnds(:,2),[COLORS(ii),'--']);

    x=clustdat(ii).bnds(:,1);
    y=clustdat(ii).bnds(:,2);
    insd=inpolygon(spkamps(:,tmp(1)),spkamps(:,tmp(2)),x,y);
    pp=find(insd);
    tmp=clustspk(ii).spkt
    if(isempty(tmp)==0)
        spkinds{ii}=tmp;
    end
        for jj = 1:size(comb,1)
        subplot(3,2,jj); hold on;
        plot(-sfct*spkamps(pp,comb(jj,1)),-sfct*spkamps(pp,comb(jj,2)),...
                    [COLORS(ii),'.'],'MarkerSize',1);
        box off;
        end
end
figure(h);
for ii=1:length(stimf(:,1))
    scfctor(1)=1.2-.2*(ii-1);
    scfctor(2)=1.2-.2*(ii);
    hght=-scfctor*spkmax;
    
    if(isempty(spkinds{ii})==0)
        for jj=1:4
           subplot(5,1,jj+1)
            drawyrasters(spkinds{ii}', hght,COLORS(ii));
            hold on;
            
        end
        end
    end
linkaxes(ax,'x');
linkaxes(ax(2:5));