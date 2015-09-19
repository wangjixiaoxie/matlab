%This plotting program is used to zoom in on individual responses in a song
%file
%plot has three columns, leftmost column is just divided into a point, with
%labels for each of the four , center column is the song, right column is
%imagevec, with red lines, to indicate syllable boundaries.

plot_dens=13
clustnum=3;
pre_tm=.5;
pos_tm=.1;
song_fs=44100;
shftfs=1000;
binsize=.005
matfile='sum.mat'
load /cobain4/twarren4/g100o55/stim/stim.mat
%stimnames={'p' 'm' 'm2' 'b' 'j' 'jc' 'j2' 'jc2'}
stimnames={'p' 'm' 'b' 'j' 'jc'}
%load -mat '/cobain4/twarren/g100o55/stim/stim.mat' imagevec rawsong xlist randvec dat_shift


load sum.mat
rawsong=corpshiftednormg{1};

mxvl=max(max(epochsum(clustnum).means));
mns=epochsum(clustnum).means;
rsp=epochsum(clustnum).rsplst;

numstim=length(mns(1,:));

for i=1:length(mns(:,1))
    plotnum=mod(i-1,plot_dens);
    if(plotnum==0)
        figure
    end
    subplot(plot_dens,2,plotnum*2+1)
    plt_tuf_hist(mns(i,2:numstim));
    axis([0 length(mns(i,:)) 0 ceil(mxvl)]);
    axis off; box off;
    
    if(plotnum==(plot_dens-1))
        axis on;
        set(gca,'xtick',1:length(mns(i,:))-1); 
        set(gca,'ytick', [0 ceil(mxvl)]);
        set(gca,'yticklabel',[0 ceil(mxvl)]);
        set(gca,'xticklabel',stimnames,'Fontsize',14); 
        YLabel('Spike Rate (Hz)','Fontsize', 14)
    end
    
    
        subplot(plot_dens,2,plotnum*2+2)
        trg_tm=(mean([rsp(i,1) rsp(i,2)])/(1/binsize));
        tm_bnds=([(trg_tm-pre_tm) (trg_tm+pos_tm)]);
        spect_bnds=round(song_fs.*tm_bnds);
        shft_bnds=round(shftfs.*tm_bnds);
     if(spect_bnds(1)>0)  
        [sm,sp,t,f]=evsmooth(rawsong(spect_bnds(1):spect_bnds(2)),song_fs,0.01);
        imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
        set(gca,'Xtick',([min(t) max(t)]));
        tm_lbl=round(tm_bnds*100)/100
        set(gca,'Xticklabel',([tm_lbl(1) tm_lbl(2)]));
        set(gca,'ycolor','w');
        hold on;
    
        %This shows the response time over the spectrogram
        xvla=[rsp(i,1)/(1/binsize)-tm_bnds(1) rsp(i,1)/(1/binsize)-tm_bnds(1)];
        yvla=[1000 8000];
        xvlb=[rsp(i,2)/(1/binsize)-tm_bnds(1)  rsp(i,2)/(1/binsize)-tm_bnds(1)];
        yvlb=[1000 8000];
        plot(xvla,yvla,'r','Linewidth',1);
        plot(xvlb,yvlb,'r','Linewidth',1);
        
        
        ind=find(randvecfin==1);
        xvals=xlist(ind);
        xvals2=xlist(ind+1);
        yarray=500;

        drawxrasters(xvals'-tm_bnds(1),xvals2'-tm_bnds(1),yarray,'r')
        axis([min(t) max(t) 0 8000])
        
        box off;
        colormap('gray');
    end
    
    end
    Xlabel('Time(s)','Fontsize',14);
%Only way I can think to debug this is to plot all the xlist points and the entire imagevec


    