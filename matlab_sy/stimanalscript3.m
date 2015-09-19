clear bt

pvls{1}='test10/'
bt{1}='batch'
prm{1}=[70 60 100]

pvls{2}='test11/'
bt{2}='batch'
prm{2}=[70 60 200]

pvls{3}='test12/'
bt{3}='batch'
prm{3}=[70 60 100]

pvls{4}='test14/'
bt{4}='batch'
prm{4}=[80 40 100]

pvls{5}='test18/'
bt{5}='batch'
prm{5}=[90 40 100]

pvls{6}='520_100ma_40ms_bilat2/'
bt{6}='batch'
prm{6}=[100 40 100]

pvls{7}='520_90ma_40ms_bilat/'
bt{7}='batch'
prm{7}=[90 40 100]

pvls{8}='522_100ma_40ms_400Hz/'
bt{8}='batch'
prm{8}=[100 40 400]

pvls{9}='522_100ma_40ms_100Hz/'
bt{9}='batch'
prm{9}=[100 40 100]

pvls{10}='522_100ma_40ms_400Hz_2/'
bt{10}='batch'
prm{10}=[100 40 400]

pvls{11}='523_100ma_400Hz_40ms/'
bt{11}='batch'
prm{11}=[100 40 400]

pvls{12}='524_100ma_400Hz_40ms/'
bt{12}='batch'
prm{12}=[100 40 400]

pvls{13}='526_100ma_40Hz_40ms/'
bt{13}='batch'
prm{13}=[100 40 40]

pvls{14}='526_100mz_100Hz_40ms/'
bt{14}='batch26'
prm{14}=[100 40 100]

pvls{15}='528_100ma_100Hz_40ms/'
bt{15}='batch28'
prm{15}=[100 40 100]

pvls{16}='528_100ma_400Hz_40ms/'
bt{16}='batch28'
prm{16}=[100 40 400]

pvls{17}='529_100Hz_100ma_newconfig/'
bt{17}='batch'
prm{17}=[100 40 100]

pvls{18}='530_newconfig_100Hz_CORRECT/'
bt{18}='batch'
prm{18}=[100 40 100]

pvls{19}='601_400Hz_40ms_origconfig/'
bt{19}='batch'
prm{19}=[100 40 400]

pvls{20}='603_400Hz_60ms_origconfig/'
bt{20}='batch'
prm{20}=[100 60 400]

pvls{21}='stimoff515/'
bt{21}='batch16.catch.keep.rand'

pvls{22}='601_400Hz_40ms_origconfig/'
bt{22}='batch'
prm{22}=[100 40 400]

pvls{23}='601_400Hz_60ms_oldconfig/'
bt{23}='batch'
prm{23}=[100 60 400]

pvls{24}='602_100Hz_60ms_origconfig/'
bt{24}='batch'
prm{24}=[100 60 100]


pvls{25}='602_400Hz_60ms_orig_2/'
bt{25}='batch'
prm{25}=[100 60 400]

pvls{26}='602_400Hz_60ms_origconfig/'
bt{26}='batch'
prm{26}=[100 60 400]

pvls{27}='602_40Hz_60ms_origconfig/'
bt{27}='batch02'
prm{27}=[100 60 40]

pvls{28}='603_400Hz_60ms_origconfig/'
bt{28}='batch'
prm{28}=[100 60 400] 

pvls{29}='603_100Hz_60ms_origconfig/'
bt{29}='batch'
prm{29}=[100 60 100]

pvls{30}='603_100Hz_60ms_origconfig/'
bt{30}='batch'
prm{30}=[100 60 100]

pvls{31}='603_40Hz_60ms_origconfig/'
bt{31}='batch'
prm{31}=[100 60 40]

avls.pvls=pvls
pathpref='/oriole6/bk57w35/'
muinds=[1:31]
% clear fvpt valsa
for ii=1:length(muinds)
        ii
    crind=muinds(ii);
    pathvl=avls.pvls{crind}
    cmd=['cd ' pathpref pathvl]
    eval(cmd);
    btcr=bt{crind};

    tbinshft=0.079;
    NFFT=512;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    fbins=[6000 8000];
    save BINS_B NFFT fbins tbinshft
    % frequency analysis just for 'b'
    load BINS_B
    NT='a';PRENT='';PSTNT='';

    fvpt{crind}=findwnote6(btcr,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
    fv=fvpt{crind}
    save sumdata.mat fv
    valsa{crind}=getvals(fvpt{crind},1,'TRIG');
end



clear notind fbind ctind
for run_num=1:length(bt)
    crind=run_num
    pathvl=avls.pvls{crind}
    cmd=['cd ' pathpref pathvl]
    eval(cmd)
    load sumdata.mat
   notind{run_num}=[]
    fbind{run_num}=[];
    ctind{run_num}=[];
    ct_time{run_num}=[];
    fb_time{run_num}=[];
    
    fvt=fv
    
    for ii=1:length(fvt)
        if(fvt(ii).STIMTRIG)
            if(fvt(ii).STIMCATCH)
                if(fvt(ii).STIMTIME<60)
                    ctind{run_num}=[ctind{run_num} ii]
                    ct_time{run_num}=[ct_time{run_num} fvt(ii).STIMTIME]
                end
            else
                if(fvt(ii).STIMTIME<60)
                    fbind{run_num}=[fbind{run_num} ii]
                    fb_time{run_num}=[fb_time{run_num} fvt(ii).STIMTIME]
                end
            end
        else
            if(run_num<6)
               if(fvt(ii).SOUNDTRIG)
                   if(fvt(ii).SOUNDCATCH)
                       if(fvt(ii).STIMTIME<60)
                        ctind{run_num}=[ctind{run_num} ii]
                        ct_time{run_num}=[ct_time{run_num} fvt(ii).STIMTIME]
                       end
                   else
                       if(fvt(ii).STIMTIME<60)
                            fbind{run_num}=[fbind{run_num} ii]
                            fb_time{run_num}=[fb_time{run_num} fvt(ii).STIMTIME]
                       end
                       end
               else
                   notind{run_num}=[notind{run_num} ii]
               end
            else
                notind{run_num}=[notind{run_num} ii]
            end
        end
    end
                    
         
end

%make hists
hstsperfig=14;
edges=[6300:50:7500]
% for ii=5:5  2:length(muinds)
for ii=1:length(bt)
    
    
    hstoutctind{ii}=histc(valsa{ii}(ctind{ii},2),edges)
    ctln{ii}=length(valsa{ii}(ctind{ii},2));
    fbln{ii}=length(valsa{ii}(fbind{ii},2));
    hsctnrm{ii}=hstoutctind{ii}./length(hstoutctind{ii})
    
    hstoutfdind{ii}=histc(valsa{ii}(fbind{ii},2),edges)
    hsfbnrm{ii}=hstoutfdind{ii}./length(hstoutfdind{ii})
    mnctlist{ii}=mean(valsa{ii}(ctind{ii},2))
    mnfblist{ii}=mean(valsa{ii}(fbind{ii},2))
    sdctlist{ii}=std(valsa{ii}(ctind{ii},2))
    sdfblist{ii}=std(valsa{ii}(fbind{ii},2))
end
plotind=[];
for ii=1:length(bt)
    if(ctln{ii}>20&fbln{ii}>20)
        plotind=[plotind ii]
    end
end


for ii=1:length(plotind)
    
    plotpos=mod(ii, hstsperfig)
    if(plotpos==0)
        plotpos=hstsperfig
    end
    
    if(plotpos==1)
        figure;
    end
    ax(ii)=subplot(hstsperfig,1,plotpos)
    axes(ax(ii))
    crind=plotind(ii)
    
    stairs(edges,hsfbnrm{crind},'r')
    hold on;
    stairs(edges,hsctnrm{crind},'k')
    box off;
   mnfb=mnfblist{crind}
   sdfb=sdfblist{crind}
   mnct=mnctlist{crind}
   sdct=sdctlist{crind}
   plot([mnfb-sdfb mnfb+sdfb],[.7 .7],'r')
   plot([mnct-sdct mnct+sdct],[.62 .62],'k')
   plot([mean(valsa{1}(:,2)) mean(valsa{1}(:,2))], [0 0.6], 'k--') 
   ampvl=prm{crind}(1);
   tmvl=prm{crind}(2);
   hzvl=prm{crind}(3);
   combstr=[num2str(ampvl) 'ma ' num2str(tmvl) 'ms ' num2str(hzvl) 'hz']
   text(6300,0.9, combstr)
   text(7200, 0.9, ['n=' num2str(ctln{crind})], 'Color', 'k')
    text(7350, 0.9, ['n=' num2str(fbln{crind})], 'Color', 'r')
   if(plotpos>1)
       axis off;
   end
end


linkaxes(ax);




