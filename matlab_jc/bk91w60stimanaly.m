% bk91w60stimanaly
% segment at 30,2,10000,2
cleandir4('batch',50000,500,5,5);
load /bulbul3/bk91w60/TemplatesAutolabel.mat
mk_tempf('batch',templaTarg,2,'obs0');
get_trigt2('batch',cntrngTarg,0.2,128,1,1);
label_trigs('batch','a','obs0',10000,1,1,10,30);




%     clear fv fbins bt smoothamp baspath conbins avls
%     avls.baspath='/cardinal2/StimBirds/bk91w60/'
%     avls.datdir='datasum'
% % Data location
% i=1;
%     avls.pvls{i}='1008_stimoff' % old stim parameters   
% % Changed stim parameters here - kept constant thereafter
% % 
%     avls.pvls{i}='1011_24L_24R_0msdelay_90msdur_50uA_20prctcatch' % WN on,inverted catch prctage! accident
%     avls.pvls{i}='1012_wnoff_catchofpt8' % WN off, corrected the catch prctage
%     avls.pvls{i}='1013_wnon_catchpt8' % WN on at dawn
%      avls.pvls{i}='1013_wnon_catchpt8_70uA' % WN on at dawn   
%      avls.pvls{i}='1013_wnon_catchpt5_50uA' % WN on at dawn        
%      avls.pvls{i}='1014_catchpt8_50uA'
%       avls.pvls{i}='1015_wnstillon'    
%       avls.pvls{i}='1018_wnoff'
%     avls.cvl{i}='batchfiles'
%     avls.NT{i}='b'
    
    
     clear fv fbins bt smoothamp baspath conbins avls
    avls.baspath='/bulbul3/bk91w60/'
    avls.datdir='datasum'
% Data location
i=1;
   
    avls.pvls{i}='0706_tmptest'
    avls.pvls{i}='0707_stim100prct_R34_L34_20msdel_30msdur_40uA'
    avls.pvls{i}='0708_stim20prct_R34_L34_0msdel_30msdur_40uA'    
    avls.pvls{i}='0708_stim20prct_R14_L14_0msdel_30msdur_40uA'
    avls.pvls{i}='0708_stim20prct_R23_L23_0msdel_30msdur_40uA'
    avls.pvls{i}='0708_stim50prct_R34_L34_0msdel_30msdur_80uA'    
    avls.pvls{i}='0709_stim50prct_R34_L34_0msdel_50msdur_80uA'    
    avls.pvls{i}='0710_stim20prct_R34_L34_0msdel_50msdur_80uA'    
   
    avls.cvl{i}='batch'
    avls.NT{i}='a'
    
    
%%%
% i=1;
%     avls.pvls{i}='0722_lead_in/'
%     avls.cvl{i}='batch.rand'
%     avls.NT{i}='c'



% Stim parameters
    avls.NOSTIM(i)=0;
    avls.stimDUR=.05;
    avls.del(i)=0;

% Pitch calc parameters
    avls.mkfv=[1]
    avls.con_tbinshft=0.05;
    %shift for pitch analysis
    avls.pt_tbinshft=.122;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    %shift for contours   
  %%%%% This gets changed? %%%%%
        avls.fbins=[2000 3000]  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avls.basruns=[]
    avls.contanal=1
    %1 or 0 depending set to determine whether to run an analysis
    avls.contfv=[1]
    avls.ptfv=[1]
    avls.analfv=[1]
    avls.catchstimfv=[1]
    avls.contnum=1;
    avls.con_tbinshft=0.05;
    %shift for pitch analysis
    avls.pt_tbinshft=.154;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    %shift for contours
    avls.mnbas=[]
    avls.stdbas=[]
    avls.basruns=1;
    
% Stim analysis settings    
    avls.SOUNDSTR='wav'
    avls.STIMSTR='rig'
    %window to look for stim pulse
    avls.STIMBNDS=[100 -50]
    avls.HST_EDGES=[6300:50:8000]
%     avls.STAN_RUNS=[1:7 13:15 17:19 23:25 27:36 39 41:52 54:56]
    avls.STAN_RUNS=[1]
    avls.REMOVEOUTLIERS=1;
    avls.NUMSTD=2;
    conbinsmid=[]	
    conbinshi=[]
    conbinsnorm=[]
    conbinspec=[]
avls.conbins{1}=[2000 3000]; % because harms=3
avls.conbins{2}=[2000 3000];

%%%%????    
    avls.mnbas=[7172]
    avls.stdbas=[118]
   % avls.basruns=1;

 avls=contouranal(avls);
avls=stimanal(avls);



% Summarize
% go into original folder
load batch.mat

            % align stim time onset to note onset and show pitch contours
                    % first point in pitch contour - in ms - 512/32 is for window overlap
                    fppc=(avls.con_tbinshft*1000)+512/32;
                    t=[fppc:1/8:(fppc+1/8*(size(contours,1)-1))];
                    clear t1 t2
                    for i=1:length(crfbind)
                        t1(:,i)=t;
                    end
                    for i=1:length(crctind)
                        t2(:,i)=t;
                    end
                    figure;hold on;
                            plot(t2,contours(:,crctind),'k')
                            plot(t1,contours(:,crfbind),'r')



                    for i=crfbind
                            tSTIMonset(i)=avls.del(1)*1000+fvst(i).STIMTIME;
                            tSTIMoffset(i)=tSTIMonset(i)+avls.stimDUR*1000;
                    end
                    FFrange1=max(max(contours))-min(min(contours));
                    yvals=[min(min(contours)):FFrange1/length(crfbind):max(max(contours))-FFrange1/length(crfbind)];
                    plot([tSTIMonset(crfbind);tSTIMoffset(crfbind)],[yvals;yvals],'b','Linewidth',3)
            % plot standard deviation contours
            % plot timecourse
                window=500:540;
                figure;hold on;
            plot(timing4(fvpt),mean(contours(window,:)),'*','Color','k')
            plot(timing4(fvpt(crfbind)),mean(contours(window,crfbind)),'*','Color','r')


% Group data into file that allows easy tracking of learning
        i=6;
        experiment(i).folder='0709_stim50prct_R34_L34_0msdel_50msdur_80uA'
        experiment(i).batch='batch'
        experiment(i).contours=contours;
        experiment(i).fv=fvpt;
        experiment(i).crfbind=crfbind;
        experiment(i).crctind=crctind;

%%%% PLOT
figure;hold on;
 wideCatches=[10 10 10 20 10 10 50 20];
 wideStims=[10 10 5 5 5 5 20 10];
for i=3:7%length(experiment)
            wideCatch=wideCatches(i);
            wideStim=wideStims(i);
            window=[500:550];
                                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wideCatch-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wideCatch),runningaverage(FFdata(startpoints(j):endpoints(j)),wideCatch),'*','Color','k')
                        end

                    end
                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wideStim-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wideStim),runningaverage(FFdata(startpoints(j):endpoints(j)),wideStim),'*','Color','r')
                        end
                    end

end