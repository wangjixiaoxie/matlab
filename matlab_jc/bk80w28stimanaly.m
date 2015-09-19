bk80w28stimanaly
% Stim off NOTE a - noisy note and Trig/measure NOTE b - low stack
% One stimulator for both sides

% load /cardinal2/StimBirds/bk80w28/templateStim1208.mat
% mk_tempf('batch.rand',templaStim,2,'obs0');
% get_trigt2('batch.rand',cntrngStim,0.2,128,1,1);
% label_trigs('batch.rand','a','obs0',10000,1,1,5,30);
load /cardinal/bk80w28/templateWN1208.mat
mk_tempf('batch',templaWN,2,'obs0');
get_trigt2('batch',cntrngWN,0.2,128,1,1);
label_trigs('batch','b','obs0',10000,1,1,5,30);



    clear fv fbins bt smoothamp baspath conbins avls
    avls.baspath='/cardinal/bk80w28/'
    avls.datdir='datasum'

i=1;
    avls.pvls{i}='1213_20msdelay_80msdur_R14_L23_80uA'
    avls.pvls{i}='1213_20msdelay_80msdur_R23_L14_80uA'    
    avls.pvls{i}='1213_20msdelay_80msdur_R14_L23_80uA'    
    avls.pvls{i}='1214_0msdelay_80msdur_R14_L23_100uA' 
    avls.pvls{i}='1214_R13_L23'
    avls.pvls{i}='1214_R23_L23' 
    avls.pvls{i}='1214_R23_L34' 
    avls.pvls{i}='1215_R23_L23'     
    avls.pvls{i}='1215_R23_L23_wnon_10prct'     
    avls.pvls{i}='1216_wnoff_stim50prct'    
    avls.pvls{i}='1216_R23_L23_stim150uA_20prct'  
    avls.pvls{i}='1217_R23_L12_stim150uA_30prct'   
    avls.pvls{i}='1217pm_R23_L23'       
    avls.pvls{i}='1218_R23_L23_100prctstim_0msdel_80msdur_150uA'      
    avls.pvls{i}='1218_wnon_stim100prct'  
    avls.pvls{i}='1219_wnoff_stim20prct'     
    avls.pvls{i}='0105_leadsin'         
    avls.pvls{i}='0106_stim50prct_dualstimulators_R23_L14_80uA_0msdel_80msdur' 
    avls.pvls{i}='0106_stim50prct_dualstimulators_R23_L23_80uA_0msdel_80msdur'   
    avls.pvls{i}='0106_stim20prct_dualstimulators_R23_L23_80uA_0msdel_80msdur'       
    avls.pvls{i}='0107_stim100prct_R23_L23' 
    avls.pvls{i}='0107_wnoff_stim20prct'    
    avls.pvls{i}='0108_100prctstim'
    avls.pvls{i}='0108_100prctstim_wnon'    
    avls.pvls{i}='0110_50prcstim_wnoff'        
    avls.pvls{i}='0110_50prctstim_200uA'  
    avls.pvls{i}='0110_24R_14L_200uA'     
    avls.pvls{i}='0111_R24_L14_80uA_0msdel_80msdur_stim20prct'   
    avls.pvls{i}='0111_R23_L23_100uA_0msdel_80msdur_stim20prct' 
    avls.pvls{i}='0113_WNON_20prctstim'  
    avls.pvls{i}='0113_WNON_20prctstim'     
    avls.pvls{i}='0115_WN_20prct_stim200uA'      
    avls.pvls{i}='0115_WN_20prct_stim200uA'  
    avls.pvls{i}='0115_WN_20prct_stim200uA'     
   %%%
      avls.pvls{i}='0119_wnoff'   
      avls.pvls{i}='0119_wnoff'   
      avls.pvls{i}='0119_wnoff'      
      avls.pvls{i}='0124_nown'  
   %%%
      avls.pvls{i}='0125_wnon'       
      avls.pvls{i}='0125_wnon'   
      avls.pvls{i}='0125_wnon'  
      avls.pvls{i}='0125_wnon' 
      avls.pvls{i}='0125_wnon' 
      avls.pvls{i}='0125_wnon'
      avls.pvls{i}='0131_wnasympt'
      avls.pvls{i}='0131_wnasympt'
      avls.pvls{i}='0131_wnasympt'    
   %%%
      avls.pvls{i}='0202_wnoff'      
      avls.pvls{i}='0202_wnoff'
      avls.pvls{i}='0202_wnoff'    
      avls.pvls{i}='0208_stim100prct'      
      avls.pvls{i}='0209_stim20prct_250uA'     
      avls.pvls{i}='0209_stim100prct'  
      avls.pvls{i}='0209_stim100prct_wnon'
      avls.pvls{i}='0210_wnoff_stim20prct'      
    avls.cvl{i}='batch'
    
    avls.NT{i}='b'
%%%
% i=1;
%     avls.pvls{i}='0722_lead_in/'
%     avls.cvl{i}='batch.rand'
%     avls.NT{i}='c'


% Stim parameters
    avls.NOSTIM(i)=0;
    avls.stimDUR=.08
    avls.del(i)=0

% Pitch calc parameters
    avls.mkfv=[1]
    avls.con_tbinshft=-0.050;
    %shift for pitch analysis
    avls.pt_tbinshft=.022;
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
    avls.con_tbinshft=-0.05;
    %shift for pitch analysis
    avls.pt_tbinshft=.054;
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
    avls.STIMBNDS=[200 -100]
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
                window=320:330;
                figure;hold on;
            plot(timing4(fvpt),mean(contours(window,:)),'*','Color','k')
            plot(timing4(fvpt(crfbind)),mean(contours(window,crfbind)),'*','Color','r')

% Group data into file that allows easy tracking of learning
load batch.mat
        i=53
        experiment(i).folder='0210_wnoff_stim20prct'
        experiment(i).batch='batch'
        experiment(i).contours=contours;
        experiment(i).fv=fvpt;
        experiment(i).crfbind=crfbind;
        experiment(i).crctind=crctind;

% NOTE: moved microphone before experiment(7:end) which could explain the
% apparent difference in when the contours start to look 'nice'.  Anyway,
% the window of 400-450 is good for everything.

figure;hold on;
for i=36:46%:length(experiment) % 27:35 for exp 1
    wide=20;
    window=[400:440];
    if size(experiment(i).crctind>wide)
        [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
        timedata=b;
        FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
        lastpoints=find(diff(timedata)>10);
        startpoints=[1 lastpoints+1];
        endpoints=[lastpoints length(timedata)];
        for j=1:length(endpoints)
            if endpoints(j)>startpoints(j)+wide-1
                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','k')
            end

        end
    end
    if size(experiment(i).crfbind>wide)
        [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
        timedata=b;
        FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
        lastpoints=find(diff(timedata)>10);
        startpoints=[1 lastpoints+1];
        endpoints=[lastpoints length(timedata)];
        for j=1:length(endpoints)
            if endpoints(j)>startpoints(j)+wide-1
                plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','r')
            end
        end
    end
end

