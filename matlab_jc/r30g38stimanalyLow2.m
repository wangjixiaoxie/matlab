r30g38stimanaly
% NOTE C - first high stack

% segment at default
cleandir4('batch',100000,500,6,10);
load /cardinal2/StimBirds/r30g38/TemplaStim1120.mat
mk_tempf('batch2',templaBstimC,2,'obs0');
get_trigt2('batch2',cntrngBstimC,0.2,128,1,1);
label_trigs('batch2','a','obs0',100000,1,1,5,30);
changelabel('batch2','-', 'a', '', 'b');

% load /cardinal2/StimBirds/r30g38/templaLowWN_1203.mat
% mk_tempf('batch2.rand',templaLow2,2,'obs0');
% get_trigt2('batch2.rand',cntrngLow2,0.2,128,1,1);
% label_trigs('batch2.rand','f','obs0',100000,1,1,5,30);
% % changelabel('batch.rand','-', 'a', '', 'b');
% vals=evtaf_freq('batch2.rand',[2000,3000],'b',128,'obs0',1,1);

    clear fv fbins bt smoothamp baspath conbins avls
    avls.baspath='/cardinal2/StimBirds/r30g38/'
    avls.datdir='datasum'
% Data location
i=1;
   % Stim test

        avls.pvls{i}='1130_catch'  
        avls.pvls{i}='1130_Lowstack_stim40msdel_80msdur'  
        avls.pvls{i}='1201_wnon'          
        avls.pvls{i}='1201_wn_stim20prct'     
        avls.pvls{i}='1202_wnoff'  
        avls.pvls{i}='1202_24L_13R_80uA'    
        avls.pvls{i}='1202_stimoff'  
        avls.pvls{i}='1204_100prctstim'     
        avls.pvls{i}='1204_100prctstim_wnon'  
        avls.pvls{i}='1205catch'    
        avls.pvls{i}='1208_wnon'     
        avls.pvls{i}='1213_asymptote_wn'      
         avls.pvls{i}='1214wn'    
         avls.pvls{i}='1216wn'    
         avls.pvls{i}='1218_wnoff' 
         avls.pvls{i}='1220_wnoff_stim20'   
         avls.pvls{i}='1220_wnoff_stim20'          
         avls.pvls{i}='0107_'   
         avls.pvls{i}='0108_stim100prct'            
         avls.pvls{i}='0108_stim100prct_wnon'      
         avls.pvls{i}='0109_stim20prct_wnoff'     
         avls.pvls{i}='0110_R14_L23_20prct'    
         avls.pvls{i}='0110_50prctstim_R14_L23_250uA'
         avls.pvls{i}='0111_100prctstim'     
         avls.pvls{i}='0111_20prctstim'
         avls.pvls{i}='0113_WNON_20prctstim_R14_L23_250uA_40msdel_80msdur'
         
      %%%%%%%%
      %%%%%%%%
      %%%%%%%%
         avls.pvls{i}='0211_stimtest_23R_23L_200uV_40msdel_80msdur'  
         avls.pvls{i}='0211_50prctstimtest_23R_23L_100uV_40msdel_80msdur'           
         avls.pvls{i}='0211_20prctstimtest_23R_23L_75uV_40msdel_80msdur' 
         avls.pvls{i}='0211_20prctstimtest_23R_23L_75uV_40msdel_80msdur'          
         avls.pvls{i}='0211_20prctstimtest_23R_23L_75uV_40msdel_80msdur'
          avls.pvls{i}='0214_100prctstim'   
          avls.pvls{i}='0214_100prct_wnon'     
          avls.pvls{i}='0215_wnoff_stim50prct'     
          avls.pvls{i}='0215_100uA_stim20prct'
          avls.pvls{i}='0215_100uA_stim20prct' 
          avls.pvls{i}='0216_100prctstim'  
          avls.pvls{i}='0216_100prct_wnon'   
          avls.pvls{i}='0217_20prct_wnoff'
          avls.pvls{i}='0217_20prct_wnoff'  
          avls.pvls{i}='0220_wnon'     
          avls.pvls{i}='0220_wnon'     
          avls.pvls{i}='0220_wnon' 
          avls.pvls{i}='0220_wnon' 
          avls.pvls{i}='0222wn'  
          avls.pvls{i}='0222wn'            
    avls.cvl{i}='batch2'
    avls.NT{i}='b'
%%%
% i=1;
%     avls.pvls{i}='0722_lead_in/'
%     avls.cvl{i}='batch.rand'
%     avls.NT{i}='c'


% Stim parameters
    avls.NOSTIM(i)=0;
    avls.stimDUR=.08
    avls.del(i)=.04

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
load batch2.mat
        i=46
        experiment(i).folder='0222wn'
        experiment(i).batch='batch2'
        experiment(i).contours=contours;
        experiment(i).fv=fvpt;
        experiment(i).crfbind=crfbind;
        experiment(i).crctind=crctind;
%
% load experimentB_1216saved.mat
% experiment(10).crctind=[experiment(10).crctind(1:end-591) experiment(13).crctind+length(experiment(10).fv)-761];
% experiment(10).crfbind=[experiment(10).crfbind(1:end-164) experiment(13).crfbind+length(experiment(10).fv)-761];
% experiment(10).contours=[experiment(10).contours(:,1:end-761) experiment(13).contours];
% experiment(10).fv=[experiment(10).fv(1:end-761) experiment(13).fv];
%
% NOTE: moved microphone before experiment(7:end) which could explain the
% apparent difference in when the contours start to look 'nice'.  Anyway,
% the window of 400-450 is good for everything.


% Sampled 100% through wn on
% Sampled 30% afterwards

figure;hold on;
for i=40:length(experiment)%2:5%5:9%17:length(experiment)%[9:10 14:15]%length(experiment)
    wide=20;
    window=[285:300];
    if size(experiment(i).crctind>wide)
        [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
        timedata=b;
        FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
        lastpoints=find(diff(timedata)>10);
        startpoints=[1 lastpoints+1];
        endpoints=[lastpoints length(timedata)];
        for j=1:length(endpoints)
            if endpoints(j)>startpoints(j)+wide-1
                if i==5 | i==6
                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','b')
                else
                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','k')
                end
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

