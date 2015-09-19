r87g80stimanalyB

% segment at default
% cleandir4('batch',10000,500,6,10);
load /cardinal2/StimBirds/r87g80/templaLowStack2_1018.mat
mk_tempf('batch2',templaLow2,2,'obs0');
get_trigt2('batch2',cntrngLow2,0.2,128,1,1);
label_trigs('batch2','b','obs0',10000,1,1,5,30);

% mk_tempf('batch2',templaLow2type1,2,'obs0');
% get_trigt2('batch2',cntrngLow2type1,0.2,128,1,1);
% label_trigs('batch2','b','obs0',10000,1,1,5,30);
% 
% mk_tempf('batch2',templaLow2type2,2,'obs0');
% get_trigt2('batch2',cntrngLow2type2,0.2,128,1,1);
% label_trigs('batch2','b','obs0',10000,1,1,5,30);


    clear fv fbins bt smoothamp baspath conbins avls
    avls.baspath='/cardinal2/StimBirds/r87g80/'
    avls.datdir='datasum'
% Data location
i=1;
   % Baseline, no stim
        avls.pvls{1}='1017_wnoff_stim50'
   % Baseline, stim
        avls.pvls{2}='1017_stimlowstack40msdelay'
        avls.pvls{3}='1018pm_40msdelay_0prctcatch'
        avls.pvls{4}='1019_40msdelay_0prctcatch_WNON'
        avls.pvls{5}='1020_test'    
        avls.pvls{6}='1020_wnoff_stimpt3'
        avls.pvls{7}='1021_nown_120msdelay'        
         avls.pvls{8}='1021pm_nown_40msdelay'      
         avls.pvls{9}='1022_stimtest_40ms'    
         avls.pvls{10}='1022pm_stim100prct_catch0'             
         avls.pvls{11}='1024pm_nostim'      
         avls.pvls{12}='1025_stimpt5'      
         avls.pvls{13}='1025pm_stim100'               
         avls.pvls{14}='1025pm_stim100_wnon'               
         avls.pvls{15}='1026pm_wnoff_stimpt2'      
         avls.pvls{16}='1027_wnon_pt2stim'     
          avls.pvls{17}='1029_wnoff'       
          avls.pvls{18}='1031_100prctstim'                 
          avls.pvls{19}='1031_100prctstim_wnon'                           
            avls.pvls{20}='1031_catch'          
           avls.pvls{21}='1101_100prctstim'             
          avls.pvls{22}='1101_100prctstim_wnon_up'  
            avls.pvls{23}='1101_catch'        
             avls.pvls{24}='1101_catch2'     
           avls.pvls{i}='1102_100prctstim'  
            avls.pvls{i}='1102_100prctstim_wnon_down'    
            avls.pvls{i}='1102_catch'   
           avls.pvls{i}='1103_100prctstim'      
             avls.pvls{i}='1103_100prctstim_wnon_up'    
             avls.pvls{i}='1103_catch'
            avls.pvls{i}='1104_100prctstim'  
               avls.pvls{i}='1104_catch'      
               avls.pvls{i}='1105_100prctstim'
               avls.pvls{i}='1105_catch' 
               
          %%% WN positive control     
            avls.pvls{i}='1110_40msdelay_20prctstim'      
            avls.pvls{i}='1110_20prctstim_wnon'
            avls.pvls{i}='1112_wnoff'            
            avls.pvls{i}='1114_100prctstim' 
            avls.pvls{i}='1114_100prctstim_wnon' 
            avls.pvls{i}='1115_catch'             
            avls.pvls{i}='1116_100prctstim'                         
            avls.pvls{i}='1117_20prctstim'   
             avls.pvls{i}='1117_20prct_stim_acausal_0msdelay'    
             avls.pvls{i}='1118_100prctstim_acausal'  
             avls.pvls{i}='1118_100prctstim_wnon'    
            avls.pvls{i}='1119_catch' 
            avls.pvls{i}='1119_stim20prct_40msdelay_catch' 
            avls.pvls{i}='1121_wnon_shiftdown_20prct' 
            avls.pvls{i}='1122_wnoff'    
              avls.pvls{i}='1124_leadsin'              
            avls.pvls{i}='1124_stim100prct'    
            avls.pvls{i}='1124_stim100prct_wnon'                
            avls.pvls{i}='1125_catch'   
             avls.pvls{i}='1126_100prctstim'    
            avls.pvls{i}='1126_100prctstim_wnon'    
            avls.pvls{i}='1127_catch'    
            avls.pvls{i}='1128_stimon_20prct'      
            avls.pvls{i}='1129_stim100prct'      
            avls.pvls{i}='1129_stim100prct_wnon'  
             avls.pvls{i}='1129_catch'    
              avls.pvls{i}='1129_R12_L34'                
              avls.pvls{i}='1130_R12_L34_stim100'                
              avls.pvls{i}='1130_catch'
              avls.pvls{i}='1201_100prctstim'   
              avls.pvls{i}='1201_100prctstim_wnon'   
              avls.pvls{i}='1202_catch'   
              avls.pvls{i}='1202_acausalstim_0msdelay_20prct'        
              avls.pvls{i}='1203_acausal_stim100prct'      
              avls.pvls{i}='1203_acausalstim100prct_wnon' 
              avls.pvls{i}='1204_catch'
              avls.pvls{i}='1204catch_causal'     
              
     %%%%%%%%%         
  %%%%%%%%%%%            
        avls.pvls{i}='1219_wnon'           
         avls.pvls{i}='1220_wn' 
         avls.pvls{i}='1220_wn'                        
        avls.cvl{i}='batch'
        
    avls.NT{i}='b'
    
    
%     clear fv fbins bt smoothamp baspath conbins avls
%     avls.baspath='/cardinal2/StimBirds/r87g80/'
%     avls.datdir='datasum'
    % Data location
%     for i=1:length(experimentREV)
%         avls.pvls{i}=experimentREV(i).folder;
%         avls.cvl{i}=experimentREV(i).batch;
%         avls.NT{i}='b';
%     end
% for i=1:length(experimentREV)
%%%
% i=1;
%     avls.pvls{i}='0722_lead_in/'
%     avls.cvl{i}='batch.rand'
%     avls.NT{i}='c'



%%%% If you change stim parameters, change these fields:
% Stim parameters
    avls.NOSTIM(i)=0;
    avls.stimDUR=.09
    avls.del(i)=0
 %%%% CRITICAL - window starting at note onset - change for diff stim parameters
    avls.STIMBNDS=[100 -40]
 %%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
 
 
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
                    t=[fppc:1/8:(fppc+1/8*(length(contours)-1))];
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
        i=66;
        experiment(i).folder='1201_100prctstim_wnon' 
        experiment(i).batch='batch'
        experiment(i).contours=contours;
        experiment(i).fv=fvpt;
        experiment(i).crfbind=crfbind;
        experiment(i).crctind=crctind;
% 2 has 10ms delay
figure;hold on;
% baseline

for i=35:38%2:length(experiment) %50:length(experiment)%18:21%[1:length(experiment)]
    wide=51;
    window=[360:400];
    if size(experiment(i).crctind>wide)
        [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
        timedata=b;
        FFdata=median(experiment(i).contours(window,experiment(i).crctind(a)));
        lastpoints=find(diff(timedata)>10);
        startpoints=[1 lastpoints+1];
        endpoints=[lastpoints length(timedata)];
        for j=1:length(endpoints)
            if endpoints(j)>startpoints(j)+wide-1
                    plot((timedata(startpoints(j):endpoints(j))),(FFdata(startpoints(j):endpoints(j))),'b.','Markersize',15)   
                    plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',6)
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
                plot((timedata(startpoints(j):endpoints(j))),(FFdata(startpoints(j):endpoints(j))),'r.','Markersize',15)
                plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',3)
            end
        end
    end

end
% 
% 
% 
% 
% figure;hold on;
% % baseline
% 
% for i=1:length(experimentCatch)
%     wide=10;
%     window=[330:400];
%     if size(experimentCatch(i).crctind>wide)
%         [b,a]=sort(timing4(experimentCatch(i).fv(experimentCatch(i).crctind)));
%         timedata=b;
%         FFdata=mean(experimentCatch(i).contours(window,experimentCatch(i).crctind(a)));
%         lastpoints=find(diff(timedata)>10);
%         startpoints=[1 lastpoints+1];
%         endpoints=[lastpoints length(timedata)];
%         for j=1:length(endpoints)
%             if endpoints(j)>startpoints(j)+wide-1
%                 if i==5 | i==6
%                     plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','b')
%                 else
%                     plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','k')
%                 end
%             end
% 
%         end
%     end
%     if size(experimentCatch(i).crfbind>wide)
%         [b,a]=sort(timing4(experimentCatch(i).fv(experimentCatch(i).crfbind)));
%         timedata=b;
%         FFdata=mean(experimentCatch(i).contours(window,experimentCatch(i).crfbind(a)));
%         lastpoints=find(diff(timedata)>10);
%         startpoints=[1 lastpoints+1];
%         endpoints=[lastpoints length(timedata)];
%         for j=1:length(endpoints)
%             if endpoints(j)>startpoints(j)+wide-1
%                 plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','r')
%             end
%         end
%     end
% 
% end
% 


load /cardinal2/StimBirds/r87g80/ExperimentB_all2.mat % second of two stack notes
% FIRST EXPERIMENT - UP
        figure;hold on;
        subplot(221);hold on;
        % baseline
        [b,a]=sort(timing4(experiment(37).fv(experiment(37).crctind)));
        starttime=min(b);
        for i=36:38%2:length(experiment) %50:length(experiment)%18:21%[1:length(experiment)]
            wide=51;
            window1=[360:400];
            if size(experiment(i).crctind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crctind(a)));
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if i==37
                        plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'r.','Markersize',15)
                    else
                        plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'b.','Markersize',15)
                    end

                    if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                    end
                end
            end
        end
        ylim([2500 2750])
        plot([-10 60],[median(mean(experiment(36).contours(window1,:))) median(mean(experiment(36).contours(window1,:)))],'b')
        xlim([-10 60])
        %
        % Chronic stim + contingent WN
        subplot(222);hold on;
        % baseline
        [b,a]=sort(timing4(experiment(40).fv(experiment(40).crfbind)));
        starttime=min(b);
        for i=38:41%2:length(experiment) %50:length(experiment)%18:21%[1:length(experiment)]
            wide=51;
            window1=[360:400];
            if size(experiment(i).crctind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crctind(a)));
                if i==38;md38=median(FFdata(116:end));end
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'b.','Markersize',15)
                    if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                    end
                end
            end
            if size(experiment(i).crfbind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crfbind(a)));
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if i==39 | i==40
                        plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'r.','Markersize',15)
                        if endpoints(j)>startpoints(j)+wide-1
                            plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                        end
                    end
                end
            end
        end
        ylim([2500 2750])
        plot([-30 50],[median(mean(experiment(39).contours(window1,:))) median(mean(experiment(39).contours(window1,:)))],'r')
        plot([-30 50],[md38 md38],'b')
        xlim([-30 50])
        %
        %
        %
        subplot(223);hold on;
        % baseline
        [b,a]=sort(timing4(experiment(42).fv(experiment(42).crfbind)));
        starttime=min(b);
        for i=41:43%2:length(experiment) %50:length(experiment)%18:21%[1:length(experiment)]
            wide=51;
            window1=[360:400];
            if size(experiment(i).crctind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crctind(a)));
                if i==41;md41=median(FFdata);end
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'b.','Markersize',15)
                    if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                    end
                end
            end
            if size(experiment(i).crfbind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crfbind(a)));
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if i==42
                        plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'r.','Markersize',15)
                        if endpoints(j)>startpoints(j)+wide-1
                            plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                        end
                    end
                end
            end
        end
        ylim([2500 2750])
        plot([-25 30],[median(mean(experiment(42).contours(window1,:))) median(mean(experiment(42).contours(window1,:)))],'r')
        plot([-25 30],[md41 md41],'b')
        xlim([-25 30])
        %
        %
        %
        subplot(224);hold on;
        % baseline
        [b,a]=sort(timing4(experiment(46).fv(experiment(46).crfbind)));
        starttime=min(b);
        for i=44:48%2:length(experiment) %50:length(experiment)%18:21%[1:length(experiment)]
            wide=51;
            window1=[360:400];
            if size(experiment(i).crctind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crctind(a)));
                if i==44;md44=median(FFdata);end
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'b.','Markersize',15)
                    if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                    end
                end
            end
            if size(experiment(i).crfbind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crfbind(a)));
                if i==44;md44s=median(FFdata);end
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if 1
                        if i==46;plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'k.','Markersize',15);end
                        if endpoints(j)>startpoints(j)+wide-1
                            plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'r','Linewidth',10)
                        end
                    end
                end
            end
        end
        ylim([2500 2750])
        plot([-25 50],[md44s md44s],'r')
        plot([-25 50],[md44 md44],'b')
        xlim([-25 50])
        
% SECOND EXPERIMENT - DOWN
        figure;hold on;
        subplot(221);hold on;
        % baseline
        [b,a]=sort(timing4(experiment(49).fv(experiment(49).crctind)));
        starttime=min(b);
        for i=48:50%2:length(experiment) %50:length(experiment)%18:21%[1:length(experiment)]
            wide=51;
            window1=[360:400];
            if size(experiment(i).crctind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crctind(a)));
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if i==49
                        plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'r.','Markersize',15)
                    else
                        plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'b.','Markersize',15)
                    end

                    if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                    end
                end
            end
        end
        ylim([2450 2750])
        plot([-10 60],[median(mean(experiment(48).contours(window1,:))) median(mean(experiment(48).contours(window1,:)))],'b')
        xlim([-10 60])
        %
        % Chronic stim + contingent WN
        subplot(222);hold on;
        % baseline
        [b,a]=sort(timing4(experiment(65).fv(experiment(65).crfbind)));
        starttime=min(b);
        for i=[62 64:67]%2:length(experiment) %50:length(experiment)%18:21%[1:length(experiment)]
            wide=51;
            window1=[360:400];
            if size(experiment(i).crctind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crctind(a)));
                if i==64;md64=median(FFdata);end
                if i==62;md62=median(FFdata);end
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'b.','Markersize',15)
                    if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                    end
                end
            end
            if size(experiment(i).crfbind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crfbind(a)));
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if i==65 | i==66
                        plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'r.','Markersize',15)
                        if endpoints(j)>startpoints(j)+wide-1
                            plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                        end
                    end
                end
            end
        end
        ylim([2500 2750])
        plot([-45 50],[median(mean(experiment(65).contours(window1,:))) median(mean(experiment(65).contours(window1,:)))],'r')
        plot([-45 50],[mean([md62 md64]) mean([md62 md64])],'b')
        xlim([-45 50])
        %
        %
        %

        %
        subplot(224);hold on;
        % baseline
        [b,a]=sort(timing4(experiment(70).fv(experiment(70).crfbind)));
        starttime=min(b);
        for i=68:72%2:length(experiment) %50:length(experiment)%18:21%[1:length(experiment)]
            wide=51;
            window1=[360:400];
            if size(experiment(i).crctind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crctind(a)));
                if i==68;md68=median(FFdata);end
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'b.','Markersize',15)
                    if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'k','Linewidth',10)
                    end
                end
            end
            if size(experiment(i).crfbind>wide)
                [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                timedata=b;
                FFdata=mean(experiment(i).contours(window1,experiment(i).crfbind(a)));
                if i==69;md69s=median(FFdata);end
                if i==68;md68s=median(FFdata);end                
                lastpoints=find(diff(timedata)>10);
                startpoints=[1 lastpoints+1];
                endpoints=[lastpoints length(timedata)];
                for j=1:length(endpoints)
                    if 1
                        if i==70;plot((timedata(startpoints(j):endpoints(j)))-starttime,(FFdata(startpoints(j):endpoints(j))),'k.','Markersize',15);end
                        if endpoints(j)>startpoints(j)+wide-1
                            plot(runningaverage(timedata(startpoints(j):endpoints(j))-starttime,wide),runningmedian(FFdata(startpoints(j):endpoints(j)),wide),'r','Linewidth',10)
                        end
                    end
                end
            end
        end
        ylim([2500 2750])
        plot([-25 50],[md69s md69s],'r')
        plot([-25 50],[md68s md68s],'g')        
        plot([-25 50],[md68 md68],'b')
        xlim([-25 50])        