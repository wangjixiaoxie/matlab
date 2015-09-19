r3g74stimanaly

% segment at default
%cleandir4('batch',10000,500,6,10);
load /bulbul3/StimBirds/r3g74/TemplateWN0726.mat
mk_tempf('batch',templaWN2,2,'obs0');
get_trigt2('batch',cntrngWN2,0.2,128,1,1);
label_trigs('batch','a','obs0',10000,1,1,5,30);



    clear fv fbins bt smoothamp baspath conbins avls
    avls.baspath='/bulbul3/StimBirds/r3g74/'
    avls.datdir='datasum'
% Data location
i=1;
   % Baseline, no stim
        avls.pvls{i}='0726_R14_L14_stim50prct_50uA_50msdel_50msdur'
        avls.pvls{i}='0726_R14_L14_stim50prct_50uA_20msdel_80msdur'
        avls.pvls{i}='0726_R23_L23_stim50prct_50uA_20msdel_80msdur'
        avls.pvls{i}='0727_R23_L23_stim50prct_50uA_20msdel_80msdur'
        avls.pvls{i}='0727_R23_L23_stim50prct_50uA_20msdel_50msdur'
        avls.pvls{i}='0727_R24_L24_stim50prct_80uA_40msdel_50msdur'
        avls.pvls{i}='0727_R13_L13_stim50prct_80uA_40msdel_50msdur'
        avls.pvls{i}='0728_R23_L23_stim50prct_80uA_30msdel_70msdur'
        avls.pvls{i}='0729_R23_L23_stim50prct_120uA_20msdel_80msdur'
        avls.pvls{i}='0729_R23_L23_stim50prct_100uA_20msdel_80msdur'
        avls.pvls{i}='0729_WNON_R23_L23_stim50prct_100uA_20msdel_80msdur'
        avls.pvls{i}='0729_WNON_R34_L34_stim20prct_100uA_20msdel_80msdur'
        avls.pvls{i}='0801_WNON_R23_L23_stim50prct_100uA_20msdel_80msdur'
       
        
    avls.cvl{i}='batch'
    avls.NT{i}='a'
%%%
% i=1;
%     avls.pvls{i}='0722_lead_in/'
%     avls.cvl{i}='batch.rand'
%     avls.NT{i}='c'


% Stim parameters
    avls.NOSTIM(i)=0;
    avls.stimDUR=.08
    avls.del(i)=.02

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
        avls.fbins=[2800 3800]  
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
    avls.STIMBNDS=[300 0]
    avls.HST_EDGES=[6300:50:8000]
%     avls.STAN_RUNS=[1:7 13:15 17:19 23:25 27:36 39 41:52 54:56]
    avls.STAN_RUNS=[1]
    avls.REMOVEOUTLIERS=1;
    avls.NUMSTD=2;
    conbinsmid=[]	
    conbinshi=[]
    conbinsnorm=[]
    conbinspec=[]
avls.conbins{1}=[2800 3800]; % because harms=3
avls.conbins{2}=[2800 3800];

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
                window=350:370;
                figure;hold on;
            plot(timing4(fvpt),mean(contours(window,:)),'*','Color','k')
            plot(timing4(fvpt(crfbind)),mean(contours(window,crfbind)),'*','Color','r')

corrcoef(contours(420,crfbind),tSTIMonset(crfbind))
% positive - later ones are higher
            
% Group data into file that allows easy tracking of learning


% % 100% stim and wn on
% plot(timing4(experiment(15).fv(iscatch15)),mean(experiment(15).contours(wind1,iscatch15)),'*','Color','r')
% % 0% stim and wn off
% plot(timing4(experiment(16).fv),mean(experiment(16).contours(wind1,:)),'*','Color','b')
% plot(timing4(experiment(17).fv(experiment(17).crctind)),mean(experiment(17).contours(wind1,experiment(17).crctind)),'*','Color','b')
% plot(timing4(experiment(17).fv(experiment(17).crfbind)),mean(experiment(17).contours(wind1,experiment(17).crfbind)),'*','Color','g')
% plot(timing4(experiment(18).fv(experiment(18).crctind)),mean(experiment(18).contours(wind1,experiment(18).crctind)),'*','Color','b')
% plot(timing4(experiment(18).fv(experiment(18).crfbind)),mean(experiment(18).contours(wind1,experiment(18).crfbind)),'*','Color','g')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% figure;hold on;
% % baseline
% 
% for i=1:7 %length(experimentCatch)
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
