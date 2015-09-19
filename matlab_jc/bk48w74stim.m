bk48w74
% Look at contouranalbk48w74 for instructions - run all code

    clear fv fbins bt smoothamp baspath conbins avls
    avls.baspath='/cardinal3/StimBirds/bk48w74/'
    avls.datdir='datasum'
% Data location
i=1;
    avls.pvls{i}='0729_90prctCatch_baseline_14L_14R'
    avls.pvls{i}='0730_wnon_90prctCatch_14L_14R'
    avls.pvls{i}='0804_wnstillon_90prctCatch_14L_14R'
     avls.pvls{i}='0804_wnoff_stimon_90prctCatch_14L_14R'   
     avls.pvls{i}='0808pm_50msWN_80msSTIM_nodelay_14L_14R'
      avls.pvls{i}='0809_stimpre_and_dur'   
      avls.pvls{i}='0810_wnoff'
            avls.pvls{i}='0812_repeat0809exp'
            avls.pvls{i}='0812_switchstims_revtest'
            avls.pvls{i}='0812_wnoff_90prctcatch'
avls.pvls{i}='0816_baseline'
avls.pvls{i}='0816pm_wnon_14L_14R_downshift'
avls.pvls{i}='0818pm_wnstillon_34L_14R'
avls.pvls{i}='0819_wnoff'
avls.pvls{i}='DOWNSHIFT821/0821_baseline_23LEFT_14RIGHT'
avls.pvls{i}='DOWNSHIFT821/0822_wnon_23L_14R'
avls.pvls{i}='DOWNSHIFT821/0822_34L_14R'
    avls.cvl{i}='batchA'
    avls.NT{i}='b'
%%%
% i=1;
%     avls.pvls{i}='0722_lead_in/'
%     avls.cvl{i}='batch.rand'
%     avls.NT{i}='c'


% Stim parameters
    avls.NOSTIM(i)=0;
    avls.stimDUR=.06
    avls.del(i)=.065

% Pitch calc parameters
    avls.mkfv=[1]
    avls.con_tbinshft=.05;
    %shift for pitch analysis
    avls.pt_tbinshft=.072;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    %shift for contours   
  %%%%% This gets changed? %%%%%
        avls.fbins=[6000 8100]  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    avls.basruns=[]
    avls.contanal=1
    %1 or 0 depending set to determine whether to run an analysis
    avls.contfv=[1]
    avls.ptfv=[1]
    avls.analfv=[1]
    avls.catchstimfv=[1]
    avls.contnum=1;
    avls.con_tbinshft=.02;
    %shift for pitch analysis
    avls.pt_tbinshft=.074;
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
    avls.STIMBNDS=[200 0]
    avls.HST_EDGES=[6300:50:8000]
%     avls.STAN_RUNS=[1:7 13:15 17:19 23:25 27:36 39 41:52 54:56]
    avls.STAN_RUNS=[1]
    avls.REMOVEOUTLIERS=1;
    avls.NUMSTD=2;
    conbinsmid=[]	
    conbinshi=[]
    conbinsnorm=[]
    conbinspec=[]
avls.conbins{1}=[6000 8100]; % because harms=3
avls.conbins{2}=[6000 8100];

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
                window=300:400;
                figure;hold on;
            plot(timing4(fvpt),mean(contours(window,:)),'*','Color','k')
            plot(timing4(fvpt(crfbind)),mean(contours(window,crfbind)),'*','Color','r')


% Group data into file that allows easy tracking of learning
        i=3;
           experiment(i).folder='DOWNSHIFT821/0822_34L_14R'
        experiment(i).contours=contours;
        experiment(i).fv=fvpt;
        experiment(i).crfbind=crfbind;
        experiment(i).crctind=crctind;
%
% countFB=0;countCT=0;indFB=[];indCT=[];
% for i=1:467
% if sum(crfbind==i)>0 & std(contours(150:200,i))<120
% countFB=countFB+1;
% indFB(countFB)=i;
% else
% if sum(crctind==i)>0 & std(contours(150:200,i))<120
% countCT=countCT+1;
% indCT(countCT)=i;
% end
% end
% end
%
        figure;hold on;
        % baseline
            i=1;
            wide=15;
            window=[300:350];
                    a=[21:length(experiment(i).fv)];
                    plot(runningaverage(timing4(experiment(i).fv(experiment(i).crfbind(experiment(i).crfbind>20))),wide),runningaverage(mean(experiment(i).contours(window,experiment(i).crfbind(experiment(i).crfbind>20))),wide),'r.','Markersize',15)
                    plot(runningaverage(timing4(experiment(i).fv(experiment(i).crctind(experiment(i).crctind>20))),wide*3),runningaverage(mean(experiment(i).contours(window,experiment(i).crctind(experiment(i).crctind>20))),wide*3),'b.','Markersize',15)
        % wn on  
        for i=[2:3]% 7 14]
                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'r.','Markersize',15)
                        end
                    end
                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide*3),runningaverage(FFdata(startpoints(j):endpoints(j)),wide*3),'b.','Markersize',15)
                        end

                    end
        end
        % wn off
        for i=[4:5 8 11 12 13]
                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'r.','Markersize',15)
                        end
                    end
                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'b.','Markersize',15)
                        end

                    end
        end
                for i=[6 9 10]
                    [b,a]=sort(timing4(experiment(i).fv));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'b.','Markersize',15)
                        end
                    end
        end
%
% 1st harmonic - downshift
figure;hold on;
        % wn on  
        
             wide=10;
            window=[180:220];       
        for i=[16:19]
                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crfbind)));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,experiment(i).crfbind(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','r')
                        end
                    end
                    [b,a]=sort(timing4(experiment(i).fv(experiment(i).crctind)));
                    timedata=b;
                    FFdata=mean(experiment(i).contours(window,experiment(i).crctind(a)));
                    lastpoints=find(diff(timedata)>10);
                    startpoints=[1 lastpoints+1];
                    endpoints=[lastpoints length(timedata)];
                    for j=1:length(endpoints)
                        if endpoints(j)>startpoints(j)+wide-1
                        plot(runningaverage(timedata(startpoints(j):endpoints(j)),wide),runningaverage(FFdata(startpoints(j):endpoints(j)),wide),'*','Color','b')
                        end

                    end
        end
