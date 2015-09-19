r30g38stimanaly
% NOTE C - first high stack

% segment at default
cleandir4('batch',10000,500,6,10);
load /cardinal2/StimBirds/r30g38/TemplaStim1120.mat
mk_tempf('batch',templaBstimC,2,'obs0');
get_trigt2('batch',cntrngBstimC,0.2,128,1,1);
label_trigs('batch','b','obs0',10000,1,1,5,30);



    clear fv fbins bt smoothamp baspath conbins avls
    avls.baspath='/cardinal2/StimBirds/r30g38/'
    avls.datdir='datasum'
% Data location
i=1;
   % Stim test
        avls.pvls{i}='1120_stimtest_110msdel_80msdur_50uAB'
        avls.pvls{i}='1121_stim_150msdel_80msdur_50uAB'
        avls.pvls{i}='1122_wnonB'   
        avls.pvls{i}='1122_wnoffB'    
        avls.pvls{i}='1123_wnonB'     
        avls.pvls{i}='1124_wnoff_stimtestB'   
        avls.pvls{i}='1125_fixed_oscillationsB' 
        avls.pvls{i}='1126_stim80prctB'  
        avls.pvls{i}='1127_stim100prctB' 
        avls.pvls{i}='1127_stim100prct_wnonB'     
        avls.pvls{i}='1128_catchB'        
        avls.pvls{i}='1129_wnon'                
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
    avls.del(i)=.15

% Pitch calc parameters
    avls.mkfv=[1]
    avls.con_tbinshft=0.17;
    %shift for pitch analysis
    avls.pt_tbinshft=.242;
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
    avls.con_tbinshft=0.17;
    %shift for pitch analysis
    avls.pt_tbinshft=.274;
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
    avls.STIMBNDS=[100 -300]
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
load batch.rand.mat

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
        i=12
       experiment(i).folder='1129_wnon'
        experiment(i).batch='batch'
        experiment(i).contours=contours;
        experiment(i).fv=fvpt;
        experiment(i).crfbind=crfbind;
        experiment(i).crctind=crctind;

% NOTE: moved microphone before experiment(7:end) which could explain the
% apparent difference in when the contours start to look 'nice'.  Anyway,
% the window of 400-450 is good for everything.

figure;hold on;
for i=1:length(experiment)
    wide=30;
    window=[400:450];
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

