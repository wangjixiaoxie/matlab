oaSTIManaly
% Look at contouranalbk48w74 for instructions - run all code


% Transfer data from bigbird to the base path in swift 
% Label song

% 1. change avls.baspath to the base path in swift
% 2. change avls.pvls{1} to the folder where you have labeled song
% 3. run code from lines 14-100
% 4. go back into the folder where you have labeled song
% 5. run code from lines 102-end

                clear fv fbins bt smoothamp baspath conbins avls
                avls.baspath='/cardinal3/bk48w74/'
                avls.datdir='datasum'
            % Data location
            i=1;
                avls.pvls{i}='727pmtest' % for the folder you want to look at
                dirf('*.cbin','batch') % create a batch file
                evsonganaly % go through each file and label 'c' only the flat stacks preceded by a red triangle
                randsamp('batch',0.2)
                evsonganaly % go through each randomly selected file and label 'c' all the flat stacks

                avls.cvl{i}='batch'
                avls.NT{i}='c'
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
                avls.contnum=3;
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
            avls.conbins{1}=[2000 2700]; % because harms=3
            avls.conbins{2}=[2000 2700];

            %%%%????    
                avls.mnbas=[7172]
                avls.stdbas=[118]
               % avls.basruns=1;



        avls=contouranal(avls);
        avls=stimanal(avls);
        %
        
        
        
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
        figure;hold on;
            plot(t2,mean(contours(:,crctind)'),'k')
            plot(t1,mean(contours(:,crfbind)'),'r')


        for i=crfbind
                tSTIMonset(i)=avls.del(1)*1000+fvst(i).STIMTIME;
                tSTIMoffset(i)=tSTIMonset(i)+avls.stimDUR*1000;
        end
        FFrange1=max(max(contours))-min(min(contours));
        yvals=[min(min(contours)):FFrange1/length(crfbind):max(max(contours))-FFrange1/length(crfbind)];
        plot([tSTIMonset(crfbind);tSTIMoffset(crfbind)],[yvals;yvals],'b','Linewidth',3)
% plot standard deviation contours
    figure;hold on;
    plot(std(contours(:,crfbind)'),'r')
    plot(std(contours(:,crctind)'),'k')
% plot timecourse
    window=300:400;
    figure;hold on;
    plot(timing4(fvpt(crfbind)),mean(contours(window,crfbind)),'*','Color','r')
    plot(timing4(fvpt(crctind)),mean(contours(window,crctind)),'*','Color','k')