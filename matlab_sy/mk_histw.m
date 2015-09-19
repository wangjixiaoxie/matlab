%clear all;

load '/cobain2/twarren2/w90o24/stim/stimexp2.mat' corpshiftednormg
rast_dens=4;
plotcount=1;
stimnames={'D:\tim\w90o24\stim\9024.wav' 'D:\tim\w90o24\stim\9024rev.wav' 'D:\tim\w90o24\stim\9024m10.wav' 'D:\tim\w90o24\stim\9024p10.wav' 'D:\tim\w90o24\stim\9024jig.wav' 'D:\tim\w90o24\stim\9024jigcomp.wav'}% 'D:\tim\w90o24e2\stim\9024alt.wav' 'D:\tim\w90o24e2\stim\9024altcomp.wav'};
structnames={'bos' 'rev' 'm10' 'p10' 'jig' 'jigc'}% 'alt' 'altc TRIAL'}

bt='batch';
stimleng=length(corpshiftednormg{1})/44100;%in seconds


%calculate stimleng

binwid=1e-3;
%stimleng=ceil(stimleng/binwid); % # of bins

histbinwid=5.0e-3;

%figure(1);clf;hold on;grid on;
%figure(2);clf;hold on;grid on;
%figure(3);clf;hold on;grid on;
stimf=[];fnames=[];cnt=0;
fid=fopen(bt,'r');
count=1;


%this loop assigns the right number of stimfs, and coutns each, and creates
%a fnames array
while (1)
	count=count+1;
    fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end
	rd=readrecf(fn);

	cnt=cnt+1;
	fnames(cnt).fn=fn;

	ind=0;
	
    for ii=1:length(stimnames)
        stimf(ii).cnt=0;
    end
    
    
    for ii = 1:length(stimnames)
		stimf(ii).outfile=stimnames{ii};
        
        if (strcmp(rd.outfile,stimnames{ii}))
			ind=ii;
			
            stimf(ind).cnt=stimf(ind).cnt+1;
		end
    end
   
	
end
fclose(fid);

% do for each cluster bound
NCLUST = 2%length(stimf);

for ii=2:2%:NCLUST
	%clear out the rasters
	for kk = 1:length(stimf)
		%stimf(kk).rast = zeros([stimf(kk).cnt,stimleng]);
		stimf(kk).cnt  = 0;
	end

	for kk = 1:length(fnames)
		rd=readrecf(fnames(kk).fn);
		cmd=['load -mat ',fnames(kk).fn,'.clust'];
		eval(cmd);
        
        for jj = 1:length(stimf)
            if (strcmp(stimf(jj).outfile,rd.outfile))
                ind=jj;
            end
        end
        stimf(ind).cnt = stimf(ind).cnt + 1;
        ind1=stimf(ind).cnt;
        if(ind1==1)
            stimf(ind).spkarray=[];
        end
        
        
        spkt=clustspk(ii).spkt;
        %spkt=ceil(spkt./binwid);
        trials=zeros(length(spkt),1);
        trials(:)=ind1;
        [spkt]=[spkt trials];
        %This line just added 1 to the right indices
        %stimf(ind).rast(ind1,spkt)=stimf(ind).rast(ind1,spkt)+1;
        stimf(ind).rast{ind1}=spkt;
        stimf(ind).spkarray=[stimf(ind).spkarray;spkt];
    
    end

    

    clr=['brkgmcy'];
    hh=figure;hold on;
 
       
        
        
    edges=0:.05:stimleng;
        %edges2=[0 startsong*32000 endsong*32000]
    for (stim=1:length(stimf))
        
                
        trials=stimf(stim).cnt;
        
        for trialnum=1:trials
            ind2=find(stimf(stim).spkarray(:,2)==trialnum);
            if(isempty(ind2))
                stimf(stim).histdist(trialnum,length(edges))=0;
            else
                stimf(stim).histdist(trialnum,:)=(histc(stimf(stim).spkarray(ind2,1),edges))
        %histdist3{stim}(trialnum,:)=(histc(spkind(spksin{clustnum}(ind(ind2)),1),edges2))
            end
        end
       
 
        ax1 = gca;
        
        stimf(stim).meanhist=mean(stimf(stim).histdist)
        maxval(stim)=20*max(stimf(stim).meanhist);
    end    
    maxplot=max(maxval)
    
    
    
    for (stim=1:length(stimf)  )   
        if(plotcount==0|plotcount>rast_dens)
            figure;
            plotcount=1;
        end

        subplot(rast_dens,1,plotcount);
        
        clear axt
        clear h1
        clear h2
        [axt,h1,h2]=plotrasters3(stimf(stim).rast,edges,20*stimf(stim).meanhist,stimleng);
        set(h1,'LineStyle','--')
        axes(axt(1))
        axis([0 stimleng 0 trials+1]) 
        ylabel(structnames(stim),'Fontsize',16,'Color','r');
        
        axes(axt(2))
        axis([0 stimleng 0 2*ceil(maxplot)])
        
        if (plotcount==rast_dens)
            xlabel('Time (s)','Fontsize',16);
            ylabel('Firing Rate (Hz)','Fontsize',16);
        end
        set(gca,'YTick',[0 ceil(maxplot) 2*ceil(maxplot)])
        if(plotcount~=rast_dens)
            set(axt(1),'XTick',[])
            set(axt(2),'XTick',[])
            set(gca,'xcolor','w')
        end
        
        
        plotcount=plotcount+1;
    end 
        
end
   
        
        
        
        
        
        
        
   
        
        
        
        
   
        
        
    
        
   