%clear all;
rast_dens=6;
binsize=.005
plotcount=1;
plotsum=1;
figure
clear stimf
%load '/cobain2/twarren2/w90o24/stim/stimexp2.mat' corpshiftednormg

load '/cobain4/twarren/g100o55/stim/stim3.mat' xlist

rawwav='/cobain4/twarren/g100o55/stim/1055.wav';
if(~exist('rawsong'))
        [rawsong,fs_sound]=wavread(rawwav);
end
axt(1:2)=subplot(rast_dens+3,1,1:2)
axcount=3;
[sm,sp,t,f]=evsmooth(rawsong,44100,0.01);
imagesc(t,f,log(abs(sp)));syn;ylim([0,1e4]);
hold on;
for ii=1:length(xlist)
        plot([xlist(ii) xlist(ii)],[0 10000],'w');hold on;
end
    set(gca,'Xtick',xlist);
    strlist={};
    strlist{14}='-'
    strlist{35}='-'
    set(gca,'Xticklabel',strlist);

 stimnames={'D:\tim\g100o55\stim\1055.wav' 'D:\tim\g100o55\stim\1055rev.wav' 'D:\tim\g100o55\stim\1055m10.wav' 'D:\tim\g100o55\stim\1055p10.wav' 'D:\tim\g100o55\stim\1055jig.wav' 'D:\tim\g100o55\stim\10055jigcomp.wav'}% 'D:\tim\w90o24e2\stim\9024alt.wav' 'D:\tim\w90o24e2\stim\9024altcomp.wav'};
%stimnames={'D:\tim\g100o55\stim\1055.wav' 'D:\tim\g100o55\stim\syl1_14_35.wav' 'D:\tim\g100o55\stim\syl2_14_35.wav' 'D:\tim\g100o55\stim\syl3_14_35.wav' 'D:\tim\g100o55\stim\syl4_14_35.wav' 'D:\tim\g100o55\stim\syl5_14_35.wav' 'D:\tim\g100o55\stim\syl6_14_35.wav'} 
 structnames={'bos' 'rev' 'm10' 'p10' 'jig' 'jic'}% 'alt' 'altc TRIAL'}
%structnames={'bos' 'dd' 'du'  'ud' 'uu' 'nd' 'nu'}; 
bt='batch';
stimleng=length(rawsong)/fs_sound;%in seconds


%calculate stimleng

binwid=1e-3;
%stimleng=ceil(stimleng/binwid); % # of bins

histbinwid=5.0e-3;

%figure(1);clf;hold on;grid on;
%figure(2);clf;hold on;grid on;
%figure(3);clf;hold on;grid on;
stimf=[];fnames=[];cnt=0;
stimf.outind=[];
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
            stimf(ind).outind=[stimf(ind).outind cnt]; 
        end
    end
   
	
end
fclose(fid);

% do for each cluster bound
NCLUST = 2%length(stimf);

for ii=1:1%:NCLUST
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
   hold on;
 
       
        
        
    edges=0:binsize:stimleng;
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
        
        maxval(stim)=(1/binsize)*max(stimf(stim).meanhist);
    end    
    maxplot=max(maxval)
    
    
    
    for (stim=1:length(stimf)  )   
        
        if(plotcount>rast_dens)
            figure;
            plotcount=1;
        end

        ax(plotcount+2)=subplot(rast_dens+3,1,plotcount+2);
        
       
        
        [axt(axcount:axcount+1),h1,h2]=plotrasters3(stimf(stim).rast,edges,(1/binsize)*stimf(stim).meanhist,stimleng);
        
        set(h1,'Color','r')
        set(h2,'Color','k')
        axes(axt(axcount))
        %axis([0 stimleng 0 trials+1]) 
        ylabel(structnames(stim),'Fontsize',16,'Color','r');
        
        axes(axt(axcount+1))
        axis([0 stimleng 0 ceil(maxplot)])
        
        if (plotcount==rast_dens)
            xlabel('Time (s)','Fontsize',16);
            ylabel('Firing Rate (Hz)','Fontsize',16);
        end
        set(gca,'YTick',[0 ceil(maxplot) 2*ceil(maxplot)])
        if(plotcount~=rast_dens)
            set(axt(axcount:axcount+1),'XTick',[])
            set(axt(axcount:axcount+1),'XTick',[])
            set(gca,'xcolor','w')
        end
        
        
        
        
        plotcount=plotcount+1;
        axcount=axcount+2;
    end 
        
end

linkaxes(axt,'x')
       

        
   
        
        
    
        
   