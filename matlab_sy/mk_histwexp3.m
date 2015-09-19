%this is to deal with data from the ra inactivation experiment.
%this is for quick look at thresholded data before clustering...which I
%need to do.


%mk_histwexp2.m
% This code puts individual clust files into stimf structure, which can
% then be analyzed.
%commented 7.16.06
%necessary input is stimulus names
%instead of using clusts, I am going to adopt this code, so that 
%the first column of stimf(XXX,.) is the number of the input channel

%clear all;
rast_dens=6;
binsize=.005
plotcount=1;
plotsum=1;
fs_sound=44100;
%clear stim
clusts=[1 2]

clustnum=length(clusts);
stimf=[];
%load '/cobain2/twarren2/w90o24/stim/stimexp2.mat' corpshiftednormg

%load '/cobain4/twarren4/g100o55/stim/stim.mat' xlist corpshiftednormg

%this needs to be updated.
rawwav='/cobain6/twarren6/pk96g69/pk96g69.wav';
if(~exist('rawsong'))
        [rawsong,fs_sound]=wavread(rawwav);
end

%TO TEST THIS I NEED TO FIX STIMNAMES!!!, so that it is in above format.

 stimnames={'D:\tim\pk96g69\stimuli\pk96g69.wav' 'D:\tim\pk96g69\stimuli\pk96g69rev.wav'};
 %stimnames={'D:\tim\g100o55e2\stim\1055.wav' 'D:\tim\g100o55e2\stim\1055rev.wav' 'D:\tim\g100o55\stim\1055m10.wav' 'D:\tim\g100o55\stim\1055p10.wav'  'D:\tim\g100o55\stim\1055jig.wav' 'D:\tim\g100o55\stim\10055jigcomp.wav'}
 structnames={'bos' 'rev'}
%structnames={'bos' 'dd' 'du'  'ud' 'uu' 'nd' 'nu'}; 
bt='batch';
stimleng=length(rawsong)/fs_sound;%in seconds

%calculate stimleng

binwid=1e-3;
%stimleng=ceil(stimleng/binwid); % # of bins

histbinwid=5.0e-3;
stimf=[];fnames=[];cnt=0;
stimf.outind=[];


fid=fopen(bt,'r');
count=1;

%this loop just initializes
for clustcnt=1:clustnum
    ii=clusts(clustcnt);
    for jj=1:length(stimnames)
        stimf(ii,jj).cnt=0;;
        stimf(ii,jj).spkarray=[];
end
end

%this loop assigns the right number of stimfs, and coutns each, and creates
%a fnames array, fnames array just gives access to the raw name of each
%file in the batch file
lpcnt=1;
while (1)
	count=count+1;
    fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end
	rd=readrecf(fn);

	cnt=cnt+1;
	fnames(cnt).fn=fn;
    
    for ii = 1:length(stimnames)
		stimf(1,ii).outfile=stimnames{ii};
        
        %This should only be the case
        if (strcmp(rd.outfile,stimnames{ii}))
			ind=ii;
			 stimtype(lpcnt)=ii;
            stimf(1,ind).cnt=stimf(1,ind).cnt+1;
            stimf(1,ind).outind=[stimf(1,ind).outind cnt]; 
            break;
        
        else
            stimtype(lpcnt)=0;
        end
        end
   
    lpcnt=lpcnt+1;
end
fclose(fid);

for ii=1:clustnum
    for jj=1:length(stimnames)
        stimf(ii,jj).cnt=stimf(1, jj).cnt;
        stimf(ii,jj).outfile=stimf(1,jj).outfile;
        stimf(ii,jj).outind=stimf(1,jj).outind;
    end
end
% do for each cluster bound
NCLUST = 2%length(stimf);


%critical loop which takes the clustered data and puts into the stimf
%structure, which is used for later analysis  (rasters, etc...)
for clustval=1:clustnum
	%clear out the rasters
	
    for kk = 1:length(fnames)
		rd=readrecf(fnames(kk).fn);
		%this is the really critical line which reads the .clust file
        cmd=['load -mat ',fnames(kk).fn,'.clust'];
        eval(cmd);

        stvl=stimtype(kk);
        if(stvl)
            numstm=find(stimf(clustval,stvl).outind==kk);
            spkt=clustspk(clustval).spkt;
            if(isempty(spkt))
                spkt=-10;
            end
        trials=zeros(length(spkt),1);
        trials(:)=numstm;
        [spkt]=[spkt trials];
        stimf(clustval,stvl).rast{numstm}=spkt;
        stimf(clustval,stvl).spkarray=[stimf(clustval,stvl).spkarray;spkt];
        end
    end

   clr=['brkgmcy'];
   hold on;
    edges=0:binsize:stimleng;
        %edges2=[0 startsong*32000 endsong*32000]
    for (stim=1:length(stimf))     
        trials=stimf(clustval,stim).cnt;
        ind2=[];
        for trialnum=1:trials
            cr_arry=stimf(clustval,stim).spkarray;
            
                ind2=find(cr_arry(:,2)==trialnum);
                stimf(clustval,stim).histdist(trialnum,:)=(histc(cr_arry(ind2,1),edges))
        %histdist3{stim}(trialnum,:)=(histc(spkind(spksin{clustnum}(ind(ind2)),1),edges2))
            
        end
        ax1 = gca;
        stimf(clustval,stim).meanhist=mean(stimf(clustval,stim).histdist)
    end    
end    
    
save stimf.mat stimf  

        
   
        
        
    
        
   