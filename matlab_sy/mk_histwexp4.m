%mk_histwexp4.m
%revised to deal with data saved in sam's .mat structure.


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
binsize=.020
plotcount=1;
plotsum=1;
fs_sound=44100;
spkchn=2;
rawsong=[];

%clear stim

stimleng=length(dat(:,1))/32000;

stimf=[];
%load '/cobain2/twarren2/w90o24/stim/stimexp2.mat' corpshiftednormg

%load '/cobain4/twarren4/g100o55/stim/stim.mat' xlist corpshiftednormg

%this needs to be updated.
% rawwav='/home/twarren/pu3bu86tw.wav';
% if(~exist('rawsong'))
%         [rawsong,fs_sound]=wavread(rawwav);
% end




%TO TEST THIS I NEED TO FIX STIMNAMES!!!, so that it is in above format.

 stimnames={'pu3bu86tw.wav' 'pu3bu86revtw.wav' 'pu3bu86m10tw.wav' 'pu3bu86p10tw.wav'}
 %stimnames={'D:\tim\g100o55e2\stim\1055.wav' 'D:\tim\g100o55e2\stim\1055rev.wav' 'D:\tim\g100o55\stim\1055m10.wav' 'D:\tim\g100o55\stim\1055p10.wav'  'D:\tim\g100o55\stim\1055jig.wav' 'D:\tim\g100o55\stim\10055jigcomp.wav'}
%  structnames={'bs4' 'rev' 'p4 ' 'm4 ' 'bs2' 'm2 ' 'p2 ' 'm8 ' 'bs8' 'p8 ' }
%structnames={'bos' 'dd' 'du'  'ud' 'uu' 'nd' 'nu'}; 
bt='batch';
% stimleng=length(rawsong)/fs_sound;%in seconds

%calculate stimleng

binwid=1e-3;
%stimleng=ceil(stimleng/binwid); % # of bins

histbinwid=5.0e-3;
stimf=[];fnames=[];cnt=0;
stimf.outind=[];


fid=fopen(bt,'r');
count=1;

%this loop just initializes

    for ii=1:length(stimnames)
        stimf(ii).cnt=0;;
        stimf(ii).spkarray=[];
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
    cmd=['load -mat ',fn,'.spk.mat'];
        eval(cmd);

	cnt=cnt+1;
	fnames(cnt).fn=fn;
    
    %go through all the stimnames, and assign the stimname
    for ii = 1:length(stimnames)
		stimf(ii).outfile=stimnames{ii};
        
        %This should only be the case
        if (strcmp(pbfile,stimnames{ii}))
			ind=ii;
			 stimtype(lpcnt)=ii;
            stimf(ind).cnt=stimf(ind).cnt+1;
            stimf(ind).outind=[stimf(ind).outind cnt]; 
            break;
        
        else
            stimtype(lpcnt)=0;
        end
        end
   
    lpcnt=lpcnt+1;
end
fclose(fid);

% do for each cluster bound


%critical loop which takes the clustered data and puts into th_e stimf
%structure, which is used for later analysis  (rasters, etc...)

	%clear out the rasters
	spkt=[]
    
    for kk = 1:length(fnames)
		rd=readrecf(fnames(kk).fn);
		%this is the really critical line which reads the .clust file
        cmd=['load -mat ',fnames(kk).fn,'.spk.mat'];
        eval(cmd);

        stvl=stimtype(kk);
        if(stvl)
            numstm=find(stimf(stvl).outind==kk);
            spkt=spktimes;
            if(isempty(spkt))
                spkt=-10;
            end
        
        trials=numstm;
        
        stimf(stvl).rast=spkt;
        stimf(stvl).spkarray{trials}=[spkt];
        end
    end

   
    edges=0:binsize:stimleng;
        %edges2=[0 startsong*32000 endsong*32000]
    for (stim=1:length(stimf))     
        trials=stimf(stim).cnt;
        ind2=[];
        for trialnum=1:trials
            cr_arry=stimf(stim).spkarray{trialnum};
            
                
                stimf(stim).histdist(trialnum,:)=(histc(cr_arry,edges))
        %histdist3{stim}(trialnum,:)=(histc(spkind(spksin{clustnum}(ind(ind2)),1),edges2))
            
        end
        ax1 = gca;
        stimf(stim).meanhist=mean(stimf(stim).histdist)
    end    
   
    
save stimf.mat stimf  

        
   
        
        
    
        
   