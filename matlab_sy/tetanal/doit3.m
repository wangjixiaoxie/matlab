%What I would like to do is add more points to existing cluster with two
%different colors.  One color for all previous points, another color for
%all new points.

%And that when I change the cluster dimension.  It does so for all the
%points.



batchfile='/cobain/g79g80/redtet/HVC/batch5'
outfilefilter={'D:\g79g80\g79g80_bos.wav' 'D:\g79g80\g79g80_revbos.wav' 'D:\g79g80\g79g80_p10bos.wav' 'D:\g79g80\g79g80_m10bos.wav' 'D:\g79g80\g79g80_p5.wav' 'D:\g79g80\g79g80_m5.wav'} 
structnames={'bos' 'revbos' 'p10' 'm10' 'p5' 'm5'}
%comb_cluster=5;
display_skip=5;
secondflag=0;
clear clust_hand;
ampbuffer=zeros(1,5)

%initialize spike and amp structures to zero%
for ii=1:size(outfilefilter) 
    strname=structnames{ii};
    eval(['spkamp.' strname '=zeros(1,5);']);
    eval(['spkind.' strname '=zeros(1,2);']);
    eval(['spkampclust.' strname '=zeros(1,5);']);
    eval(['spkindclust.' strname '=zeros(1,2);']);

end

clear fnm;

%if (secondflag==0)
  %  clear indvec
   % clear ampvec
   %clear spkind
    %clear spkamp
%end


%get all the files in the batch file
cnt=0;
badfilecount=1;
anfilecount=1;
filtfilecount=1;
clustnum=1;

fid=fopen(batchfile,'r');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	if (~exist(fn,'file'))
		continue;
	end
	
	cnt=cnt+1;
	fnm{cnt}=fn;
end
fclose(fid);

% this goes through and gets the main spike cluster for each file in the
%batch file
tet_chans=[3:6];song_chan=1;

%chans=0;
%x=0;
%y=0;

%This loops through all the files in the batch file
for ifn=1:length(fnm)
    
    %Every display_skip number of trials, 
    if(mod(filtfilecount,display_skip)==0)
        filtfilecount=filtfilecount+1;
        skiptrial=0;
    else
        skiptrial=1;
        filtfilecount=filtfilecount+1;
    end
    
    fn=fnm{ifn};
    %run tetanal on the file get spike times and spkamps
	
    rd=readrecf(fn);
    %read in data
    [data,fs,spkindnew,spkampnew]=tetanal(fn,-2000,song_chan,tet_chans);
    
    spkindnew(:,2)=ifn;
    spkampnew(:,5)=ifn;
    
    %Add to buffer of <display_skip> trials worth of data>
    ampbuffer=[ampbuffer;spkampnew];    
    
    %This loop assigns data to appropriate data structure
    for ii=1:size(outfilefilter) 
        if (strcmp(rd.outfile,outfilefilter{ii}))
            strname=structnames{ii};
            eval(['spkamp.' strname '=[spkamp.' strname ';spkampnew];']);
            eval(['spkind.' strname '=[spkind.' strname ';spkindnew];']);
            break;
        end
    end

    if(skiptrial==0)
       pltscat2;
    
    
        while(1)
            R=input('Keep Cluster? 1 for Yes')
            %strname is the correct name for clustered spikes%
            if (R==1)
                eval(['spkampclust.' strname '=[spkampclust.' strname ';inspikes];']);
                eval(['spkindclust.' strname '=[spkindclust.' strname ';indbuffer(ppind)];']);
                break;
            else
                [x,y]=ginput();
                pltscat2;
                clustnum=clustnum+1;
            end
        end
    
        
    end
    %clear the buffer and initialize it to 0%
    if(skiptrial==0)
        clear ampbuffer;
        ampbuffer=zeros(1,5);
        clustnum=clustnum+1;
    end
end


        %pltdat;pltscat;
    
    
        if(secondflag)
            insd=inpolygon(spkamp(:,chans(1)),spkamp(:,chans(2)),x,y);
            pp=find(insd==1);
        end       
    
 %{ 
spkind=indvec;
    spkamp=ampvec;


    pltdat;pltscat;
    
    %get chans, x, and y from this .mat file
	%load CLUSTBND
    if(ifn>1)
        insd=inpolygon(spkamp(:,chans(1)),spkamp(:,chans(2)),x,y);
        pp=find(insd==1);

    %this finds the spike times inside the 
    %cluster
	spkind_clust = spkind(pp)/fs; %in seconds
	
    rd=readrecf(fn);
    
    %save this info 
	eval(['save ',fn,'.spkind spkind_clust rd']);
	disp(fn);

%}	
 