%What I would like to do is add more points to existing cluster with two
%different colors.  One color for all previous points, another color for
%all new points.

%And that when I change the cluster dimension.  It does so for all the
%points.


%modified from doit5 to cluster n clusters at once.

%initialize variables
%batchfile='//cobain2/exp2b3b4/bluetet/batch8'
%outfilefilter={'D:\tim\b3b4\exp2stimuli\34.wav' 'D:\tim\b3b4\exp2stimuli\34rev.wav' 'D:\tim\b3b4\exp2stimuli\34l.wav' 'D:\tim\b3b4\exp2stimuli\34m10.wav' 'D:\tim\b3b4\exp2stimuli\34p10.wav' 'D:\tim\b3b4\exp2stimuli\34m5.wav' 'D:\tim\b3b4\exp2stimuli\34p5.wav' 'D:\tim\b3b4\exp2stimuli\34m10tl.wav' 'D:\tim\b3b4\exp2stimuli\34p10tl.wav' 'D:\tim\b3b4\exp2stimuli\34q.wav'}
colorlist={'r' 'g' 'b'} 
structnames={'bos' 'rev' 'lou' 'm10' 'p10' 'm5 ' 'p5 '  'm10t' 'p10t' 'qui'}


display_skip=5;
secondflag=0;
clust_hand={};
ampbuffer=[];
indbuffer=[];
%spksinclust={};
spksoutclust=[];
spkindclust=[];
spkamp=[];
spkind=[];
inspikes={};
ppin=[];
insd={};
cnt=0;
offset=0;
totclusts=1;
badfilecount=1;
anfilecount=1;
filtfilecount=1;
clustnum=1;
clear fnm;
tet_chans=[2:5];song_chan=1;
totalpolygons=0;

%get all the files in the batch file
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

%This loops through all the files in the batch file
for ifn=1:length(fnm)
    
    %Every display_skip number of trials, skip a trial(skiptrial=1)
    if(mod(filtfilecount,display_skip)==0||ifn==length(ifn))
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
    if(isempty(rd)==0)
        [data,fs,spkindnew,spkampnew]=tetanaldc(fn,-2000,song_chan,tet_chans);
    
    %This loop assigns data to appropriate data structure
        inputnum=0;
    
        for ii=1:length(outfilefilter) 
            if (strcmp(rd.outfile,outfilefilter{ii}))
                inputnum=ii;
                break;
            end
        end

    %Next-to-last column is trial number, last column is input type
        spkindnew(:,2)=ifn;
        spkampnew(:,5)=ifn;
        spkindnew(:,3)=inputnum;
        spkampnew(:,6)=inputnum;
    
    
    %pltdat;
    
    %R=input('Keep data? 1 for Yes')
    
    
    %Add to buffer of <display_skip> trials worth of data>
    
    %if (R==1)
        ampbuffer=[ampbuffer;spkampnew];    
        indbuffer=[indbuffer;spkindnew];
    %end
        if(skiptrial==0|ifn==length(fnm))
        if (skiptrial==0)
                pltscat6;
        else
                pltscat6;
        end
        if(ifn>1)
            pltdat;
            figure(DATPLT);
                    %for ii = 1:Nusechans
                       
                      %  
                       % subplot(Nusechans+1,1,ii+1);
                        %hold on;grid on;
                        %plot(spkt(pp),dat(spki,usechans(ii-1)),[clr,'o']);
                        %plot(indbuffer(ppin(indices2),1)/fs,-inspikes(indices,ii),[clr,sym]);
                    %end
       
       
            figure(CLUSTPLT);
        end
        clustnum=clustnum+1;
    
        while(1)
            %first time through
            if(clustnum==2)
                totclusts=input('How many clusters?')
                clustnum=clustnum+1;
                for jj=1:totclusts
                    spksinclust{jj}=[];
                end
            end
            R=input('Keep Cluster? 1 for Yes')
            %strname is the correct name for clustered spikes%
            if (R==1)
               for jj=1:totclusts
                    spksinclust{jj}=[spksinclust{jj};ppin{jj}+offset];
  %                  spksoutclust{jj}=[spksoutclust{jj};ppout{jj}+offset];
                    totalpolygons=totalpolygons+1;
                    clusterinfo{jj,totalpolygons}.dims=[x{jj},y{jj}];
                    clusterinfo{jj,totalpolygons}.chans=[chans{jj}];
                    
               end
               break; 
            else
                figure(CLUSTPLT);
                for jj=1:totclusts
                    [x{jj},y{jj}]=ginput();
                    clust_hand{jj}=gca;
                    
                end
                pltscat6;
                pltdat;
                figure(DATPLT);
                for ii = 2:Nusechans+1
                    for jj=1:totclusts
                    
                        ax(ii)=subplot(Nusechans+1,1,ii);
                        hold on;grid on;
                        %plot(spkt(pp),dat(spki,usechans(ii-1)),[clr,'o']);
                        plot(0:1/32000:((length(data)-1)/fs),data(:,ii))
                        %indices=find(spkind(spksinclust{jj},2)==ifn);
                        indices=find(inspikes{jj}(:,5)==ifn);
                        indices2=find(indbuffer(ppin{jj},2)==ifn);
                        xrow=[indbuffer(ppin{jj}(indices2),1)/fs indbuffer(ppin{jj}(indices2),1)/fs];
                        yrow=[ones(length(xrow),1)-5000 ones(length(xrow),1)-7000];
                        plot(xrow',yrow','LineWidth',2,'Color',colorlist{jj})
                    end    
                end
            end
                
   % if(indices)
    %plot(spkind(spksin{1}(indices,1))/fs,-7000,'r.');
    %end
    
    %indices=find(spkind(spksin{3},2)==75);
    %plot(spkind(spksin{3}(indices,1))/fs,-5000,'g.');  
        end
        
        spkind=[spkind; indbuffer];
        spkamp=[spkamp;ampbuffer];
        offset=offset+length(ampbuffer(:,1));
        ampbuffer=[];
        indbuffer=[];
        clustnum=clustnum+1;
 
  
        end
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
 