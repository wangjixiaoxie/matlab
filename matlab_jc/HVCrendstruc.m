function [Rendition]=HVCrendstruc(neuron,ptwindows,theseclusters,numsylls,wins)
%
count=0;
     for thisfolder=1:length(neuron)  
     for thissyllable=1:numsylls
         FFmeans(thisfolder,thissyllable).data=[];
         % pitch window measurement time relative to onset of this rendition, milliseconds
            adjfactor=16+mean(ptwindows(thissyllable).data)/8; 
     for thissong=1:length(neuron(thisfolder).song)
     if length(neuron(thisfolder).song(thissong).pitchdata)>thissyllable-1
         % spike times for this song
             thisclust=theseclusters{thisfolder};
             storespikes=[];
             for i=1:length(thisclust)
                 cclust=thisclust(i);
                 if cclust<length(neuron(thisfolder).song(thissong).spiketimes)+1
                    storespikes=[storespikes neuron(thisfolder).song(thissong).spiketimes{cclust}];
                 end
             end
     if length(neuron(thisfolder).song(thissong).pitchdata)>thissyllable-1
     for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
         count=count+1;
         % Basic information
             % .neuron
               Rendition(count).neuron=thisfolder;          
             % .sylltype
               Rendition(count).sylltype=thissyllable;
             % .song    % which song
               Rendition(count).song=thissong;
             % .rendition  find(Rendition(count).spiketimes>wins(j) % which rendition in this song
               Rendition(count).rendition=thisrendition;         
         % Pitch information
             % .pitchcontour
               Rendition(count).pitchcontour=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
             % .ptwindow
               Rendition(count).ptwindow=ptwindows(thissyllable).data;
             % .pitchmean
               Rendition(count).pitchmean=mean(Rendition(count).pitchcontour(Rendition(count).ptwindow));
               FFmeans(thisfolder,thissyllable).data=[FFmeans(thisfolder,thissyllable).data Rendition(count).pitchmean];
         % Context information
               thistime=neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition)+adjfactor;
                clear sylltimes
                for i=1:length(neuron(thisfolder).song(thissong).noteons) % for each syllable
                    % .relativesylltimes       % onsets relative to onset of this rendition, milliseconds 
                    Rendition(count).relativesylltimes{i}=neuron(thisfolder).song(thissong).noteons{i}()-thistime;
                    mns(i)=min(neuron(thisfolder).song(thissong).noteons{i}()); % earliest onset
                end
            % .sylltime        % pitch window measurement time relative to onset of first syll in song, milliseconds
               Rendition(count).sylltime=thistime-min(mns);
         % Spikes information
               Rendition(count).spiketimes=storespikes-thistime; % spikes relative to pitch measurement window
               for j=1:length(wins)-1
                   Rendition(count).fr_rate(j)=length(find(Rendition(count).spiketimes>wins(j) & Rendition(count).spiketimes<wins(j+1)));
                   Rendition(count).fr_time(j)=wins(j)+0.5;
               end
               
     end
     end
     end
     end
     end
     end
     
     % Pitch residuals information (calculated separately for separate neurons)
         % avg pitch for each syllable/neuron
             for thisfolder=1:length(neuron)
                 for thissyllable=1:numsylls
                     Fmeans(thisfolder,thissyllable)=mean(FFmeans(thisfolder,thissyllable).data);
                 end
             end
         % calculate pitch residuals for each
             for i=1:length(Rendition)
                 Rendition(i).pitchresiduals=Rendition(i).pitchmean-Fmeans(Rendition(i).neuron,Rendition(i).sylltype);
             end
    % Bursts information
            for i=1:length(neuron)
                for j=1:numsylls
                    allrates(i,j).data=[];
                end
            end
            for i=1:length(Rendition)
                allrates(Rendition(i).neuron,Rendition(i).sylltype).data=[allrates(Rendition(i).neuron,Rendition(i).sylltype).data;Rendition(i).fr_rate];
            end
    % Define bursts based on 
        ravgwin=5;
        tallpeakheight=10; % burst peak firing rate must be 10* mean firing rate
        smallpeakwidth=2; % burst min firing rate must be 2* mean firing rate
        burstparameters=[];
            for i=1:length(neuron)
            for j=1:numsylls
            % Find burst centers
            tallpeaks=find(runningaverage(mean(allrates(i,j).data),ravgwin)>tallpeakheight*mean(mean(allrates(i,j).data)));
            smallpeaks=find(runningaverage(mean(allrates(i,j).data),ravgwin)>smallpeakwidth*mean(mean(allrates(i,j).data)));
            if length(smallpeaks)>0
                smalldiffs=[smallpeaks(1) diff(smallpeaks) ];
            else
                smalldiffs=[];
            end
            countstart1=0;countend1=0;burststartcandidate=[];burstendcandidate=[];
            for k=2:length(smalldiffs)-1
                if smalldiffs(k)==1 & smalldiffs(k-1)>1
                    countstart1=countstart1+1;
                    burststartcandidate(countstart1)=smallpeaks(k)+floor(ravgwin/2);
                end
                if smalldiffs(k)==1 & (smalldiffs(k+1)>1 | k+1==length(smalldiffs))
                    countend1=countend1+1;
                    burstendcandidate(countend1)=smallpeaks(k)+floor(ravgwin/2);
                end
            end
            counters=0;
            for k=1:length(burststartcandidate)
                starting=burststartcandidate(k);
                a=burstendcandidate(find(burstendcandidate>burststartcandidate(k)));
                if ~isempty(a)
                    ending=a(1);
                    if ~isempty(find(tallpeaks>starting & tallpeaks<ending))
                        counters=counters+1;
                        burstparameters(i,j).times(1,counters)=wins(burststartcandidate(k))-1; % convert back to spiketimes
                        burstparameters(i,j).times(2,counters)=wins(a(1))+1;
                    end
                end
            end
            end
            end
      % Save burst data
      for i=1:length(Rendition)
          if Rendition(i).neuron>size(burstparameters,1) | Rendition(i).sylltype>size(burstparameters,2)
                Rendition(i).burstparams=[];
          else
          Rendition(i).burstparams=[burstparameters(Rendition(i).neuron,Rendition(i).sylltype).times];
          for j=1:size(burstparameters(Rendition(i).neuron,Rendition(i).sylltype).times,2) % for each burst
              Rendition(i).burstspiketimes{j}=Rendition(i).spiketimes(find(Rendition(i).spiketimes>burstparameters(Rendition(i).neuron,Rendition(i).sylltype).times(1,j)...
                  & Rendition(i).spiketimes<burstparameters(Rendition(i).neuron,Rendition(i).sylltype).times(2,j)));
              Rendition(i).burstspikecounts(j)=length(Rendition(i).burstspiketimes{j});
              if isempty(Rendition(i).burstspiketimes{j})
                  Rendition(i).mnburstspiketime(j)=-1000;
              else
                Rendition(i).mnburstspiketime(j)=mean(Rendition(i).burstspiketimes{j})-burstparameters(Rendition(i).neuron,Rendition(i).sylltype).times(1,j);
              end
          end
          end
      end
      % Get burst info as residuals
      counters=zeros(length(neuron),numsylls);
      Allresidinfo=[];
      for i=1:length(Rendition)
          counters(Rendition(i).neuron,Rendition(i).sylltype)=counters(Rendition(i).neuron,Rendition(i).sylltype)+1;
          if ~isempty(Rendition(i).burstspikecounts)
              Allresidinfo(Rendition(i).neuron,Rendition(i).sylltype).data(counters(Rendition(i).neuron,Rendition(i).sylltype),:)=...
                  Rendition(i).burstspikecounts;
              AllresidinfoTIMING(Rendition(i).neuron,Rendition(i).sylltype).data(counters(Rendition(i).neuron,Rendition(i).sylltype),:)=...
                  Rendition(i).mnburstspiketime;
          end
      end
      for i=1:length(Rendition)
          if ~isempty(Rendition(i).burstspikecounts)
              Rendition(i).mnburstspikecounts=mean(Allresidinfo(Rendition(i).neuron,Rendition(i).sylltype).data);
              Rendition(i).residburstspikecountsDIFF=Rendition(i).burstspikecounts-Rendition(i).mnburstspikecounts;
              Rendition(i).residburstspikecountsFACTOR=Rendition(i).burstspikecounts./Rendition(i).mnburstspikecounts; 
              for j=1:length(Rendition(i).mnburstspiketime)
                  if Rendition(i).mnburstspiketime(j)==-1000
                      Rendition(i).residburstspiketime(j)=-1000;
                  else
                      Rendition(i).residburstspiketime(j)=Rendition(i).mnburstspiketime(j)-...
                          mean(AllresidinfoTIMING(Rendition(i).neuron,Rendition(i).sylltype).data(...
                          find(AllresidinfoTIMING(Rendition(i).neuron,Rendition(i).sylltype).data(:,j)~=-1000),j));
                  end
              end
          end
      end
      
   % MOTIF NUMBER
   for i=1:length(Rendition)
       if Rendition(i).sylltime<340
           Rendition(i).motif=0;
       else
           if Rendition(i).sylltime<650
               Rendition(i).motif=1;
           else if Rendition(i).sylltime<1260
                   Rendition(i).motif=2;
               else if Rendition(i).sylltime<1850
                       Rendition(i).motif=3;
                   else
                       Rendition(i).motif=100;
                   end
               end
           end
       end
   end