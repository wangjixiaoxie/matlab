function [outvars,outdata]=analyzeHVCallspikes2(neuron,clusters,wins,harmwin,indfolders,model)
neuron=neuron(indfolders);
for k=1:length(indfolders)
    clustersnew{k}=clusters{indfolders(k)}();
end


for mls=1:length(wins) % for each perimotor window
    winstart=wins(mls);
    winstop=wins(mls)+50;
mls
    for numsylls=1:size(harmwin,1) % for each syllable
        % Analyze - this loop returns:
        % spikecount{unit #}(rendition)   --- number of spikes in premotor window of the harmonic stack
        % pitchmean{unit #}(rendition) --- mean pitch in the harmonic stack
        % for the first syllable, the harmonic part goes from 45-70ms (includes 16ms adjustment - so really 29-54ms)
        pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
        harmstart=harmwin(numsylls,1);
        harmend=harmwin(numsylls,2);
        neuroncount=0;
        clear spikecount pitchmean % clear things that get incremented in the loop
        for thisfolder=1:length(clustersnew) % for the neuron/recording wire
            goodclust=clustersnew{thisfolder};
            thissyllable=numsylls; % for each syllable
            for thisclust=goodclust     % for the cluster
                neuroncount=neuroncount+1; % each unit
                count=0; % each rendition with this unit
                for thissong=1:length(neuron(thisfolder).song)
                    if length(neuron(thisfolder).song(thissong).noteons)>thissyllable-1 % almost always true (not sure why it sometimes isn't)
                        for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                            timewindow=harmstart+winstart:harmend+winstop;
                            spikevals=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                            spikecounter=length(find(spikevals>harmstart+winstart & spikevals<harmend+winstop));
                            if model==1
                                if ~isempty(spikecounter)
                                    count=count+1;
                                    spikecount{neuroncount}(count)=spikecounter;
                                    pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                    pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                end
                            else if model==2
                                    if ~isempty(spikecounter)
                                        count=count+1;
                                        if spikecounter==0
                                            spikecount{neuroncount}(count)=0;
                                        else
                                            spikecount{neuroncount}(count)=1;
                                        end
                                        pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                        pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                    end
                            else if model==3
                                    if ~isempty(spikecounter) & spikecounter>0
                                        count=count+1;
                                        spikecount{neuroncount}(count)=spikecounter;
                                        pitchvals=neuron(thisfolder).song(thissong).pitchdata{thissyllable}.datt(thisrendition,:);
                                        pitchmean{neuroncount}(count)=mean(pitchvals(find(pitchtimes>harmstart & pitchtimes<harmend)));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        end
        if model==1
            % Only analyze results for the neurons that actually spike in the vicinity
            % of the harmonic stack note - the variable 'doesspike' evaluates this
            clear doesspike
            for neuralunit=1:length(spikecount) % for each 'neuron'
                doesspike(neuralunit)=sum(spikecount{neuralunit}())>0; % nine units
            end

            % Calculate correlation coefficients and corresponding t-stats
            rval=[];tval=[];
            for neuralunit=1:length(spikecount) % for each 'neuron'
                if doesspike(neuralunit)
                    m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                    r=m(2);
                    rval(neuralunit)=r;
                            if r>0.999 % This is to avoid huge tvals
                                tval(neuralunit)=3;
                            else
                                if r<-0.999
                                    tval(neuralunit)=-3;
                                else
                                    tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                                end
                            end
                end
            end
        else
            if model==2;
                clear doesspike
                for neuralunit=1:length(spikecount) % for each 'neuron'
                    doesspike(neuralunit)=sum(spikecount{neuralunit}())>0; % nine units
                end

                % Calculate correlation coefficients and corresponding t-stats
            rval=[];tval=[];
                for neuralunit=1:length(spikecount) % for each 'neuron'
                    if doesspike(neuralunit)
                        m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                        r=m(2);
                        rval(neuralunit)=r;
                            if r>0.999
                                tval(neuralunit)=3;
                            else
                                if r<-0.999
                                    tval(neuralunit)=-3;
                                else
                                    tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                                end
                            end
                    end
                end
            else
                if model==3;
                    % Only analyze results for the neurons that actually spike in the vicinity
                    % of the harmonic stack note - the variable 'doesspike' evaluates this
                    clear doesspike
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        doesspike(neuralunit)=length(spikecount{neuralunit}())>2;
                    end

                    % Calculate correlation coefficients and corresponding t-stats
            rval=[];tval=[];
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        if doesspike(neuralunit)
                            m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                            r=m(2);
                            rval(neuralunit)=r;
                            if r>0.999
                                tval(neuralunit)=3;
                            else
                                if r<-0.999
                                    tval(neuralunit)=-3;
                                else
                                    tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                                end
                            end
                        end
                    end

                end
            end
        end
        rval(doesspike);
        tval(doesspike);
        Syll.rval=rval;
        Syll.tval=tval;
        Syll.spikecount=spikecount;
        Syll.pitchmean=pitchmean;
        outvars.Tmn(numsylls,mls)=mean((Syll.tval(Syll.tval~=0 & ~isnan(Syll.tval))));
        outvars.Tse(numsylls,mls)=std((Syll.tval(Syll.tval~=0 & ~isnan(Syll.tval))))/sqrt(sum(Syll.tval~=0 & ~isnan(Syll.tval)));
        outvars.Rmn(numsylls,mls)=mean((Syll.rval(Syll.rval~=0 & ~isnan(Syll.rval))));
        outvars.Rse(numsylls,mls)=std((Syll.rval(Syll.rval~=0 & ~isnan(Syll.rval))))/sqrt(sum(Syll.rval~=0 & ~isnan(Syll.rval)));        
        outvars.neuroncountT(numsylls,mls)=sum(Syll.tval~=0 & ~isnan(Syll.tval));
        outvars.neuroncountR(numsylls,mls)=sum(Syll.rval~=0 & ~isnan(Syll.rval));        
        outdata(numsylls,mls).tvaldata=tval;
        outdata(numsylls,mls).rvaldata=rval;
    end
end


                
         