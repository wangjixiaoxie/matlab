function [outvars]=analyzeHVCallspikes_y32g92(neuron,clusters,wins,harmwin,model)


%harmwin=[50 65;50 65;110 120;60 80];


for mls=1:length(wins);
    winstart=wins(mls);
    winstop=wins(mls)+50;

    for numsylls=1:size(harmwin,1)

        % Analyze - this loop returns:
        % spikecount{unit #}(rendition)   --- number of spikes in premotor window of the harmonic stack
        % pitchmean{unit #}(rendition) --- mean pitch in the harmonic stack
        % for the first syllable, the harmonic part goes from 45-70ms (includes 16ms adjustment - so really 29-54ms)
        pitchtimes=[1/8:1/8:1745/8]+16; % 16ms because of window size
        harmstart=harmwin(numsylls,1);
        harmend=harmwin(numsylls,2);
        neuroncount=0;
        clear spikecount pitchmean
        for thisfolder=1:length(clusters) % for the neuron/recording wire
            goodclust=clusters{thisfolder};
            for thissyllable=numsylls % for the syllable
                for thisclust=1:length(goodclust)     % for the cluster
                    neuroncount=neuroncount+1; % each unit
                    count=0; % each rendition with this unit
                    for thissong=1:length(neuron(thisfolder).song)
                        for thisrendition=1:length(neuron(thisfolder).song(thissong).noteons{thissyllable})
                            timewindow=harmstart+winstart:harmend+winstop;
                            spikevals=neuron(thisfolder).song(thissong).spiketimes{thisclust}...
                                -neuron(thisfolder).song(thissong).noteons{thissyllable}(thisrendition);
                            spikecounter=length(find(spikevals>harmstart+winstart & spikevals<harmend+winstop));
                            if ~isempty(spikecounter)
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

        if model==1
            % Only analyze results for the neurons that actually spike in the vicinity
            % of the harmonic stack note - the variable 'doesspike' evaluates this
            clear doesspike
            for neuralunit=1:length(spikecount) % for each 'neuron'
                doesspike(neuralunit)=sum(spikecount{neuralunit}())>0; % nine units
            end

            % Calculate correlation coefficients and corresponding t-stats
            clear rval tval
            for neuralunit=1:length(spikecount) % for each 'neuron'
                if doesspike(neuralunit)
                    m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                    r=m(2);
                    rval(neuralunit)=r;
                    tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                end
            end
        else
            if model==2;
                clear doesspike
                for neuralunit=1:length(spikecount) % for each 'neuron'
                    doesspike(neuralunit)=sum(spikecount{neuralunit}())>0; % nine units
                end

                % Calculate correlation coefficients and corresponding t-stats
                clear rval tval
                for neuralunit=1:length(spikecount) % for each 'neuron'
                    if doesspike(neuralunit)
                        m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                        r=m(2);
                        rval(neuralunit)=r;
                        tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                    end
                end
            else
                if model==3;
                    % Only analyze results for the neurons that actually spike in the vicinity
                    % of the harmonic stack note - the variable 'doesspike' evaluates this
                    clear doesspiketwice
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        doesspiketwice(neuralunit)=length(spikecount{neuralunit}())>2;
                    end

                    % Calculate correlation coefficients and corresponding t-stats
                    clear rval tval
                    for neuralunit=1:length(spikecount) % for each 'neuron'
                        if doesspiketwice(neuralunit)
                            m=corrcoef(pitchmean{neuralunit}(),spikecount{neuralunit}());
                            r=m(2);
                            rval(neuralunit)=r;
                            tval(neuralunit)=r/sqrt((1-r^2)/(length(spikecount{neuralunit}()-2)));
                        end
                    end

                end
            end
        end
        rval(doesspiketwice);
        tval(doesspiketwice);
        Syll.rval=rval;
        Syll.tval=tval;
        Syll.spikecount=spikecount;
        Syll.pitchmean=pitchmean;
        outvars.mn(numsylls,mls)=median((Syll.tval(Syll.tval~=0 & ~isnan(Syll.tval))));
        outvars.se(numsylls,mls)=std((Syll.tval(Syll.tval~=0 & ~isnan(Syll.tval))))/sqrt(sum(Syll.tval~=0 & ~isnan(Syll.tval)));
        outvars.neuroncount(numsylls,mls)=sum(Syll.tval~=0 & ~isnan(Syll.tval));
    end
end


                
         