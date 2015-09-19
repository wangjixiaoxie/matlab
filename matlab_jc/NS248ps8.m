% PROBLEM 1 - Code to simulate a spike train given a rate function (i.e. model)
        function spiketimes=simulatespiketrain(times,rates)
            % Generate ISIs

            i=0;
            index=1;
            % until you reach the end of your time window (t>tmax)
            while index < length(times)
                i=i+1;
                r=exprnd(1);
                % Given your current time (t), how long must you wait for the next spike?
                % Calculate this by integrating the firing rate function from the
                % time of the previous spike until you get to an area > r
                j=0;
                % Keep integrating until you predict a spike or until you pass the
                % edges of the time domain of your rate function
                areaundercurve=0;
                while areaundercurve < r && index+j < length(times)
                    j=j+1;
                    areaundercurve=trapz(times(index:index+j),rates(index:index+j));
                end
                % This is the estimated time of the next spike
                index=index+j;
                spiketimes(i)=time1(index);
            end

% PROBLEM 2 - Simulate a spike train from your model using code in problem 1
    % Do 100 simulations 
        rate1=[];
        time1=[0:0.001:199.999];
        for i=1:100
            rate1=[rate1 rate];
        end
    % Calculate
        spiketimes1=simulatespiketrain(time1,rate1);
    % Correct all the times for the 100 simulations so that you can combine them into a PSTH
                % - remainder after division by 2 seconds
        spiketimes1=rem(spiketimes1,2);
    % Make PSTH
        starting=0;
        stopping=1.9999;
        binsize=0.010;
        % for each bin
        for i=1:200
            binstart=binsize*(i-1);
            binend=binsize*i;
            PSTH(i)=sum(spiketimes1>binstart & spiketimes1<binend);
        end
        figure;plot([0:0.01:1.999],PSTH)
        hold on;plot(time,rate,'r') 


% PROBLEM 3 - Calculate how well the original model fits the data generated
% in problem 2.  Also calculate how well the PSTH generated in problem 2
% fits that data.

% Function that does the time-rescaling theorem and determines z's that
% can be used to test for goodness of fit.
        function z=computeksstat(times,spiketimes,rates,binwidth)
            lastspiketime=0;
            % Integrate between adjacent spikes
            for i=1:length(spiketimes)-1
                % Which part of the firing rate function is between the last spike and
                % the next spike?
                relevant_rates=rates(times >= spiketimes(i) & times < spiketimes(i+1));
                % Integrate over the relevant part of the firing rate:
                z(i)=trapz(relevant_rates)*binwidth;
            end

% See how well the original model does
        newtrials=[1 find(diff(spiketimes1)<0)+1 length(spiketimes1)];
        % for each trial
        z=[];
        binwidth=0.001;
        for i=1:100
            spikes=spiketimes1(newtrials(i):newtrials(i+1)-1);
            z=[z computeksstat(time,spikes,rate,binwidth)];
        end
        znew=1-exp(-z);
        
        % Test 1 - should equal one (approximately)
            mean(z)=0.9542;
        % Test 2 -  should be uniform (approximately)
            figure;hist(znew) 
        % Test 3 - should equal zero (approximately)
            corrcoef(z(1:end-1),z(2:end))=0.0133
        % Test 4 - ksplot should be uniform
            figure;getKSPlot(znew)

% See how well the PSTH does.
        % for each trial
        time2=[0:0.01:1.99];
        z2=[];
        binwidth=0.01;
        for i=1:100
            spikes=spiketimes1(newtrials(i):newtrials(i+1)-1);
            z2=[z2 computeksstat(time2,spikes,PSTH,binwidth)];
        end
        z2new=1-exp(-z2);
        % Test 1
            mean(z2)=0.6934;   % less than a perfect fit
        % Test 2
            figure;hist(z2new) % way too many are around zero 
                 % This means that the model has too many small ISIs, possibly
                 % because it fails to take into account the refractory
                 % period.
        % Test 3
            corrcoef(z2(1:end-1),z2(2:end))=0.0176
        % Test 4
            figure;getKSPlot(z2new)
% PROBLEM 4 - write code that simulates a spike train by taking a model
% (i.e. rate function - PSTH from problem 2 in this case) and 
% incorporating a refractory period by using -1*an exponential decay term
        function spiketimes=simulatespiketrainEXP(times,rates)
            % Generate ISIs
            i=0;
            index=1;
            % until you reach the end of your time window (t>tmax)
            while index < length(times)
                i=i+1;
                r=exprnd(1);
                % Given your current time (t), how long must you wait for the next spike?
                % Calculate this by integrating the firing rate function from the
                % time of the previous spike until you get to an area > r
                j=0;
                % Keep integrating until you predict a spike or until you pass the
                % edges of the time domain of your rate function
                areaundercurve=0;
                while areaundercurve < r && index+j < length(times)
                    j=j+1;
                    refractoryperiodfactor(index:index+j)=exp(-0.002./(times(index:index+j)-times(index)));
                    areaundercurve=trapz(times(index:index+j),rates(index:index+j).*refractoryperiodfactor(index:index+j));
                end
                % This is the estimated time of the next spike
                index=index+j;
                spiketimes(i)=times(index);
            end
% Calculate the PSTH generated by this exp simulation (as above)
        % Simulate
                rate2=[];
                time2=[0:0.01:199.99];
                for i=1:100
                    rate2=[rate2 PSTH];  % bincount IS the PSTH 
                end
            spiketimes2=simulatespiketrainEXP(time2,rate2);
        % Correct all the times - remainder after division by 2 seconds
            spiketimes2=rem(spiketimes2,2);
        % Make PSTH
            starting=0;
            stopping=1.9999;
            binsize=0.010;
            % for each bin
            for i=1:200
                binstart=binsize*(i-1);
                binend=binsize*i;
                PSTH2(i)=sum(spiketimes2>binstart & spiketimes2<binend);
            end
            figure;plot([0:0.01:1.999],PSTH2)
            hold on;plot(time,rate,'r') 
% Compute z-vals (as above)
        newtrials2=[1 find(diff(spiketimes2)<0)+1 length(spiketimes2)];
        % for each trial
        z3=[];
        binwidth=0.01;
        for i=1:100
            spikes2=spiketimes2(newtrials2(i):newtrials2(i+1)-1);
            z3=[z3 computeksstat(time2,spikes2,PSTH2,binwidth)];
        end
        z3new=1-exp(-z3);
        mean(z3)=0.7374; 

        figure;hist(z3new) 
        figure;getKSPlot(z3new)

% Problem 5 - write code for a model that includes a bursting factor
                   % What do these terms look like?
                    t=[0:0.0001:0.1];
                    burstingfactor=(1+exp(2.5-((0.005-(t))/0.002).^2));
                    refractoryperiodfactor=exp(-0.002./(t));
        %%% only thing that changes in the code is this:
        refractoryperiodfactor(index:index+j)=exp(-0.002./(times(index:index+j)-times(index)));
        burstingfactor(index:index+j)=(1+exp(2.5-((0.005-(times(index:index+j)-times(index)))/0.002).^2));
        areaundercurve=trapz(times(index:index+j),rates(index:index+j).*refractoryperiodfactor(index:index+j).*burstingfactor(index:index+j));
% Simulate (as above)
        spiketimes4=simulatespiketrainBURST(time2,rate2);
        % Correct all the times - remainder after division by 2 seconds
        spiketimes4=rem(spiketimes4,2);
        % % Make PSTH
        binsize=0.010;
        % for each bin
        for i=1:200
            binstart=binsize*(i-1);
            binend=binsize*i;
            PSTH4(i)=sum(spiketimes4>binstart & spiketimes4<binend);
        end
        figure;plot([0:0.01:1.999],PSTH4)
        hold on;plot(time,rate,'r') 
% Test (as above)
        newtrials4=[1 find(diff(spiketimes4)<0)+1 length(spiketimes4)];
        % for each trial
        z4=[];
        binwidth=0.01;
        for i=1:100
            spikes4=spiketimes4(newtrials4(i):newtrials4(i+1)-1);
            z4=[z4 computeksstat(time2,spikes4,PSTH4,binwidth)];
        end
        z4new=1-exp(-z4);
        mean(z4)=0.7319; 

        % A
        figure;hist(z4new) 
        % B
        figure;getKSPlot(z4new)
