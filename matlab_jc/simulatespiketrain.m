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
    spiketimes(i)=times(index);
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
