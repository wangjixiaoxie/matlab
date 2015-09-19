function [ratio1]=jc_flucMANYifd(index)

    period=1:1:200; % 1 to 200ms period
    fs = 32000; % cbin sampling rate
    t = 0:1/fs:1; %32000pts
    
for i=1:length(index) %1:50

    % Create artificial signal with the appropriate period using vco
        modu=cos(round((1000/index(i))*2*pi)*t);
        est=modu*100+3000;
        x = vco(modu,[2900 3100],fs); % 2400 and 3600 are the minimum and maximum frequency values that the pitch oscillates between
    % Zero pad to allow for sonogram window
        xfinal=[zeros(1,512) x zeros(1,512)];
    % Calculate pitch with gaussian windows of various widths
        pitch1=jc_PitchData208(xfinal,1024,1020,2,2600,3400,1,'obs0');
    % "resample" estimate to match overlap of gaussian windows
        estfinal=est([1:4:length(est)]);
    % Calculate psd with matlab periodogram - fs/4 because of overlap of
    % gaussian windows
        [x,psdPRED1]=jcpsd4(pitch1(1000:7000),fs/4);
        [x,psdACT]=jcpsd4(estfinal(1000:8000),fs/4);
    % Determine the point at which the gaussian windows smooth too much
        ratio1(i)=max(psdPRED1)/max(psdACT);
end
