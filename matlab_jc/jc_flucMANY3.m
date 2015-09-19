function [ratio1,ratio2,ratio3]=jc_flucMANY3

    period=1:1:200; % 1 to 200ms period
    fs = 32000; % cbin sampling rate
    t = 0:1/fs:1; %32000pts
    
for i=1:50
    % Create artificial signal with the appropriate period using vco
        modu=cos((1000/period(i))*2*pi*t);
        est=modu*100+3000;
        x = vco(modu,[2900 3100],fs); % 2400 and 3600 are the minimum and maximum frequency values that the pitch oscillates between
    % Zero pad to allow for sonogram window
        xfinal=[zeros(1,512) x zeros(1,512)];
    % Calculate pitch with gaussian windows of various widths
%         pitch1=jc_pitchmat1024(xfinal,1024,1020,0.5,2600,3400,1,'obs0',1);
        pitch2=jc_pitchmat1024(xfinal,1024,1020,1,2600,3400,1,'obs0',1);
%         pitch3=jc_pitchmat1024(xfinal,1024,1020,2,2600,3400,1,'obs0',1);
    % "resample" estimate to match overlap of gaussian windows
        %estfinal=est([1:4:length(est)]);
    % Calculate psd with matlab periodogram - fs/4 because of overlap of
    % gaussian windows
%         [x,psdPRED1]=jcpsd(pitch1,fs/4);
        %[x,psdPRED2]=jcpsd(pitch2,fs/4);
%         [x,psdPRED3]=jcpsd(pitch3,fs/4);
        %[x,psdACT]=jcpsd(estfinal,fs/4);
    % Determine the point at which the gaussian windows smooth too much
%         ratio1(i)=max(psdPRED1)/max(psdACT);
        %ratio2(i)=max(psdPRED2)/max(psdACT);
        ratio2(i)=(max(pitch2)-min(pitch2))./(max(est)-min(est));
%         ratio3(i)=max(psdPRED3)/max(psdACT);
end
g=7;
