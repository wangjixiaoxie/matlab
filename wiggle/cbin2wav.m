function cbin2wav(cbin,targetstr)
%
% writes the song waveform in a cbin file to a wav file specified by
% targetstr, or the cbin file's name + .wav if targetstr is empty

% cbin2wav('r4w53_060612_1828.234.cbin','r4w53afterWN')
% go to the folder that the file resides and do above 
% cbin = cbin file name with ''; targetstr is name of wav file with ''

if(nargin < 2)
    targetstr = cbin(1:(end-4));
    targetstr = [targetstr 'wav'];
end

[rawsong FS] = ReadCbinFile(cbin);
rawsong = rawsong(:,1);
rawsong = rawsong / max(rawsong);

wavwrite(rawsong,FS,16,targetstr);

return;