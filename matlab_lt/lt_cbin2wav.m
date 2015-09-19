function cbin2wav(cbin,targetstr)
%
% writes the song waveform in a cbin file to a wav file specified by
% targetstr, or the cbin file's name + .wav if targetstr is empty
%
%

if(nargin < 2)
    targetstr = cbin(1:(end-4));
    targetstr = [targetstr 'wav'];
end

[rawsong FS] = ReadCbinFile(cbin);
rawsong = rawsong(:,1);
rawsong = rawsong / max(rawsong);

wavwrite(rawsong,FS,16,targetstr);

return;