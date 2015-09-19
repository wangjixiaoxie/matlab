%THIS SCRIPT IS ONLY FOR o57pu67 in 1/06

colordef black
figure
clear normfactor;
[infile,fs]=wavread('bu94bu15.20060807.0003');

outfile='bu94bu15.wav';
samp_fs=fs;
fs=44100;
rawsong=resample(infile,fs,samp_fs);



startsong=2.8;
endsong=5.8;
lowfiltval=75
highfiltval=20000

finlofiltval=150;
finhifiltval=12000;




datsong=rawsong(floor(fs*startsong):floor(fs*endsong));

%highpassfiltertoremovelow
[datsong]=bandpass(datsong,fs,lowfiltval,highfiltval,'butter');

normfactor(1)=max(datsong);
normfactor(2)=abs(min(datsong));
normfactor(1)=max(normfactor);

datsong=(1/normfactor(1))*datsong;

%writeaiff_file



%FINALLY, WRITE TO WAV FILES*/


    %corpshiftednormg{i}=[initbuff;corpshiftednormg{i};finalbuff];
    %resample(corpshiftednormg{i},44100,44150);
   wavwrite(datsong, 44100,16,outfile);
%for i=2:3
 %   %sngshift{i}=[initbuff;sngshift{i};finalbuff];
  %  wavwrite(sngshift{i}, fs,16,wavstring{i});
%end
