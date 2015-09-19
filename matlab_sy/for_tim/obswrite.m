function makesongfile(infile,outfile);

%makesongfile(infile,outfile);
%
%	MAKESONGFILE dumps vector INFILE to a observer-
%	playable file OUTFILE
%  (OUTFILE should be name in single quotes)ma
%   (Sample rate is not specified in file: observer default is normally 32000Hz).

%infile=(infile/max(infile))*100;  %may not be necc!

filename=[outfile];
[fp,m1]=fopen(filename,'w','b');
fwrite(fp,infile,'int16');
fclose(fp);
