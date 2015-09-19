function [y,Fs,Format]=wavread16(wavefile)
%
%   [y]=WAVREAD(wavefile) loads a .WAV format file specified by "wavefile", 
%       returning the sampled data in variable "y". The .WAV extension 
%       in the filename is optional.
%
%   [y,Fs]=WAVREAD(wavefile) loads a .WAV format file specified by  
%       "wavefile", returning the sampled data in variable "y" and the 
%       sample rate in variable "Fs".
%   
%   [y,Fs,Format]=WAVREAD(wavefile) loads a .WAV format file specified by 
%       "wavefile",returning the sampled data in variable "y", the sample 
%       rate in variable "Fs", and the .WAV file format information in 
%       variable "Format". The format information is returned as a 6 element
%       vector with the following order:
%
%               Format(1)       Data format (always PCM) 
%               Format(2)       Number of channels
%               Format(3)       Sample Rate (Fs)
%               Format(4)       Average bytes per second (sampled)
%               Format(5)       Block alignment of data
%               Format(6)       Bits per sample


if nargin~=1
	error('WAVREAD takes one argument, which is the name of the .WAV file');
end

%if findstr(wavefile,'.')==[]
%	wavefile=[wavefile,'.wav'];
%end

[fid,msg]=fopen(wavefile,'r','l');
%fid=fopen(wavefile,'r');
if fid ~= -1 
	% read riff chunk
	header=fread(fid,4,'uchar');		% 'RIFF'
	header=fread(fid,1,'ulong');		% length of file
	header=fread(fid,4,'uchar');		% 'WAVE'
	
	% read format sub-chunk
	header=fread(fid,4,'uchar');		% 'fmt '
	header=fread(fid,1,'ulong');		% length of fmt chunk
	
	Format(1)=fread(fid,1,'ushort');	% 1 == PCM
	Format(2)=fread(fid,1,'ushort');	% 1 == mono ; 2 == stereo
	Fs=fread(fid,1,'ulong');		% samples per second
	Format(3)=Fs;
	Format(4)=fread(fid,1,'ulong');		% average bytes per second
	Format(5)=fread(fid,1,'ushort');	% bytes per sample
	Format(6)=fread(fid,1,'ushort');	% bits per sample
	
	
	% read data sub-chunck
	header=fread(fid,4,'uchar');		% 'data'
	nsamples=fread(fid,1,'ulong');		% length of data
	if Format(5) == 1
		y=fread(fid,nsamples,'uchar');	% 8-bits unsigned
	else
		y=fread(fid,nsamples,'short');	% 16-bits signed
	end
	fclose(fid);
end     

if fid == -1
  disp(msg);
  error(['Can''t open .WAV file for input: ' wavefile]);
end;
