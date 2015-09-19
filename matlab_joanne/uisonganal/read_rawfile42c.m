function [outdata, ad_freq] = read_rawfile42c(file,datatype);

%  [outdata, ad_freq] = read_rawfile(file,datatype);
%
%  READ_RAWFILE42c is version of read_rawfile that runs with 4.2c
%  if there are problems, fiddle with 'int16' vs 'short' as datatype
%  
%
%  READ_RAWFILE can read old-style foosong files and songfilt-
%  produced stimulus files and also msong datafiles. FILE is 
%  the name of the datafile to be read. DATATYPE is a string
%  (e.g. 'int16' or 'int32' or 'double') describing the type of
%  binary data stored in FILE. Data read is returned as
%  OUTDATA, and AD_FREQ is the A/D sampling rate of OUTDATA.
%  





PROPORT_FREQ = 32000;  % sample rate for proport data.

[fid,errmsg] = fopen(file,'r','b');
if (fid == -1),
  error(errmsg);
end

%adstring = fgets(fid,8);
adstring = fread(fid,8,'uchar');
adstring=char(adstring);
adstring=adstring';
if strcmp(adstring,'AD_FREQ:'),   
  % it's an msong file or a songfilt stim file

  frewind(fid);  adstring = fgetl(fid);
  [ad_freq,n,errmsg]=sscanf(adstring, 'AD_FREQ: %d Hz');
  if n ~= 1,
    error(errmsg)
  end

else   % it's an old-style foosong file
  ad_freq=PROPORT_FREQ;
  frewind(fid);
end 

len = fread(fid,1,'int');
[outdata(:,1),count] = fread(fid,datatype);

if (count ~= len),
	fclose(fid);
	error('read_rawfile: # of data points read and # specified in file do not agree!');
end

fclose(fid);

