function StrOut=ReadStringFromLVBinFile(fid);

%first read the string length
StrLen=fread(fid,1,'int32');

%setup a string
StrOut=blanks(StrLen);

%read into an array of char
tmp=fread(fid,StrLen,'char');
for ii=1:StrLen
    StrOut(ii)=tmp(ii);
end
return