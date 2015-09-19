function StrOut=ReadPathFromLVBinFile(fid);

%labview first write PTH0
fread(fid,4,'char');

%first read the string length
StrLen=fread(fid,1,'int32');

%path type?
PType=fread(fid,2,'char');

%number of components?
NComp=fix(fread(fid,1,'int16'));

%setup a string
StrOut='';

for ii=1:NComp
    StrLen=fix(fread(fid,1,'char'));
    tmp=fread(fid,StrLen,'char');
    for jj=1:StrLen
        StrOut=[StrOut,tmp(jj)];
    end
    if (ii==1)
        StrOut=[StrOut,':\'];
    elseif (ii<NComp)
        StrOut=[StrOut,'\'];
    end
end
return
