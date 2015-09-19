function WritePathToLVBinFile(fid,pathname);

% split up the string
pp=findstr(pathname,':\');
if (~isempty(pp))
    maindir=pathname(1:pp(1)-1);
    pathname=pathname(pp(1)+2:end);
else
    maindir='C'; % default to avoid an error';
end


pathparts={maindir};
pp=findstr(pathname,'\');
while (1)
    pp=findstr(pathname,'\');
    if (~isempty(pp))
        pathparts{length(pathparts)+1}=pathname(1:pp(1)-1);
        pathname=pathname(pp(1)+1:end);
    else
        if (isempty(pathname))
            break;
        else
            pathparts{length(pathparts)+1}=pathname;
            break;
        end
    end
end

%labview first write PTH0
fwrite(fid,'PTH0','char');

StrLen=0;
for ii=1:length(pathparts)
    StrLen=StrLen + length(pathparts{ii}) + 1;
end
StrLen = StrLen + 4;

%first read the string length
fwrite(fid,StrLen,'int32');

%path type?
fwrite(fid,0,'char');fwrite(fid,0,'char');

%number of components?
fwrite(fid,length(pathparts),'int16');

for ii=1:length(pathparts)
    fwrite(fid,length(pathparts{ii}),'char');
    fwrite(fid,pathparts{ii},'char');
end
return
