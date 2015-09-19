function WriteStringToLVBinFile(fid,Str);
%WriteStringToLVBinFile(fid,Str);

%first write the string length
fwrite(fid,length(Str),'uint32');
for ii=1:length(Str)
    fwrite(fid,Str(ii),'char');
end
return