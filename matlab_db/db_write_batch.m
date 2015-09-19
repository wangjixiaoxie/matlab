function [ ] = db_write_batch( batch )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fid = fopen(batch,'w');

filenames = dir('*.cbin');

for i = 1:length(filenames)
	fprintf(fid,'%s\n',filenames(i).name);
end

fclose(fid);

end

