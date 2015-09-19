function [] = db_batch_convert( batchfile, new_batchfile )
%db_batch_convert Takes a cbin.not.mat batch file and converts it to a
%.cbin batch file
%   Reads a cbin.not.mat file line by line and takes out the '.not.mat'
%   part. To be used for calculating entropy with db_contour

fprintf('\nConverting batch file...')

%opens original batch file to read
fid_old = fopen(batchfile, 'r');
%opens a new batch file to write
fid_new = fopen(new_batchfile, 'w');

%gets the first line of the old batch file
tline = fgetl(fid_old);

while ischar(tline)
    %creates a new line that removes the '.not.mat'
    new_tline = tline(1:strfind(tline, '.not.mat')-1);
    %prints the new line to the new batch file
    fprintf(fid_new, '%s\n', new_tline);
    %goes to the next line in the old batch file
    tline = fgetl(fid_old);
end

fclose(fid_old);
fclose(fid_new);

fprintf('done!\n') 

end

