function [] = get_rec(fileroot, indmin, indmax, outprefix);
%
% function [] = get_rec(fileroot, indmin, indmax, outprefix);
%
% GET_REC  Read *.rec files with the rootname FILEROOT and
%    with trial indeces between and including INDMIN and
%    INDMAX.  Creates batch files in the form:
%
%       outprefix.stimulus_type
%       
%

TRUE = 1;
FALSE = 0;
fbname = '';
fb_ind = 1;

% Loop through *.rec files and scan for the feedback names and 
%  start times.
for i=indmin:indmax,
   if (i < 10),
      insert = '00';
   elseif i < 100,
      insert = '0';
   else
      insert = '';
  end

   if (fileroot(length(fileroot)) == '.'),
      cur_base = sprintf('%s%s%d.', fileroot, insert, i);
   else
      cur_base = sprintf('%s.%s%d.', fileroot, insert, i);
   end
   recfile = [cur_base, 'rec'];
   bbinfile = [cur_base, 'bbin'];
   
	[fid,errmsg] = fopen(recfile, 'r');
   if (fid == -1),
      %disp(['Unable to open file: ' recfile]);
   else
      scan_it = FALSE;
      
		% Get stimulus type
      while 1
        line = fgetl(fid);
        if ~isstr(line), break, end
        if (findstr(line,'Stimulus:')),                 % find appropriate line
           spaces = findstr(line,' ');                  % find spaces
           rec_type = line(spaces(1)+1 : spaces(2)-1);  % pull out stimulus type
       
           outfile = [outprefix,'_',rec_type];               
           batchid=fopen(outfile,'a');                  % open appropriate batch
           fprintf(batchid,'%s\n',bbinfile);             % write in current .rec file
           fclose(batchid);                             % close it
        end        
          
      end
		  fclose(fid);
   end
end
