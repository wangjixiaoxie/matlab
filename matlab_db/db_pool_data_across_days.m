function [ a ] = db_pool_data_across_days( start_date, end_date, batchfile, data_type, variable_name  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%If you entered a number for date_type (ex: start_date = 1,
%end_date = 10), it will load the first 10 entries in your batchfile (made
%using day_order function)
if isnumeric(start_date) == 1
    fid = fopen(batchfile,'r');
    current_line = fgetl(fid);
    line_number = 1;
    i = 1;
    
    while ischar(current_line)
        if line_number >= start_date && line_number <= end_date
            date_of_interest{i} = current_line;
            i = i+1;
        else
        end
        
        line_number = line_number+1;
        current_line = fgetl(fid);
    end
    
%If you entered a date for date_type (ex: start_date = 01Mar2013, end_date = 10Mar2013),
%it will load the entries with that date from the batchfile
elseif ischar(start_date) == 1
    fid = fopen(batchfile,'r');
    current_line = fgetl(fid);
    
    while ischar(current_line)
        if ~isempty(strfind(current_line,start_date)) == 1
            date_of_interest{1} = current_line;
            i = 2;
            current_line = fgetl(fid);
            while isnumeric(i)
                if isempty(strfind(current_line,end_date)) == 1
                    date_of_interest{i} = current_line;
                    i = i+1;
                    current_line = fgetl(fid);
                else
                    date_of_interest{i} = current_line;
                    i = 'end_of_date_range';
                end
            end
        end
        current_line = fgetl(fid);
    end
end

for
    

load(data)
a = who;

end

