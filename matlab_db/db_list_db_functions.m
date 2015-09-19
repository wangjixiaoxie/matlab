function [ listed_functions ] = db_list_db_functions( function_name )
%db_list_db_functions Lists the functions by DB (db_function_name) in a
%program.
%   Detailed explanation goes here

load_name = [function_name '.m'];

fid = fopen(load_name,'r');
current_line = fgetl(fid);
i=1;
j = 1;

while ischar(current_line)
    if ~isempty(strfind(current_line,'db_')) == 1 && strcmpi(current_line(1),'%') == 0 ...
        && isempty(strfind(current_line,function_name)) == 1

        spaces = strfind(current_line,' ');
        parens = strfind(current_line,'(');
        
        spaces = min(spaces(spaces > strfind(current_line,'db_')));
        parens = min(parens(parens > strfind(current_line,'db_')));
        
        if isempty(spaces) == 1
            function_end = parens;
        elseif isempty(parens) == 1
            function_end = spaces;
        elseif spaces > parens
            function_end = parens;
        else
            function_end = spaces;
        end
         
        listed_functions{j} = current_line(strfind(current_line,'db_'):function_end-1);
        
        j = j+1;
    end
    i=i+1;
    current_line = fgetl(fid);
end

fclose(fid);

if exist('listed_functions','var') == 0
    display('No db functions in this script')
    listed_functions = {};
    return
end

for i = 1:length(listed_functions)
    copies = find(strcmp(listed_functions{i},listed_functions) == 1);
    if copies == i
    else
        listed_functions{i} = '';
    end
end

listed_functions = listed_functions(~cellfun('isempty',listed_functions));


end

