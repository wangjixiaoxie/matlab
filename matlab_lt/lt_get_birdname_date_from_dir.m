%% LT 3/15/14 
% start from within a day folder (e.g.) /bluejay3/lucas/birds/pu11wh87/030114_CPseq_durWN_day1
% and get birdname and date.  paths have to be structured as above.

function [varargout]=lt_get_birdname_date_from_dir(input)

% input =0 in the bird folder
% input =1 Day folder
% input =0 have [birdname bluejaynum] as output; 
% input =1, have [birdname bluejaynum date(cell array) phrase] as outputs

dirpath=pwd;
slashes=findstr(dirpath,'/');
bluejay=findstr(dirpath,'bluejay');

if input==0;
    bird_name=dirpath(slashes(end)+1:end);
    
    bluejay_num=dirpath(bluejay+7);
    
    varargout{1}=bird_name;
    varargout{2}=bluejay_num;
end

if input==1;
    bird_name=dirpath(slashes(end-1)+1:slashes(end)-1);
    
    date_str_mmddyy=dirpath(slashes(end)+1:slashes(end)+6);
    date_num=datenum(date_str_mmddyy,'mmddyy');
    date_str_ddmmmyyyy=datestr(date_num,'ddmmmyyyy');
    date{1}=date_num;
    date{2}=date_str_ddmmmyyyy;
        
    bluejay_num=dirpath(bluejay+7);

    underscore_after_slashes=findstr(dirpath(slashes(end)+1:end),'_');
    if length(underscore_after_slashes)==1;
    phrase=dirpath(slashes(end)+8:end);
    elseif length(underscore_after_slashes)>1;
        phrase = dirpath(slashes(end)+8:slashes(end)+underscore_after_slashes(2)-1);
    end
    
    varargout{1}=bird_name;
    varargout{3}=date;    
    varargout{2}=bluejay_num;
    varargout{4}=phrase;

end



end
