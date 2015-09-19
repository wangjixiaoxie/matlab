%this will make an nx2 matrix which can then be used for binning data into
%time bins.
function [matrixvals]=make_time_indices (start_day, end_day, start_times, end_times) 
    
    dayvalscomb=[];
    dayvals={};

    for ii=1:length(start_times)

        dayvals{ii}=datenum(start_day, 'yyyy-mm-dd'):1:datenum(end_day,'yyyy-mm-dd');
        dayvals2{ii}(:,1)=dayvals{ii}+start_times(ii)/24
        dayvals2{ii}(:,2)=dayvals{ii}+end_times(ii)/24
    end
    
    %now combine the dayvals columns and sort them into one.
    for ii=1:length(start_times)
        dayvalscomb=[dayvalscomb; dayvals2{ii}]
    
    end
    matrixvals=sort(dayvalscomb)
end