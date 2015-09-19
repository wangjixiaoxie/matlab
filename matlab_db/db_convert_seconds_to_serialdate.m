function [ output_vector ] = db_convert_seconds_to_serialdate( time_vector, input_sec_or_ser )
%db_convert_seconds_to_serialdate Will convert a vector in seconds to
%serial date time or vice versa. input_sec_or_ser specifies whether the
%input vector is in seconds 'sec' or serial date time 'ser'.


if strcmpi(input_sec_or_ser, 'sec')
    
    %Converts input vector in seconds to serial date time (86400 seconds in
    %one day)
    output_vector = time_vector .* (1/86400);

elseif strcmpi(input_sec_or_ser, 'ser')
    
    %Converts input vector in serial date time to seconds
    output_vector = time_vector .* 86400;
    
end


end

