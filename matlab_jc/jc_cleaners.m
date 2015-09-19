function array=jc_cleaners(input_array)
for i=1:10
    b=round(rand*(length(input_array)-1))+1;
    array(i)=input_array(b);
end