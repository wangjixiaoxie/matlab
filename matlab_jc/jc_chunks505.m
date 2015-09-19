function [pitcher]= jc_chunks505(input_array,on1,off1,on2,off2,on3,off3)
for i=1:length(input_array)
    left(i)=median(input_array(i).pitches(on1:off1));
    middle(i)=median(input_array(i).pitches(on2:off2));
    right(i)=median(input_array(i).pitches(on3:off3));
end
pitcher(1)=mean(left);
pitcher(2)=mean(middle);
pitcher(3)=mean(right);