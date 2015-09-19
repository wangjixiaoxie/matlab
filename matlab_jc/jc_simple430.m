function [new_HistoMatrix]=jc_simple430(input_matrix,lower_limit,higher_limit)
new_HistoMatrix=0;
for i=1:length(input_matrix) 
    if (input_matrix(i)>lower_limit & input_matrix(i)<higher_limit) 
        new_HistoMatrix=[new_HistoMatrix input_matrix(i)];
    end
end
new_HistoMatrix=new_HistoMatrix(1:length(new_HistoMatrix));