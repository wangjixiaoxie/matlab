function lt_goto_bluejay(argument)

if argument==0;
    cd(['/bluejay/lucas/birds'])
else
    try
        cd(['/bluejay' num2str(argument) '/lucas/birds']);
    catch err
        disp('error: that bluejay number does not exist as directory')
    end
end
end
