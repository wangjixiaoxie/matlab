function [] = db_just_catch( yes_or_no )
%db_just_catch Similar to db_catch_and_clean but does not use cleandirAuto


% Find catch and get rid of noise files

batchname = 'batch';

%makes batch file
db_write_batch(batchname)

% batchname = input('What is the name of the batch file?  ', 's');

% only adds catch trials to the batch file
findcatch(batchname)
batchname = [batchname '.catch'];

% so you can examine how many songs are in the final batch file
edit(batchname)

%clears batchname
clear batchname

%run evsonganaly
if yes_or_no == 1
    evsonganaly
else
end

end

