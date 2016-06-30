%% LT 1/6/14 - changed to save batch as batch.CATCH.keep, so it does not overlap with batch.catch.keep, which is made by db_transfer_calls
% use 1

function [] = db_catch_and_clean( yes_or_no )
%db_catch_and_clean Makes a batch of the cbin files, finds the catch
%trials, uses cleardirAuto to get rid of noise files, and loads evsonganaly
%if you enter 1
%   Detailed explanation goes here

% Find catch and get rid of noise files

batchname = 'batch';
% batchname = input('What is the name of the batch file?  ', 's');

%makes a batch file
db_write_batch(batchname)

% only adds catch trials to the batch file
findcatch(batchname)
batchnameOLD = [batchname '.catch'];
batchname = [batchname '.CATCH'];

eval(['!mv ' batchnameOLD ' ' batchname]);  % renames batch.catch to batch.CATCH
   

% cleans batch file to contain only suspected songs
lt_cleandirAuto(batchname,1000,4,4)
batchname = [batchname '.keep'];

% so you can examine how many songs are in the final batch file
edit(batchname)

%clears variables
clear batchname

%runs evsonganaly
if yes_or_no == 1
    evsonganaly
else
end

end

