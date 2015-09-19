%% LT 6/1/14 - open any 1_NOTE_TO_SELF.txt for any bird, at any location
% e.g. lt_open_any_bird_notefile('gr66gr43')
function lt_open_any_bird_notefile(bname)


% RUN thru 4 hard disks to try to open note file
try
    open(['/bluejay1/lucas/birds/' bname '/1_NOTE_TO_SELF.txt'])
catch err
    try
        open(['/bluejay2/lucas/birds/' bname '/1_NOTE_TO_SELF.txt'])
    catch err
        try
            open(['/bluejay3/lucas/birds/' bname '/1_NOTE_TO_SELF.txt'])
        catch err
            try
                open(['/bluejay4/lucas/birds/' bname '/1_NOTE_TO_SELF.txt'])
            catch err
                disp('note file/bird/something does not exist')
            end
        end
    end
end