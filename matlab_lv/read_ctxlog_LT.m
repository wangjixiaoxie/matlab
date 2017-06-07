function [alltransitions, allnotegroup] = read_ctxlog()
%% 6/6/17 - lt modified to give datenums for alltransitions (insted of datetime structure), for older version ofm atlab

%%
alltransitions = [];
allnotegroup = [];

ctxlog_files =  dir('*.ctxlog');

if length(ctxlog_files)>1
    for i=1:length(ctxlog_files)
        startdate = ctxlog_files(i).name(16:21);
        startdates_dn(i) = datenum(startdate,'MMDDYY');
    end
    [~,ix] = sort(startdates_dn);
else
    ix = 1;
end
%working for one ctxlog file for now
for i=1:length(ctxlog_files)
    filename=ctxlog_files(ix(i)).name;
    %sprtintf('processing %s',filename)

    fid1=fopen(filename);
    
    data=textscan(fid1,'%d %d %d %s %d %d %d %s %s %s %d');
    try
    transitions = datetime(data{3},data{1},data{2},data{5},data{6},data{7});
    catch err
        % < 2014b, doesn't have this function
        datearray = [data{3},data{1},data{2},data{5},data{6},data{7}];
        datearray = double(datearray);
        transitions = datenum(datearray);
    end
    notegroup = data{11};
    
    alltransitions = [alltransitions;transitions];
    allnotegroup = [allnotegroup; notegroup];
end

fclose(fid1);