function [dtnum datestring]=lt_neural_fn2datenum(fname)
% fname = 'bk7_1860_161005_125912.rhd';
% dtnum = 736608.541111111
% datestring = '161005_125912'

uscores=findstr(fname, '_');

assert(strcmp(fname(uscores(end)+7:end), '.rhd'), 'MISTAKE?');

datestring=fname(uscores(end-1)+1:end-4);
% disp(['date_time = ' datestring]);

dtnum=datenum(datestring, 'yymmdd_HHMMSS');

