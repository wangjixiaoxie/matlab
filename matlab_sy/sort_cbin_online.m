function foo
run_every_t=300;

tic
first=1;
while 1
    if first | toc>run_every_t
run_func
        pause(.1)
        tic
        first=0;
    end
end

function run_func

!dir /B **.cbin > batchfoo
%cleandir_spect('batchfoo',10^4,500,4,1)
%cleandir_spect('batchfoo',10^6,500,4,1)
%cleandir_spect('batchfoo',10^4,500,4,1)

cleandir_spect_SMART('batchfoo',10^5,500,4,3,)
%cleandir_spect('batchfoo',10^5,500,4,3)

delbat('batchfoo.dcrd')
%movebat('batchfoo.dcrd','not_song')
movebat('batchfoo.keep','keep_files')
