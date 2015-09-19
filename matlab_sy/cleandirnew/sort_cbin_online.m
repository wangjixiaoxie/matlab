function foo
run_every_t=10;

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

!dir /B *.cbin > batchfoo
!dir /B *.rec > batchfoo_rec
pause(5);  % to let cbin files close
%cleandir_spect_SMART('batchfoo',1000,500,4,1) % right-hand behav rig
cleandir_spect_SMART('batchfoo',1000,500,4,4) % new (phys) rig - trying for R behav rig too
delbat('batchfoo.dcrd')
%movebat('batchfoo.dcrd','delete_files')
movebat('batchfoo.keep','keep_files')
