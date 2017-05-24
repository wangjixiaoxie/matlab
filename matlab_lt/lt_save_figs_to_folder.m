function lt_save_figs_to_folder(foldername, closeall)
%% lt 5/4/17 - save all figs to folder, labeled with time of save and foldername

if ~exist('foldername', 'var') % name to tag folder
    foldername = '';
end

if ~exist('closeall', 'var') % close all figs once saved?
    closeall = 0;
end

%%
currdir = pwd;
savedir = '/bluejay5/lucas/FIGS/save_figs_to_folder';
cd(savedir);

tstamp = lt_get_timestamp(0);
newfoldername = [foldername '_' tstamp];

mkdir(newfoldername)
cd(newfoldername);

lt_save_all_figs;

cd(currdir)

if closeall==1
    close all;
end

