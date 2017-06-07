function plot_seq_context_LT(exptname, birdname, transitions)
%% LT - modified 6/6/17 to run not just for lena specific bird
% REQUIREMENTS:
% songs go in birdname/exptname/day folder (e.g. /bluejay5/lucas/birds/wh8pk40/contextseq/170605
% transitions go in birdname/birdname folder: e.g. /bluejay5/lucas/birds/wh8pk40/wh8pk40/all_days_transition_matrix_contextseq
% ctxlog go in any and multiple day folders (will compile all and take unique)

%%
% exptname = 'contextseq';
% dayrange = 05:06;
% cdwhere = '170606'; %name of last day folder (i.e. folder that holds ctxlog)
% birdname = 'wh8pk40';

trn1 = transitions{1};
trn2 = transitions{2};

%%
%plot context changes and sequence transition probabilies


if (0)
load lenacols2
lenacols = cols;
end
cols = [1 .7 .2; .8 .5 .1; .5 .8 .3; .1 .6 .3];
colos = [1 .9 .4; 1 .9 .3; .5 1 .5; .5 1 .5];

%% get all_trans
% 
cd(birdname)
cd(['all_days_transition_matrix_' exptname '/']);

%% COMPILE ALL DAYS INTO ONE STRUCT

% --- 1) get first day to initiate struct
TargStr = [birdname '_*'];
Dates_first = '';
Dates_last = '';
FileType = 'mat';
make_batch = 0;

[~, cellofnames, Dates_first, Dates_last, ~, arrayofdatenums]=lt_write_all_folder_contents_to_batch_v2(...
    TargStr,Dates_first,Dates_last, FileType, make_batch)

% -- load first day
load(cellofnames{1});

% day = num2str(dayrange(1));
%     if length(day)==1
%         day = ['0' day];
%     end
%     
% 
%     
%     load(['rd82wh13_' day 'Sep2016_test.mat'])

aat = all_trans;
aat.songday = ones(size(aat.syl_times));

% ----- 2) put syllables for all days  and all songs together\
for i = 2:length(cellofnames)
load(cellofnames{i});
aat.syl_order_compiled_across_songs = [aat.syl_order_compiled_across_songs all_trans.syl_order_compiled_across_songs];
aat.syl_times = [aat.syl_times all_trans.syl_times];
aat.syl_order = [aat.syl_order all_trans.syl_order];
aat.syl_order_WithDashes = [aat.syl_order_WithDashes all_trans.syl_order_WithDashes];
aat.songday = [aat.songday ones(size(all_trans.syl_times)).*i];
end


% for i = 2:length(dayrange)
%     day = num2str(dayrange(i));
%     if length(day)==1
%         day = ['0' day];
%     end
% load(['rd82wh13_' day 'Sep2016_test.mat']);
% aat.syl_order_compiled_across_songs = [aat.syl_order_compiled_across_songs all_trans.syl_order_compiled_across_songs];
% aat.syl_times = [aat.syl_times all_trans.syl_times];
% aat.syl_order = [aat.syl_order all_trans.syl_order];
% aat.syl_order_WithDashes = [aat.syl_order_WithDashes all_trans.syl_order_WithDashes];
% aat.songday = ones(size(all_trans.syl_times)).*i; % NOTE: mistake here (LT)
% end
% 
all_trans = aat;


%% get transitions (lt modified to get ctxlog from all days
cd ..
cd ..
cd(exptname);
% go to last day to get ctxlog
transtime = [];
switchtonotegroup = [];
for i=1:length(arrayofdatenums)
   cd(datestr(arrayofdatenums(i), 'yymmdd'));
   [tr, sw] = read_ctxlog();
   
   transtime = [transtime; tr];
   switchtonotegroup = [switchtonotegroup; sw];

   cd ..
end

% pull out unique times (in case same ctxlog in multiple folders)
[~, inds] = unique(transtime);
transtime = transtime(inds);
switchtonotegroup = switchtonotegroup(inds);

% cd(cdwhere)
% [transtime, switchtonotegroup] = read_ctxlog();

%%
firstday_dn = datenum(aat.parameters.date_of_exp);

%convert Lucas time to real time, have to guess first day?
% % % first_day_datenum=datenum(startdate,'ddmmmyyyy');
% % % songtime_dn = alldata1(:,2)+first_day_datenum-1;
%songtime:
firstsong_dn = all_trans.syl_times{1}(1);
lastsong_dn = all_trans.syl_times{end}(end);
transtime_dn = datenum(transtime);

%get transitions for these songs
startix = find(transtime_dn<firstsong_dn,1,'last'); % first trn
endix = find(transtime_dn>lastsong_dn,1,'first'); % last trn
if isempty(startix)
    startix = 1;
end
if isempty(endix)
    endix = length(transtime_dn);
end


transtime_dn = transtime_dn(startix:endix);
switchtonotegroup = switchtonotegroup(startix:endix);
%convert to work with time


%find actual switches (leftover from many 1-1 switches before had
%NGspecific timing
realswitchix = [1; find(diff(switchtonotegroup))+1];
% transtime = transtime(realswitchix);
transtime_dn = transtime_dn(realswitchix);
switchtonotegroup = switchtonotegroup(realswitchix);

transtime_t0 = (transtime_dn-firstday_dn)*24;

%% assign notegroups to syllable times
%for each transtime_dn, find all following syllables that are smaller than the next, take that
%switchtonotegroup
%syltimes; notegroups; 
syltimes = [all_trans.syl_times{:}];

notegroups = nan(size(syltimes));
for k = 1:length(transtime_dn)-1
notegroups(syltimes>transtime_dn(k) & syltimes<transtime_dn(k+1)) = switchtonotegroup(k);
end

%% assign notegroups to song times
songtimes = cellfun(@mean, all_trans.syl_times);

notegroups_songs = nan(size(songtimes));
for k = 1:length(transtime_dn)-1
    notegroups_songs(songtimes>transtime_dn(k) & songtimes<transtime_dn(k+1)) = switchtonotegroup(k);
end
%% copied this from plot_ctxtchanges
%% get by song info from allalldata_bysong
firstday_dn = datenum(aat.parameters.date_of_exp);
datenum_song = songtimes - firstday_dn +1; %should be same as alldata1(:,2) now
dv_song = datevec(datenum_song);
day_song = dv_song(:,3)+(dv_song(:,2)-1)*31;
time_song = (datenum_song-1)*24*60*60*1000; %time of song in ms


%% SMOOTHED PROB OVER TIME
% get transition probability of BJ
% times; transngs; transprob
sylorder = all_trans.syl_order_compiled_across_songs;
% bjix = find(sylorder(1:end-1)=='b' & sylorder(2:end)=='j');
bjix = find(sylorder(1:end-1)==trn1(1) & sylorder(2:end)==trn1(2));
bhix = find(sylorder(1:end-1)==trn2(1) & sylorder(2:end)==trn2(2));

% -- fill out each syl time with 0 or 1 or nan depending on if correspond
% to trn of interest.
transprob = nan(1,length(syltimes));
transprob(bjix) = 1;
transprob(bhix) = 0; 
transprob(isnan(transprob))=[];
xix = sort([bjix bhix]);
times = syltimes(xix);
transngs = notegroups(xix);

% --- plot smoothed trn prob
figure
plot(times, lt_smooth(transprob, 10),'k-')
% plot(times,smooth_st(transprob,30,'gauss'),'k-')

hold on
yl = get(gca,'ylim');
for i = 1:length(transtime_dn)-1
    line([transtime_dn(i) transtime_dn(i)],yl,'color','k')
    fill([transtime_dn(i) transtime_dn(i+1) transtime_dn(i+1) transtime_dn(i)],[yl(1) yl(1) yl(2) yl(2)],colos(switchtonotegroup(i)+1,:),'edgecolor','none','facealpha',.5)
    text(transtime_dn(i),yl(1)+0.5*(yl(2)-yl(1)),num2str(switchtonotegroup(i)))
end

xl = get(gca,'xlim');
    fill([transtime_dn(i) xl(2) xl(2) transtime_dn(i)],[yl(1) yl(1) yl(2) yl(2)],colos(switchtonotegroup(i)+1,:),'edgecolor','none','facealpha',.5)
    ylabel('p(B-J) / p(B-H|B-J)')
xlabel('time')
% plot(times,lt_smooth(transprob, 10),'k-')
set(gca,'xtick',transtime_dn(1:10:end),'xticklabel',datestr(transtime_dn(1:10:end)));
box off
axis tight



%% block stats (COLLECT OVER ALL BLOCKS)
[allmean allsem] = grpstats(transprob,transngs,{'mean' 'sem'})
figure
subplot(1,2,1)
title('mean/sem')
hold on
errorbar(1:length(allmean),allmean,allsem,'linestyle','none','linewidth',2)
tmp = unique(transngs);
tmp = mat2cell(tmp(~isnan(tmp)), 1, ones(1, length(tmp(~isnan(tmp)))));
set(gca,'xtick',1:length(allmean),'xticklabel',tmp)
% ylabel('p(BJ/ BJ+BH)')
ylabel([trn1 '/(' trn1 '+' trn2 ')'])


%vectors of all transitions for each probe
p11 = transprob(transngs==0);
p21 = transprob(transngs==2);
%blockstats figure
%assign block numbers
ctxt_starts = [1 find(diff(transngs))+1];
ctxt_dur = [diff(ctxt_starts) length(transngs)-ctxt_starts(end)];
blocknumber = [];
for i = 2:length(ctxt_starts)
        blocknumber = [blocknumber; ones(ctxt_starts(i)-ctxt_starts(i-1),1)*(i-1)];
end
    blocknumber = [blocknumber; ones(length(transprob)-ctxt_starts(i)+1,1)*i];

meantimes = grpstats(times,blocknumber);
[fmean, fnumel, fstd, fsem] = grpstats(transprob,blocknumber,{'mean','numel','std','sem'});
ngblock = grpstats(transngs,blocknumber);

%trans probs for each block
minblocksize = 3;
p1 = fmean(ngblock==0);
p2 = fmean(ngblock==2);
p1(fnumel(ngblock==0)<minblocksize)=[];
p2(fnumel(ngblock==2)<minblocksize)=[];
plot(ones(size(p1)),p1,'o','color',cols(1,:),'linewidth',2)
plot(ones(size(p2))+2,p2,'o','color',cols(3,:),'linewidth',2)
text(3,0.7,num2str(ranksum(p1,p2)))

subplot(1,2,2)
hold on
for thisng = 0:3
    ix = ngblock==thisng;
    errorbar(meantimes(ix),fmean(ix),fsem(ix),'linestyle','none','color',cols(thisng+1,:),'linewidth',4-2*rem(thisng,2))
end
axis tight

%% proportion test
up = mean(p11)-mean(p21);
phat = (sum(p11)+sum(p21))/(length(p11)+length(p21));
down = sqrt(phat*(1-phat)*(1/length(p11)+1/length(p21)));
up/down


%% song by song
clear thissongtime propbj
for i = 1:length(all_trans.syl_order)
    thissong = all_trans.syl_order{i};
    thissongtime(i) = mean(all_trans.syl_times{i});
    
    bjn = sum(thissong(1:end-1)==trn1(1) & thissong(2:end)==trn1(2));
    bhn = sum(thissong(1:end-1)==trn2(1) & thissong(2:end)==trn2(2));
    propbj(i) = bjn/(bjn+bhn);
    
end

figure
plot(thissongtime,smooth_lv(propbj,10),'k')
hold on
yl = get(gca,'ylim');
for i = 1:length(transtime_dn)-1
    line([transtime_dn(i) transtime_dn(i)],yl,'color','k')
    fill([transtime_dn(i) transtime_dn(i+1) transtime_dn(i+1) transtime_dn(i)],[yl(1) yl(1) yl(2) yl(2)],colos(switchtonotegroup(i)+1,:),'edgecolor','none','facealpha',.5)
    text(transtime_dn(i),yl(1)+0.5*(yl(2)-yl(1)),num2str(switchtonotegroup(i)))

end
xl = get(gca,'xlim');
fill([transtime_dn(i) xl(2) xl(2) transtime_dn(i)],[yl(1) yl(1) yl(2) yl(2)],colos(switchtonotegroup(i)+1,:),'edgecolor','none','facealpha',.5)
ylabel('p(B-J) / p(B-H|B-J)')
xlabel('time')
plot(thissongtime,smooth_lv(propbj,10),'kx')
set(gca,'xtick',transtime_dn(1:10:end),'xticklabel',datestr(transtime_dn(1:10:end)));
box off
axis tight

disp('no dashes')
grpstats(propbj,notegroups_songs)


%% song by song with dashes
clear thissongtime propbj_dash
for i = 1:length(all_trans.syl_order_WithDashes)
    thissong = all_trans.syl_order_WithDashes{i};
    thissongtime(i) = mean(all_trans.syl_times{i});
    
    bjn = sum(thissong(1:end-1)== trn1(1) & thissong(2:end)==trn1(2));
    bhn = sum(thissong(1:end-1)==trn2(1) & thissong(2:end)==trn2(2));
    propbj_dash(i) = bjn/(bjn+bhn);
    
    
end

figure
plot(thissongtime,smooth_lv(propbj_dash,10),'k')
hold on
yl = get(gca,'ylim');
for i = 1:length(transtime_dn)-1
    line([transtime_dn(i) transtime_dn(i)],yl,'color','k')
    fill([transtime_dn(i) transtime_dn(i+1) transtime_dn(i+1) transtime_dn(i)],[yl(1) yl(1) yl(2) yl(2)],colos(switchtonotegroup(i)+1,:),'edgecolor','none','facealpha',.5)
    text(transtime_dn(i),yl(1)+0.5*(yl(2)-yl(1)),num2str(switchtonotegroup(i)))

end
xl = get(gca,'xlim');
fill([transtime_dn(i) xl(2) xl(2) transtime_dn(i)],[yl(1) yl(1) yl(2) yl(2)],colos(switchtonotegroup(i)+1,:),'edgecolor','none','facealpha',.5)
ylabel('p(B-J) / p(B-H|B-J)')
xlabel('time')
plot(thissongtime,smooth_lv(propbj_dash,10),'kx')
set(gca,'xtick',transtime_dn(1:10:end),'xticklabel',datestr(transtime_dn(1:10:end)));
box off
axis tight
title('with dash')

disp('dashes')
grpstats(propbj_dash,notegroups_songs)


%% song by song with dashes proper probs
% NOTE !! - NEED TO MODIFY TO DO CONVERGENT AS WELL

for i = 1:length(all_trans.syl_order_WithDashes)
    thissong = all_trans.syl_order_WithDashes{i};
    thissongtime(i) = mean(all_trans.syl_times{i});
    
    bpreidx = thissong(1:end-1)=='b';
    postsyls = unique(thissong(bpreidx+1));
%     disp(postsyls)
    
    bjn = sum(bpreidx & thissong(2:end)=='j');
    bhn = sum(bpreidx & thissong(2:end)=='h');
    bxn = sum(bpreidx & thissong(2:end)=='-');
    alln = bjn+bhn+bxn;
    allprop_bj(i) = bjn/alln;
    allprop_bh(i) = bhn/alln;
    allprop_bx(i) = bxn/alln;
end

figure
hold on

plot(thissongtime,smooth_lv(allprop_bj,10),'k-x')
plot(thissongtime,smooth_lv(allprop_bh,10),'r-x')
plot(thissongtime,smooth_lv(allprop_bx,10),'g-x')
legend({'BJ' 'BH' 'BX'})

yl = get(gca,'ylim');
for i = 1:length(transtime_dn)-1
    line([transtime_dn(i) transtime_dn(i)],yl,'color','k')
    fill([transtime_dn(i) transtime_dn(i+1) transtime_dn(i+1) transtime_dn(i)],[yl(1) yl(1) yl(2) yl(2)],colos(switchtonotegroup(i)+1,:),'edgecolor','none','facealpha',.1)
    text(transtime_dn(i),yl(1)+0.5*(yl(2)-yl(1)),num2str(switchtonotegroup(i)))

end
xl = get(gca,'xlim');
fill([transtime_dn(i) xl(2) xl(2) transtime_dn(i)],[yl(1) yl(1) yl(2) yl(2)],colos(switchtonotegroup(i)+1,:),'edgecolor','none','facealpha',.1)
% ylabel('p(B-J) / p(B-H|B-J)')
xlabel('time')
% plot(thissongtime,smooth_lv(propbj_dash,10),'kx')
set(gca,'xtick',transtime_dn(1:10:end),'xticklabel',datestr(transtime_dn(1:10:end)));
box off
axis tight
title('with dash')

disp('BJ')
grpstats(allprop_bj,notegroups_songs)
disp('BH')
grpstats(allprop_bh,notegroups_songs)
disp('BX')
grpstats(allprop_bx,notegroups_songs)



%% align song at light switches

switches = makeswitchstructure(notegroups_songs,transtime_t0,switchtonotegroup,day_song,time_song);


%% nOTE !!!! -  NEED TO GENERALIZE (LINK NG/COLOR WITH TARG TRANSITION)


figure
smoothie = 3;
nmin = 0;
subplot(2,2,1)
%transition to orange training from green training (target BJ)
thisix = find(switches.preng==3 & switches.postng==1);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix))>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),2)
line([0 0],get(gca,'ylim'),'color',lenacols(2,:));
box off
axis tight
hold on
%transition to orange probe from green training
thisix = find(switches.preng==3 & switches.postng==0);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix))>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),1)
line([0 0],get(gca,'ylim'),'color',lenacols(1,:));
%transition to green probe from green training
thisix = find(switches.preng==3 & switches.postng==2);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix),1)>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),3)
line([0 0],get(gca,'ylim'),'color',lenacols(3,:));
title('start in green training target BJ')

subplot(2,2,2)
%transition to orange training from green probe
thisix = find(switches.preng==2 & switches.postng==1);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix),1)>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),2)
line([0 0],get(gca,'ylim'),'color',lenacols(2,:));
box off
axis tight
hold on
%transition to orange probe from green probe
thisix = find(switches.preng==2 & switches.postng==0);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix),1)>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),1)
% line([presize presize],get(gca,'ylim'),'color',lenacols(1,:));
%transition to green training from green probe
thisix = find(switches.preng==2 & switches.postng==3);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix),1)>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),4)
% line([presize presize],get(gca,'ylim'),'color',lenacols(4,:));
title('start in green probe')
subplot(2,2,3)
%transition to green training from orange training
thisix = find(switches.preng==1 & switches.postng==3);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix))>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),2)
line([0 0],get(gca,'ylim'),'color',lenacols(2,:));
box off
axis tight
hold on
%transition to green probe from orange training
thisix = find(switches.preng==1 & switches.postng==2);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix))>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),1)
% line([presize presize],get(gca,'ylim'),'color',lenacols(1,:));
%transition to orange probe from orange training
thisix = find(switches.preng==1 & switches.postng==0);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix),1)>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),3)
% line([presize presize],get(gca,'ylim'),'color',lenacols(3,:));
title('start in orange training target BH')

subplot(2,2,4)
%transition to green probe from orange probe
thisix = find(switches.preng==0 & switches.postng==2);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix),1)>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),1)
line([0 0],get(gca,'ylim'));
box off
axis tight
hold on
%transition to orange training from orange probe
thisix = find(switches.preng==0 & switches.postng==1);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix),1)>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),4)
%transition to green training from orange probe
thisix = find(switches.preng==0 & switches.postng==3);
[bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smoothie);
takeix = sum(~isnan(bigmatrix),1)>nmin;
presize = presize-find(takeix,1,'first')+1;
shadeplot(1-presize:size(bigmatrix(:,takeix),2)-presize,bigmatrix(:,takeix),2)
%%
cd ..
cd ..