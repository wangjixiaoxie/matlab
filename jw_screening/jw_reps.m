function jw_reps(batch)
%format: jw_reps('batch') where batch is the batchfile of .not.mat files
%that contain the repeats to quantify.
%Modification of cara_reps. Calls ev_repeats internally.
%note that All_Lower argument of ev_repeats is 1

[tran_mat,reps,notes,nums]=ev_repeats(batch,1);

for repnumb = 3:length(notes)
vec = reps{repnumb};
mean=[1:length(vec)]*vec.'./sum(vec);
disp(['There are about ',num2str(mean),' reps per train for ',notes(repnumb)])
end
format compact
searchfor = input(['\nWhich syllable do you want to analyze? Enter one of: \n',notes, '\nYour syllable choice = '], 's');
foundsyl=findstr(searchfor, notes);
vec = reps{foundsyl};
mean=[1:length(vec)]*vec.'./sum(vec);
stdev=(sum(([1:length(vec)].^2).*(vec./sum(vec)))-((sum([1:length(vec)].*(vec./sum(vec)))).^2));
N=(sum([1:length(vec)]));
sterr=(stdev/sqrt(N));
disp(['Mean   = ',num2str(mean)])
disp(['StdDev = ',num2str(stdev)])
disp(['StdErr = ',num2str(sterr)])
%below is for plotting
figure;hold on;

bar(vec./sum(vec));
%from cara_reps: legend('Pre','1d post','wk4','wk6');
xlabel('Number of Repeats');ylabel('Probability');
