function [dp, FRmatDist] = lt_neural_v2_PseudoPop_SameDiff(FRmatMotifByNeur, AllMotifRegexp, ...
    AllNeurLocation, thisloc, plotOn, distmetric)
%  OUTPUTS dprime (dist5rtibution of distances for diff type syls subtract
%  distrubitoin for same-type syls)

% distmetric = 'correlation';
% thisloc = 'LMAN';

temp_pu69_same = 0; % if 1, then for pu69, only calls same-type those with different following syls
% this allows me to ask whether the bump after syl offset is a lagged
% represntation of current syl, or predict subsequent syl


%% ========== CORRELATION BETWEEN ALL PAIRS 

% ============ 
inds = strcmp(AllNeurLocation, thisloc);

if ~any(inds)
    disp('THIS BRAIN AREA NO DATA');
    return
end

Y = pdist(FRmatMotifByNeur(:,inds), distmetric);
FRmatDist = squareform(Y);
nummotifs = size(FRmatMotifByNeur,1);

% ============ GET SAME-TYPE, DIFF-TYPE PAIRS
% --- get list of single syls (so can determine whether is same-type)
[~, ~, ~, AllMotifSinglesyl] = ...
    lt_neural_v2_extractSameType(AllMotifRegexp, {AllMotifRegexp{1}});

AllDist = [];
AllDistSametype = [];
for i=1:nummotifs
    for ii=i+1:nummotifs
        
        % ========= get distance
        distthis = FRmatDist(i, ii);
        
        % ========= is this same or diff pair of syls?
            issame = strcmp(AllMotifSinglesyl{i}, AllMotifSinglesyl{ii});

        if temp_pu69_same==1

            % --- these are cases where post syl is same, throw them out.
            % ....
            
            if (i==3 & ii==9) | (i==4 & ii==10) | (i==5 & ii==11)
                issame=2;
            end
                    
        end
        
        AllDist = [AllDist; distthis];
        AllDistSametype = [AllDistSametype; issame];
   
        if distthis<0.7 & issame==0
           disp([AllMotifRegexp{i} '-' AllMotifRegexp{ii}]); 
        end
            
    end
end

%% ================= output the dprime between distributions

% -- same
dat =  AllDist(AllDistSametype==1);

y1 = mean(dat);
y1var = var(dat);

% -- diff
dat =  AllDist(AllDistSametype==0);

y2 = mean(dat);
y2var = var(dat);

% -- dprime
dp = (y2 - y1)/sqrt(0.5*(y1var + y2var));


%% ================ PLOT
if plotOn==1
lt_figure; hold on;

x = [1 2];
Y = cell(1,2);

Y{1} = AllDist(AllDistSametype==1);
Y{2} = AllDist(AllDistSametype==0);

lt_plot_MultDist(Y, x, 1, 'k');
lt_plot_zeroline;
lt_plot_annotation(1, ['dprime(diff-same)=' num2str(dp)], 'r')
end



