function cara_reps(vec,vec1,vec4,vec6)

%vec is the vector from ev_repeats, for the specific repeat to analyze
%values are the rp{#} from ev_repeats for each week
%example p=find(pre_notes=='c');vec=pre_reps{p};
figure;hold on;

%plot and get statistics for pre
plot(vec./sum(vec),'b-','Linewidth',2);
mn=[1:length(vec)]*vec.'./sum(vec);mn
stdev=(sum(([1:length(vec)].^2).*(vec./sum(vec)))-((sum([1:length(vec)].*(vec./sum(vec)))).^2));stdev
N=(sum([1:length(vec)]));N
sterr=(stdev/sqrt(N));sterr

%plot and get statistics for 1d post
plot(vec1./sum(vec1),'g--','Linewidth',2);
mn1=[1:length(vec1)]*vec1.'./sum(vec1);mn1
stdev1=(sum(([1:length(vec1)].^2).*(vec1./sum(vec1)))-((sum([1:length(vec1)].*(vec1./sum(vec1)))).^2));stdev1
N1=(sum([1:length(vec1)]));N1
sterr1=(stdev1/sqrt(N1));sterr1

%plot and get statistics for week4
plot(vec4./sum(vec4),'r--','Linewidth',2);
mn4=[1:length(vec4)]*vec4.'./sum(vec4);mn4
stdev4=(sum(([1:length(vec4)].^2).*(vec4./sum(vec4)))-((sum([1:length(vec4)].*(vec4./sum(vec4)))).^2));stdev4
N4=(sum([1:length(vec4)]));N4
sterr4=(stdev4/sqrt(N4));sterr4

%plot and get statistics for week6
plot(vec6./sum(vec6),'k--','Linewidth',2);
mn6=[1:length(vec6)]*vec6.'./sum(vec6);mn6
stdev6=(sum(([1:length(vec6)].^2).*(vec6./sum(vec6)))-((sum([1:length(vec6)].*(vec6./sum(vec6)))).^2));stdev6
N6=(sum([1:length(vec6)]));N6
sterr6=(stdev6/sqrt(N6));sterr6

legend('Pre','1d post','wk4','wk6');xlabel('Number of Repeats');ylabel('Probability');