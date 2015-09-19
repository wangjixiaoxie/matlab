%stability_analysis
%takes input from ~/ml/seqbirdmfiles/synbirdstruct
function [sqout,ss,sumout]=sqmst(ss)

%first find the indices which have stability analysis.
[sumout]=stabil_anal(ss)
PLOTFLAG=0
[outmat,ss]=stabil_plot(sumout,ss,PLOTFLAG)
[sqout]=smsq_anal(ss);

    