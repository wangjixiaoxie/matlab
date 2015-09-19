guide to inactivation analysis

DATA PREPARATION
*****************************
*****************************

basic analysis.
First each bird has a dat file, which initializes variables.

Then PITCHSYNTAXANAL<X>(AVLS)  does basic pitch analysis and organization across files.
    pitchsyntaxanal2 is standard, pitchsyntaxanal3 has added capacity to realign notes according to template.
load pathvals<X>-analdata.mat
INACTIVANAL3(AVLS,GRAPHVALS) does additional plotting.

********************************
META-ANALYSIS OF PHARMACOLOGY DATA.
birdstructlist.m - list of birds
shiftanal11 - latest version, outputs phsumbs struct
mkdynstruct.m - adds additional functionality.

shf

[shsall,sumshs,sumbs]=shiftanal8(bs,1:9)
[sumplot]=selectinbasruns5(sumbs,1 ,1:5 ,1.5)
[sumdyn]=selectdynamicruns2- find the diminishment runs.

PLOTDYNMETASCRIPT CALLS BELOW
plotdynamics3(sumdyn).



[stats]=calcsumplotstats(sumplot)

2.  %baseline variability effect

INACTIVEXPU34BAS, plots raw data for pu34.

PLOTBAS - plotsbaseline scatter
function [sumstats]=plotscatterbas(sumbs);
%this creates the scatter plot for baseline figure


3. %Example figure - bk61w42
INACTIVEXBK61W42.


%this creates the summary figures for plotbarinitind3
PLOTBARINITIND3(sumplot,'bar')
plotdynamics3(sumdyn)...
plotscatteffpctchange

PLOTSCATTERXY3


INACTIV_DYNAMICS, plots time course.

contouranal - plots an example contour for bk61w42


%STIMANAL added 6.27, modified 9.07
standard analysis.  
set parameters, i.e. contouranalbk48w74b.m
contouranal(avls)
stimanal(avls)
stimanalvardelay
plotstimcontour2
plotstimexample.m
metstimanal takes sumstructlist struct to make sumbs.
plotbarstim-takes sumbs, to make basic bar plot.


contourana2l.m - makes contours.
stimanalscript2 - does analysis of contours
plotstimcrountour.m - plots contours
variabledelaystimanal.m -- analysis of variable stim exps.

****getresiduals2*** takes sps (short stim struct set up in stimshortstructlist.m)





%FIGURE GUIDE..
%MASTERSCRIPTS.



PLOTHISTORYFIGSCRIPT.M

