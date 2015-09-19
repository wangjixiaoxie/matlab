function [pitchcurves,residuals]=FFanalysis725(batch,note,f_low,f_high)
[fvals]=findwnote4(batch,note,'','',0,[f_low f_high],6500,1,'obs0');
shifted=jc_Align(fvals);

pitchcurves=jc_PitchData729(shifted,1024,1020,1,f_low,f_high);
residuals=jc_plotresiduals725(pitchcurves);