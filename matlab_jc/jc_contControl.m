function jc_contControl(ZFnorm,BFnorm,window)

for i=0:3
    sds=0+0.5*i;
    jc_contAlldataAV(ZFnorm,BFnorm,sds,'below',window);
    jc_contAlldataAV(ZFnorm,BFnorm,sds,'above',window);
end