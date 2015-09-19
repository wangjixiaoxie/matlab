function tvals=jc_IOtrigger(cntrng)
tvals=[];
for TH=5:-.1:4.4
	cntrng(1).TH=TH;
    for dist=0.3:0.1:0.8
        get_trigt2('batchTEMPtest',cntrng,0.3,128,1,1);
        [vals,trigs]=triglabel('batchTEMPtest','a',1,1,0,1);
        toff=[];
        for ii=1:length(trigs)
            toff=[toff;trigs(ii).toffset];
        end
        std(toff);

        tvals=[tvals;TH,sum(vals),std(toff)];
    end
end