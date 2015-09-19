function jc_520(ab)
tvals=[];
for TH=2.7:-.2:1.2
	cntrng(1).TH=TH;
	get_trigt('batchFinal.rand',ab,0.08,64,1,1);
	[vals]=jc_triglabel('batchFinal.rand','a',1,1,0,1);
	tvals=[tvals;TH,sum(vals)];
end

plot(tvals(:,1),tvals(:,2)./tvals(:,3),'bs-')