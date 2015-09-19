function testout=testtrunc(test)
todel=[];
for rs=2:length(test)
	if test(rs)==test(rs-1)
		todel=[todel rs]
	end
end
test(todel)=[];
testout=test;
