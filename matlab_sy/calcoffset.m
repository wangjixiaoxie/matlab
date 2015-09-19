function [offset]=calcoffset(templ, curdata)


       
        xcorrarray=xcorr(templ,curdata,'coeff');
        offset=find(xcorrarray==max(xcorrarray));


