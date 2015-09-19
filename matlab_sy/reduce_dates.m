function [vlot]=reduce_dates(indays,firstday)
%takes dates from jesus format and puts them in a format relative to the
%first day listed

indays=indays-firstday;
indays=datevec(indays);
vlot=indays(:,3)+1;
vlot=ceil(vlot);
