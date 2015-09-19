function [ output_args ] = untitled5(birdname, date_label)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

filename= [birdname '_' date_label '_entropy_b.mat'];

load ['/home/lucas/data/song/pu13bk43/all_days_entropy_randomWN/' filename]

figure; plot(data.b)


end

