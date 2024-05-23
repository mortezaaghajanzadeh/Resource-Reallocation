clear all
clc
cd('E:\Dropbox\pollution_productivity_trade\replicationFiles');
addpath('E:\Dropbox\ado\matlab');
addpath('E:\Dropbox\pollution_productivity_trade\replicationFiles\codeMATLAB')
options = optimset('Display','off','MaxFunEvals',60000,'MaxIter',4500,'TolFun',1e-14,'TolX',1e-14,'Algorithm','trust-region-dogleg');
    



for mainpoll = 1:7 
    mainpoll
    clearvars -except options tol mainpoll;

    load 'dataMATLAB/rawFile.mat'
    
    run p2.m
    run p3.m
    
    if mainpoll == 2 run r1.m; end
                     run r2.m
    if mainpoll == 7 run r3.m; end

end;


