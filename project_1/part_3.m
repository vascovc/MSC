%% 3.1
clear all
close all
clc

n = [10 100 1000];
d_y = 0.005;

n_i = n(1);
y = nan([length(1:10E6) 1]);

for N=1:10E6
y(N) = 1/n_i*(sum(rand(n_i,1)));
end

[M,y_m,bins] = histcounts(y,(1/d_y)-1);

P_yk=nan(1,1/d_y-1);
y_k = nan(1,1/d_y-1);

for k = 1:1/d_y-1
    P_yk(k)=M(k)/N;
    y_k(k)=(k+0.5)*d_y;
end
normalize=sum(P_yk) % tem que dar 1

mean_y = mean(y);
var_y = var(y);
