clc
clear all
close all

delta = 0.15;

q = 0.5;
start = 0;
tf = [400 401];

N = 5E4;

end_point_even = zeros(1,N);
end_point_odd = zeros(1,N);
for n=1:N
    end_point_odd(n) = random_walk(q,delta,start,tf(1));
    end_point_even(n) = random_walk(q,delta,start,tf(2));
end

[counts_even,number_evens] = groupcounts(end_point_even');
[counts_odd,number_odd] = groupcounts(end_point_odd');

counts_even = counts_even/N;
counts_odd = counts_odd/N;

events = [number_evens;number_odd];
count_events = 0.5*[counts_even;counts_odd];

plot(events,count_events,'.')
hold on
x = linspace(-tf(2),tf(2),5000);
plot(x,1./sqrt(2*pi.*mean(tf)).*exp(-(x-2*mean(tf)*delta).^2./(2.*mean(tf))))

ans_2_1 = sum(count_events);
disp(['2.1 ',num2str(ans_2_1)])

ans_2_2 = sum(events.*count_events);
disp(['2.2 ',num2str(ans_2_2)])
disp(['2.2_real ',num2str(2*tf(1)*delta )])

ans_2_3 = sum((events-ans_2_2).^2.*count_events);
disp(['2.3 ',num2str(ans_2_3)])

function start = random_walk(q,delta,start,tf)
for i=2:tf
if rand < q+delta
   add = 1;
else
   add = -1;
end
start = start + add;
end
end