clear all
close all
clc

%t_f = [50000 50001];
t_f = [5000 5001];
%t_f = [500 501];
%barrier = nan;

N = 1E5;
num_reps = 50;
for barrier = [-10,10,nan]
A_all = nan(num_reps,4);
for i=1:num_reps
    A_even = nan(N,1);
    A_odd = A_even;
    
    tic
    parfor rep=1:N
        point = zeros(1,2);
        A_even(rep) = random_walk_2D(point,t_f(1),barrier);
        A_odd(rep)  = random_walk_2D(point,t_f(2),barrier);
    end
    toc
    A = [A_even;A_odd];
    [counts, number] = groupcounts(A);
    for j=1:length(counts)
        if ~isnan(number(j))
            A_all(i,number(j)) = counts(j);
        end
    end
end

A_all_mean = mean(A_all);
error = std(A_all)./sum(A_all_mean);
X = categorical({'Bottom', 'Left', 'Top', 'Right'});
X = reordercats(X,{'Bottom', 'Left', 'Top', 'Right'});
Y = A_all_mean./sum(A_all_mean);
figure
bar(X,Y)
hold on
errorbar(X,Y,error,'k');

saveas(gcf,['task_2_1_extra_',num2str(barrier),'.png'],'png')
end
function A = random_walk_2D(point,tf,boundary)
for i=1:tf
    A = randi([1 4]);
    mov_x = 0; mov_y = 0;
    if A==1
        mov_y = 1;
    elseif A==2
        mov_x = 1;
    elseif A==3
        mov_y = -1;
    else
        mov_x = -1;
    end
    point(1) = point(1) + mov_x;
    point(2) = point(2) + mov_y;
    if point(1) == boundary
        A = nan;
        return
    end
end
if point ~= zeros(1,2)
    A = nan;
end
end

