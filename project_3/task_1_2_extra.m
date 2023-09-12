clear all
close all
clc

% t_f = [50000 50001];
% t_f = [5000 5001];
% t_f = [500 501];

N = 1E5;
repeticoes = 50;
t_f_cells = [[50000 50001] [5000 5001] [500 501]];
for t_f_cell = 1:length(t_f_cells)/2
    t_f = [t_f_cells(t_f_cell*2-1) t_f_cells(t_f_cell*2)];
    normalization = nan(1,repeticoes);
    average_x = nan(1,repeticoes);
    average_y = nan(1,repeticoes);
    variance = nan(1,repeticoes);

    for num_reps = 1:repeticoes
        end_point_even = nan(N,2);
        end_point_odd  = nan(N,2);
        
        %tic
        parfor rep=1:N
            point = zeros(1,2);
            end_point_even(rep,:) = random_walk_2D(point,t_f(1));
            end_point_odd(rep,:)  = random_walk_2D(point,t_f(2));
        end
        %toc
        end_point = [end_point_even;end_point_odd];
        [counts, number] = groupcounts(end_point);
        
        counts = counts/N;
        
        events_x = number{1};
        events_y = number{2};
        count_events = 0.5*counts;
        normalization(num_reps) = sum(count_events);        
        average_x(num_reps) = sum(events_x.*count_events);
        average_y(num_reps) = sum(events_y.*count_events);
        variance(num_reps) = sum((events_x.^2+events_y.^2).*count_events);
    end
    disp(num2str(t_f(1)))
    disp(['Normalization mean: ',num2str(mean(normalization)),' +/- ',num2str(std(normalization))])
    disp(['x mean:             ',num2str(mean(average_x)),' +/- ',num2str(std(average_x))])
    disp(['y mean:             ',num2str(mean(average_y)),' +/- ',num2str(std(average_y))])
    disp(['variance mean:      ',num2str(mean(variance)),' +/- ',num2str(std(variance))])
end

function point = random_walk_2D(point,tf)
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
end
end 