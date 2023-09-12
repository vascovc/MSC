%% general
clear 
close all
clc

N = 10000;
m = 100;
for c=[50,100,500]
%c = 100;
L = N/2*c;

q_all = nan(m,N);
q_mean_all = nan(m,1);
B_all = nan(m,1);
show_that_all = nan(m,1);
q_mean_2_all = nan(m,1);
q_mean_3_all = nan(m,1);

number_ocurrences_all = nan(m,N);

C=nan(m,1);
C_approx = nan(m,1);
rho = nan(m,1);
parfor iter=1:m
    
    %% q normal
    A = create_adjacency_matrix(N,L);
    q = sum(A)';
    q_mean = 1/N.*sum(q);
    B = 1/N/q_mean.*sum(q.*(q-1));
    show_that = B/q_mean;
    q_mean_2 = 1/N*sum(q.^2);
    q_mean_3 = 1/N*sum(q.^3);
    
    % [number_ocurrences,value] = groupcounts(q);
    % P = number_ocurrences/N;
    q_all(iter,:) = q;
    q_mean_all(iter) = q_mean;
    B_all(iter) = B;
    show_that_all(iter) = show_that;
    q_mean_2_all(iter) = q_mean_2;
    q_mean_3_all(iter) = q_mean_3;
    
    %% q probabilidade
    A = create_adjacency_matrix(N,L);
    q = sum(A)';
    q_mean = 1/N.*sum(q);
    B = 1/N/q_mean.*sum(q.*(q-1));
    show_that = B/q_mean;
    q_mean_2 = 1/N*sum(q.^2);
    q_mean_3 = 1/N*sum(q.^3);
    
    [number_ocurrences,value] = groupcounts(q);
    number_ocurrences = number_ocurrences/N;
    line = nan(1,N);
    line(value) = number_ocurrences;
    number_ocurrences_all(iter,:) = line;
    
    B_all(iter) = B;
    show_that_all(iter) = show_that;
    q_mean_2_all(iter) = q_mean_2;
    q_mean_3_all(iter) = q_mean_3;
    
    %%task 2
    %ntr = 0;
    
%     for i=1:N
%         for j=1:N
%             if A(i,j)
%                 for k=1:N
%                     if A(j,k)==1 && A(k,i)==1
%                         ntr = ntr+1;
%                     end
%                 end
%                 numerator_task_3 = numerator_task_3+A(i,j)*(q(i)-Q_task_3)*(q(j)-Q_task_3);
%             end
%         end
%     end
    ntr = trace(A^3);
    ntr = ntr*1/6;
    B_times_q_mean = 1/N.*sum(q.*(q-1));
    npt = 1/6*N*B_times_q_mean;
    
    C(iter) = ntr/npt;
    C_approx(iter) = q_mean/N;
    
    %%task 3
    numerator_task_3 = 0;
    Q = q_mean_2/q_mean;
    for i=1:N
        for j=1:N
            if A(i,j) == 1
                numerator_task_3 = numerator_task_3+(q(i)-Q)*(q(j)-Q);
            end
        end
    end
    sigma_square = q_mean_3/q_mean-(q_mean_2/q_mean).^2;

    denominator = N*q_mean*sigma_square;
    
    rho(iter) = numerator_task_3/denominator;

end
%% task 1 - prob
P = mean(number_ocurrences_all,1);

nonNaN_indices = find(~isnan(P));
first_index = nonNaN_indices(1);
last_index = nonNaN_indices(end);

x_values = max(1,first_index-20):1:last_index+20;
P_theory = exp(-c)*c.^x_values./factorial(x_values);

figure
plot(x_values,P_theory)
hold on
plot(x_values,P(x_values),'o')
xlabel("q")
ylabel("P(q)")
legend("P_{theory}","P_{obtained}")
saveas(gcf,['task_1_',num2str(N),'_',num2str(c),'.png'],'png')

disp('\n\nTASK 1 - Using P')
disp(['c = ',num2str(c)])
disp(['<q> = ',num2str(mean(q_mean_all))])
disp(['std = ',num2str(std(q_mean_all))])
disp(['B = ',num2str(mean(B_all))])
disp(['std = ',num2str(std(B_all))])
disp(['<q^2> = ',num2str(mean(q_mean_2_all))])
disp(['std = ',num2str(std(q_mean_2_all))])
disp(['<q^3> = ',num2str(mean(q_mean_3_all))])
disp(['std = ',num2str(std(q_mean_3_all))])
disp(['B/<q> = ',num2str(mean(show_that_all))])
disp(['std = ', num2str(std(show_that_all))])

%% task 2 analysis
disp('\n\nTASK 2')
disp(['C = ',num2str(mean(C))])
disp(['std = ', num2str(std(C))])
disp(['C_approx = ',num2str(mean(C_approx))])
disp(['std = ', num2str(std(C_approx))])

figure
plot(1:m,C,'.')
hold on
plot(1:m,C_approx,'.')
legend("C","C_{approx}")
xlabel('Number of repetition')
saveas(gcf,['task_2_',num2str(N),'_',num2str(c),'.png'],'png')


%% task 3 analysis
disp('\n\nTASK 3')
disp(['rho = ',num2str(mean(rho))])
disp(['std = ', num2str(std(rho))])

figure
plot(1:m,rho,'.')
hold on
plot(1:m,zeros(1,m),'k-.')
legend('\rho','0')
xlabel('Number of repetition')
ylabel('Pearson coefficient')
saveas(gcf,['task_3_',num2str(N),'_',num2str(c),'.png'],'png')
end