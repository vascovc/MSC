% Vasco Costa - 97746

%% part 1
%% task 1.1
close all
clear all
clc
trajectories = 3;
num_jumps = 100;

position = zeros(2,num_jumps,trajectories);
for times=1:trajectories
    for jump=2:num_jumps+1
        dir = rand();
        mov_x = 0;
        mov_y = 0;
        if dir < 0.25
            mov_x = 1;
        elseif dir < 0.5
            mov_x = -1;
        elseif dir < 0.75
            mov_y = 1;
        else
            mov_y = -1;
        end
        position(1,jump,times) = position(1,jump-1,times) + mov_x;
        position(2,jump,times) = position(2,jump-1,times) - mov_y;
    end
end

figure
hold on
for i=1:trajectories
    plot(position(1,:,i),position(2,:,i))
end
%% task 1.1 version 2
close all
clear all
clc
trajectories = 3;
num_jumps = 100;

position = zeros(3,num_jumps+1,trajectories);
for times=1:trajectories
    for jump=2:num_jumps+1
        dir = rand();
        mov_x = 0;
        mov_y = 0;
        if dir < 0.25
            mov_x = 1;
        elseif dir < 0.5
            mov_x = -1;
        elseif dir < 0.75
            mov_y = 1;
        else
            mov_y = -1;
        end
        position(1,jump,times) = position(1,jump-1,times) + mov_x;
        position(2,jump,times) = position(2,jump-1,times) - mov_y;
        position(3,jump,times) = jump-1;
    end
end

figure
Legend = cell(trajectories,1);
for i=1:trajectories
    plot3(position(1,:,i),position(2,:,i),position(3,:,i),'.-')
    Legend{i} = strcat('Trajectory ',num2str(i));
    hold on
end
legend(Legend)
xlabel("X")
ylabel("Y")
zlabel("Jumps")
view(3)
saveas(gcf,"task_1_1_3d.png",'png')
view(2)
axis equal
saveas(gcf,"task_1_1_plano.png",'png')
%% task 1.2
clear all
close all
clc

tic
t_f = [50000 50001];
t_f = [5000 5001];
t_f = [500 501];

N = 1E5;
end_point_even = nan(N,2);
end_point_odd  = nan(N,2);

tic
parfor rep=1:N
    point = zeros(1,2);
    end_point_even(rep,:) = random_walk_2D(point,t_f(1));
    end_point_odd(rep,:)  = random_walk_2D(point,t_f(2));
end
toc
end_point = [end_point_even;end_point_odd];
[counts, number] = groupcounts(end_point);

counts = counts/N;

events_x = number{1};
events_y = number{2};
count_events = 0.5*counts;

plot3(events_x,events_y,count_events,'r.')
hold on
x_points = (-t_f(1):t_f(1))./10';
y_points = x_points;
[X,Y] = meshgrid(-2.5*sqrt(t_f(1)):1:2.5*sqrt(t_f(1)));
Z = 1/(pi*mean(t_f))*exp(-(X.^2 +Y.^2)/mean(t_f));
surface(X,Y,Z)

xlabel("X")
ylabel("Y")
zlabel("Probability")
legend('Experimental','Theory')
view(3)
saveas(gcf,['task_1_2_3d_',num2str(t_f(1)),'.png'],'png')
view([20 15])
saveas(gcf,['task_1_2_3d_2_',num2str(t_f(1)),'.png'],'png')
%contas
normalization = sum(count_events);
disp(['normalização: ',num2str(normalization)])
disp(['Esperado: ','1'])

average_x = sum(events_x.*count_events);
average_y = sum(events_y.*count_events);
disp(['average x: ',num2str(average_x)])
disp(['average y: ',num2str(average_y)])
disp(['Esperado: ','0'])

variance = sum((events_x.^2+events_y.^2).*count_events);
disp(['Varience: ',num2str(variance)])
disp(['Esperado: ',num2str(mean(t_f))])


%% task 2
%% task 2.1
clear all
close all
clc

x_b = -30;
%x_b = -50;
delta_t = 10;
N = 50E3;
%N = 50E4;
limit = N+1;

step = nan(N,1);
tic
parfor i=1:N
    point = zeros(1,2);
    step(i) = random_walk_2D_x_boundary(point,x_b,limit);
end
toc

[num_times, edges] = histcounts(step,'BinWidth',delta_t);

f = num_times/N;

s = nan(1,length(edges)-1);
abc = 1-sum(f);
for i=1:length(edges)-1
    s(i) = abc + sum(f(i:end));
end
x_0=0;
D=1/4;
t=edges(1:end-1)+delta_t/2;

S = erf(abs(x_b-x_0)./(2.*sqrt(D.*t)));
S_approx = abs(x_b-x_0)./sqrt(pi.*D.*t);

F = abs(x_b-x_0)./(2.*sqrt(pi.*D.*t.^3)).*exp(-(x_b-x_0).^2./(4.*D.*t));
F_approx = abs(x_b-x_0)./(2.*sqrt(pi.*D.*t.^3));

figure
plot(log(t),log(f),'.')
hold on
plot(log(t),log(F_approx))
plot(log(t),log(F))
xlabel("ln(t)")
ylabel("ln(F)")
legend("Experimental","Approximated","Theoretical")

saveas(gcf,['task_2_1_f_',num2str(x_b),'_',num2str(N),'.png'],'png')

figure
plot(log(t),s,'o')
hold on
plot(log(t),S_approx)
plot(log(t),S)
xlabel("ln(t)")
ylabel("ln(S)")
ylim([0 1.5])
legend("Experimental","Approximated","Theoretical")
saveas(gcf,['task_2_1_s_',num2str(x_b),'_',num2str(N),'.png'],'png')

figure
plot(t(1:100),log(f(1:100)),'.')
hold on
plot(t(1:100),log(F_approx(1:100)))
plot(t(1:100),log(F(1:100)))
xlabel("t")
ylabel("ln(F)")
legend("Experimental","Approximated","Theoretical")
[a,b] = max(F);
disp(num2str(t(b)))
%% part 3
%% part 3.1
clear all
close all
clc

mu = [1.6 2 2.6];
l_max = 1000;
N = 1000;
trajectories = length(mu);

position = zeros(3,N+1,trajectories);
for i=1:trajectories
    for j=2:N+1
        r = rand();
        l = (1-r*(1-l_max^(1-mu(i))))^(1/(1-mu(i)));
        phi = 2*pi*rand();
        position(1,j,i) = position(1,j-1,i)+l*cos(phi);
        position(2,j,i) = position(2,j-1,i)+l*sin(phi);
        position(3,j,i) = j-1;
    end
end

figure
Legend = cell(trajectories,1);
for i=1:trajectories
    plot3(position(1,:,i),position(2,:,i),position(3,:,i),'.-')
    Legend{i} = strcat('mu=',num2str(mu(i)));
    hold on
end
legend(Legend,'location','bestoutside')
xlabel("X")
ylabel("Y")
zlabel("Jumps")
view(3)
saveas(gcf,"task_3_1_3d.png",'png')
view(2)
saveas(gcf,"task_3_1_plano.png",'png')

position_random_walk = zeros(3,N+1,trajectories);
for i=1:trajectories
    for j=2:N+1
        dir = rand();
        mov_x = 0;
        mov_y = 0;
        if dir < 0.25
            mov_x = 1;
        elseif dir < 0.5
            mov_x = -1;
        elseif dir < 0.75
            mov_y = 1;
        else
            mov_y = -1;
        end
        position_random_walk(1,j,i) = position_random_walk(1,j-1,i) + mov_x;
        position_random_walk(2,j,i) = position_random_walk(2,j-1,i) - mov_y;
        position_random_walk(3,j,i) = j-1;
    end
end
figure
Legend = cell(2*trajectories,1);
counter=1;
for i=1:trajectories
    plot3(position(1,:,i),position(2,:,i),position(3,:,i),'.-')
    Legend{counter} = strcat('mu=',num2str(mu(i)));
    counter = counter+1;
    hold on
end
for i=1:trajectories
    plot3(position_random_walk(1,:,i),position_random_walk(2,:,i),position_random_walk(3,:,i),'.-')
    Legend{counter} = strcat('R',num2str(i));
    counter = counter+1;
    hold on
end
view(2)
legend(Legend,'location','bestoutside')
xlabel("X")
ylabel("Y")
zlabel("Jumps")
view(3)
saveas(gcf,"task_3_1_3d_r.png",'png')
view(2)
saveas(gcf,"task_3_1_plano_r.png",'png')

figure
Legend = cell(2*trajectories-1,1);
counter = 1;
for i=2:trajectories
    plot3(position(1,:,i),position(2,:,i),position(3,:,i),'.-')
    Legend{counter} = strcat('mu=',num2str(mu(i)));
    counter = counter+1;
    hold on
end
for i=1:trajectories
    plot3(position_random_walk(1,:,i),position_random_walk(2,:,i),position_random_walk(3,:,i),'.-')
    Legend{counter} = strcat('R',num2str(i));
    counter = counter+1;
    hold on
end
view(2)
legend(Legend,'location','bestoutside')
xlabel("X")
ylabel("Y")
zlabel("Jumps")
view(3)
saveas(gcf,"task_3_1_3d_r_2.png",'png')
view(2)
saveas(gcf,"task_3_1_plano_r_2.png",'png')

figure
Legend = cell(2*trajectories-2,1);
counter = 1;
for i=3:trajectories
    plot3(position(1,:,i),position(2,:,i),position(3,:,i),'.-')
    Legend{counter} = strcat('mu=',num2str(mu(i)));
    counter = counter+1;
    hold on
end
for i=1:trajectories
    plot3(position_random_walk(1,:,i),position_random_walk(2,:,i),position_random_walk(3,:,i),'.-')
    Legend{counter} = strcat('R',num2str(i));
    counter = counter+1;
    hold on
end
view(2)
legend(Legend,'location','bestoutside')
xlabel("X")
ylabel("Y")
zlabel("Jumps")
view(3)
saveas(gcf,"task_3_1_3d_r_3.png",'png')
axis equal
view(2)
saveas(gcf,"task_3_1_plano_r_3.png",'png')
%% task 3.2
clear all
close all
clc

mu = [1.6 2 2.6];
l_max = 1000;
N = 10E5;
num_mu = length(mu);

l = zeros(N+1,num_mu);
for i=1:num_mu
    for j=1:N+1
        r = rand();
        l(j,i) = (1-r*(1-l_max^(1-mu(i))))^(1/(1-mu(i)));
    end
end

l_theory = linspace(1,l_max,N+1);
C = zeros(1,num_mu);
P = zeros(N+1,num_mu);
for i=1:num_mu
    C(i) = (mu(i)-1)/(1-l_max^(1-mu(i)));
    P(:,i) = C(i).*l_theory.^(-mu(i));
end

for i=1:num_mu
figure
sgtitle(['\mu=', num2str(mu(i)), newline,'l_{max}=',num2str(l_max)])

%subplot(1,2,1)
h = histogram(l(:,i),'BinWidth',1,'Normalization','probability');
hold on
plot(l_theory,P(:,i),'r.')
xlabel('l')
ylabel('Probabilities')
xlim([0 10]),ylim([0 1])
legend('P_{experimental}','P_{Lévy}')

width=800;
height=400;
set(gcf,'position',[10,10,width,height])

saveas(gcf,['task_3_2_',num2str(mu(i)),'_',num2str(l_max),'.png'],'png')
end

%% task 3.3
clear all
close all
clc

mu = [1.6 2 2.6];
l_max = 1000;
N = 1E5;
for i=1:length(mu)
    position = zeros(N,2);
    c_mu=mu(i);
    parfor rep=1:N
        point = zeros(1,2);
        for j=1:N
            r = rand();
            l = (1-r*(1-l_max^(1-c_mu)))^(1/(1-c_mu));
            phi = 2*pi*rand();
            point(1) = point(1)+l*cos(phi);
            point(2) = point(2)+l*sin(phi);
        end
        position(rep,:)=round(point,-1);
    end
    [counts, number] = groupcounts(position);

    counts = counts/N;
    events_x = number{1};
    events_y = number{2};
    figure
    plot3(events_x,events_y,counts,'.')
    hold on
    mesh(zeros(2))
    xlabel('X')
    ylabel('Y')
    zlabel('Probability')
    axis tight
    view([260 10])
    saveas(gcf,['task_3_3_1_',num2str(mu(i)),'.png'],'png')

    figure
    plot3(events_x,events_y,counts,'.')
    hold on
    mesh(zeros(2))
    xlabel('X')
    ylabel('Y')
    zlabel('Probability')
    axis tight
    view([20 10])
    saveas(gcf,['task_3_3_2_',num2str(mu(i)),'.png'],'png')
end

%% end
close all
%% functions

%% function for task 1.1
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

%% function for task 2.1
function step = random_walk_2D_x_boundary(point,boundary,limit_tf)
step = nan;
for i=1:limit_tf
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
    step = i;
    return
end
end
end