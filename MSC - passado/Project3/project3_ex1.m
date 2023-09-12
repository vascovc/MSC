%% task 1.1

close all,clear all,clc

tf_list=[5000, 5001];
xx=-500:500;
yy=xx;
N=10000; %nr of repetitions

Nxy=zeros(length(xx),length(yy),2);
Pmean=zeros(length(xx),length(yy));
sum1=0;
sum2=0;
sum3=0;
sum4=0;

for n=1:N
    for j=1:length(tf_list)
        tf=tf_list(j);
        t=0:tf;
        x=0;
        y=0;
        for i=2:length(t)    % i is the number of steps from 1 to 10
            A=randi([1,4]);   % generates a random value drawn from the standard uniform distribution on the interval(0,1)
            mov_x=0;
            mov_y=0;
            if A==1
                mov_y=1;
            elseif A==2
                mov_x=1;
            elseif A==3
                mov_y=-1;
            elseif A==4
                mov_x=-1;
            end
            x=mov_x+x; 
            y=mov_y+y;
        end
        id_x=x+abs(min(xx));
        id_y=y+abs(min(yy));
        Nxy(id_x,id_y,j)=Nxy(id_x,id_y,j)+1;
    end
    
end

P_theor=zeros(length(xx),length(yy));
t_theor=mean(tf_list);
constant=1/(pi*t_theor);
for i=1:length(xx)
    for j=1:length(yy)
        P_theor(i,j)=constant*exp(-((xx(i)^2+yy(j)^2))/t_theor);
        Nmean=(Nxy(i,j,1)+Nxy(i,j,2))./2;
        if Nmean==0
            Pmean(i,j)=NaN;
        else
            Pmean(i,j)=Nmean./N;
            sum1=sum1+Pmean(i,j);
            sum2=sum2+Pmean(i,j)*xx(i);
            sum3=sum3+Pmean(i,j)*yy(i);
            sum4=sum4+Pmean(i,j)*(xx(i)^2+yy(i)^2);
        end
    end
end


[x,y]=meshgrid(xx,yy);
figure
subplot(221)
plot3(x,y,Pmean,'ko')
xlabel('x'),ylabel('y'),zlabel('P')
grid on
xlim([-400,400]),ylim([-400,400])
subplot(222)
mesh(x,y,P_theor)
xlabel('x'),ylabel('y'),zlabel('P')
grid on
xlim([-400,400]),ylim([-400,400])
subplot(223)
plot3(x,y,Pmean,'ko'),view(0,0)
xlabel('x'),ylabel('y'),zlabel('P')
grid on
xlim([-400,400]),ylim([-400,400])
subplot(224)
mesh(x,y,P_theor),view(0,0)
xlabel('x'),ylabel('y'),zlabel('P')
grid on
xlim([-400,400]),ylim([-400,400])
sgtitle('Probability distribution over (x,y): Experimental and Theoretical')

% checking some math properties of the generated random walks:
sum1
sum2
sum3
sum4

%% task 1.2

clear all, close all, clc

tf=100;
n=3; %number of random walks
t=0:tf;
r=zeros(length(t),2,n); %initialization the particle's positions 

col=['k','b','r'];

for k=1:n
    x=0; y=0;
    for i=2:length(t)    % i is the number of steps from 1 to 10
        A=randi([1,4]);   % generates a random value drawn from the standard uniform distribution on the interval(0,1)
        mov_x=0;
        mov_y=0;
        if A==1
            mov_y=1;
        elseif A==2
            mov_x=1;
        elseif A==3
            mov_y=-1;
        elseif A==4
            mov_x=-1;
        end
        x=mov_x+x; 
        y=mov_y+y;
        r(i,:,k)=[x y]; % new position of the particle
    end
    figure(1)
    plot3(t,r(:,1,k),r(:,2,k),[col(k) '-o'])
    hold on
    figure(2)
    plot(t,r(:,1,k),[col(k) '-o'])
    hold on
    figure(3)
    plot(t,r(:,2,k),[col(k) '-o'])
    hold on
end

figure(1)
title('Positions vs Time')
legend('random walk 1','random walk 2','random walk 3')
xlabel('t'),ylabel('x'),zlabel('y'), grid on

figure(2)
title('Positions vs Time (projection plane: Y=0)')
legend('random walk 1','random walk 2','random walk 3')
xlabel('t'),ylabel('x'), grid on

figure(3)
title('Positions vs Time (projection plane: X=0)')
legend('random walk 1','random walk 2','random walk 3')
xlabel('t'),ylabel('y'),grid on

