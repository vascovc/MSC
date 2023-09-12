%% ex 3.1

clear all
close all
clc

tf=10000;
l=1;
t=0:tf;

%fixed jumps of l=1

r_fixed=zeros(length(t),2); %initialization the particle's positions 

x=0;
y=0;
for i=2:length(t)    % i is the number of steps from 1 to 10
    A=randi([1,4]);   % generates a random value drawn from the standard uniform distribution on the interval(0,1)
    mov_x=0;
    mov_y=0;
    if A==1
        mov_y=l;
    elseif A==2
        mov_x=l;
    elseif A==3
        mov_y=-l;
    elseif A==4
        mov_x=-l;
    end
    x=mov_x+x; 
    y=mov_y+y;
    r_fixed(i,:)=[x y]; % new position of the particle
end


%varying jumps

mu = [1.6, 2, 2.6];
lmax=1000;
r_var=zeros(length(t),2,3); %initialization the particle's positions 
l=zeros(length(t),2);
col=['r' 'b' 'k'];
for k=1:length(mu)
    x=0;
    y=0;
    for i=2:length(t)    % i is the number of steps from 1 to 10
        x_rand=rand(1);
        x_rand2=rand(1);
        l(i,k)=lmax/(((lmax^(mu(k)-1)-1)*x_rand+1)^(1/(mu(k)-1)));
        angle=2*pi*x_rand2;
        mov_x=l(i,k)*cos(angle);
        mov_y=l(i,k)*cos(angle);
        x=mov_x+x; 
        y=mov_y+y;
        r_var(i,:,k)=[x y]; % new position of the particle

    end
    figure(1)
    plot3(t,r_var(:,1,k),r_var(:,2,k),[col(k) '.'])
    hold on
    figure(2)
    plot(t,r_var(:,1,k),[col(k) '.'])
    hold on
    figure(3)
    plot(t,r_var(:,2,k),[col(k) '.'])
    hold on
end

figure(1)
plot3(t,r_fixed(:,1),r_fixed(:,2),'g.')
title('Positions vs Time')
legend('\mu = 1.6','\mu = 2','\mu = 2.6','l=1')
xlabel('t'),ylabel('x'),zlabel('y'), grid on
xlim([0 1200])
figure(2)
plot(t,r_fixed(:,1),'g.')
title('Positions vs Time (projection plane: Y=0)')
legend('\mu = 1.6','\mu = 2','\mu = 2.6','l=1')
xlabel('t'),ylabel('x'), grid on

figure(3)
plot(t,r_fixed(:,2),'g.')
title('Positions vs Time (projection plane: X=0)')
legend('\mu = 1.6','\mu = 2','\mu = 2.6','l=1')
xlabel('t'),ylabel('y'),grid on


% ex 3.2

l1=l(:,1);
l2=l(:,2);
l3=l(:,3);
c1=(mu(1)-1)/(1-lmax^(1-mu(1)));
c2=(mu(2)-1)/(1-lmax^(1-mu(2)));
c3=(mu(3)-1)/(1-lmax^(1-mu(3)));
Pl1=c1./(l1.^mu(1));
Pl2=c2./(l2.^mu(2));
Pl3=c3./(l3.^mu(3));
lt=0:1000;

figure
subplot(311)
h = hist(l1,100); 
h = h/sum(h); % normalize to unit length
bar(h); 
hold on
plot(l1,Pl1,'r.'), xlim([0 10]),ylim([0 1])
title('Probabilities, P(l) (\mu=1.6)')
legend('P_{experimental}','P_{levy}')
xlabel('l'),ylabel('Probabilities')

subplot(312)
h = hist(l2,100); 
h = h/sum(h); % normalize to unit length
bar(h); 
hold on
plot(l2,Pl2,'r.'), xlim([0 10]),ylim([0 1])
title('Probabilities, P(l) (\mu=2)')
legend('P_{experimental}','P_{levy}')
xlabel('l'),ylabel('Probabilities')

subplot(313)
h = hist(l3,100); 
h = h/sum(h); % normalize to unit length
bar(h); 
hold on
plot(l3,Pl3,'r.'), xlim([0 10]),ylim([0 1])
title('Probabilities, P(l) (\mu=2.6)')
legend('P_{experimental}','P_{levy}')
xlabel('l'),ylabel('Probabilities')

