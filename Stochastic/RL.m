close all
clear
clc

p12 = 1/15;
p11 = 1-p12;
p23 = 1/15;
p22 = 1-p23;

P1 = [(1-2*(p12*p11)-(p12*p12)) (p12*p11)+(p12*p11) 0 (p12*p12) 0 0 ;
        0 (1-2*(p12*p11)-(p12*p12)) (p23*p11) (p22*p12) (p12*p23) 0 ;
        0 0 (p11) 0 (p12) 0;
        0 0 0 (1-2*(p22*p23)-(p23*p23)) (p23*p22)+(p23*p22) (p23*p23);
        0 0 0 0 (p22) (p23);
        1 0 0 0 0 0];

P2 = [1 0 0 0 0 0;
      1 0 0 0 0 0;
      1 0 0 0 0 0;
      1 0 0 0 0 0;
      1 0 0 0 0 0;
      1 0 0 0 0 0];

R1 = [2 1.5 0 1 0 0;
      0 1.5 1 1 0.5 0;
      0 0 1 0 0.5 0;
      0 0 0 1 0.5 0;
      0 0 0 0 0.5 0;
      -10 0 0 0 0 0];
    
    
ct = -2;
rw = -1.5;
rb = -4;

R2 = [ct 0 0 0 0 0;
      ct+rw 0 0 0 0 0;
      ct+rb 0 0 0 0 0;
      ct+rw+rw 0 0 0 0 0;
      ct+rw+rb 0 0 0 0 0;
      ct+rb+rb 0 0 0 0 0];
    
P = P1;
P(:,:,2) = P2;

R = R1;
R(:,:,2) = R2;
    
P = P1;
P(:,:,2) = P2;

R = R1;
R(:,:,2) = R2;



number_of_actions = 2;
number_of_states = length(P);
k = 1;
J = zeros(number_of_states,1);
Jn = zeros(number_of_states, 1);
Q = zeros(number_of_states, number_of_actions);
lambda = 0.995;
A = 1500;
B = 3000;

max_number_of_iterations = 6;
states = zeros(max_number_of_iterations,1);
actions = zeros(max_number_of_iterations,1);

states(1) = 1;
actions(1) = 1;

Q=zeros(number_of_states,2);
tic;
for k=2:max_number_of_iterations
    
    alpha = A/(B+k);
    actions(k) = randi(number_of_actions);
    states(k) = find(mnrnd(1,P(states(k-1),:,actions(k))));
    
    r = R(states(k-1), states(k), actions(k));
    q = max(Q(states(k),:));
    Q(states(k-1), actions(k)) = (1-alpha)*Q(states(k-1),actions(k))...
                                    + alpha*(r+lambda*q);
    
%     plot(max(Q'))
%     hold on
end
toc;

for i=1:number_of_states
    [argvalue,argmax] = max(Q(i,:));
    Jn(i) = max(Q(i,:));
    mu(i) = argmax;
end

mu
Jn
Q

[GC,GR] = groupcounts(states);
