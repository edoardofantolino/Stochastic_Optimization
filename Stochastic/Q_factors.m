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
   
sizes = size(P);
number_of_actions = sizes(3);
number_of_states = sizes(1);
k = 1;
J = zeros(number_of_states,1);
Jn = zeros(number_of_states, 1);
Q = zeros(number_of_states, number_of_actions);
epsilon = 0.001;
sp = 2*epsilon;
lambda = 0.995;

while sp > (((epsilon*(1-lambda))/(2*lambda)))
% for i=1:1e4
    for z=1:number_of_states

        partial = zeros(number_of_states, 1);
        for j=1:number_of_actions
            rbar = P(z,:,j)*R(z,:,j)';
            summation = lambda * sum(P(z,:,j)*J(:));

            partial(j) = rbar + summation;
            Q(z,j) = partial(j);          
        end
        Jn(z) = max(partial);

    end

%     Jn
%     figure(1)
%     plot(Jn)
%     hold on
%     grid on
%     sp = max(Jn-J) - min(Jn-J);
    sp = norm(Jn-J);
    J = Jn;
end


for i=1:number_of_states
    [argvalue,argmax] = max(Q(i,:));
    mu(i) = argmax;
end

mu
Jn
Q
