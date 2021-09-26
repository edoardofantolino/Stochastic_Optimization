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
    
    
R2 = [-2 0 0 0 0 0;
      -3.5 0 0 0 0 0;
      -6 0 0 0 0 0;
      -5 0 0 0 0 0;
      -7.5 0 0 0 0 0;
      -10 0 0 0 0 0];
    
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

vfqf = [261.9722
  257.1623
  254.6623
  255.6623
  253.1623
  250.6623];

outer_iter = 10;

vfrl = zeros(number_of_states, outer_iter);

for iter=1:outer_iter
    k = 1;
    J = zeros(number_of_states,1);
    Jn = zeros(number_of_states, 1);
    Q = zeros(number_of_states, number_of_actions);
    lambda = 0.995;
    A = 1500;
    B = 3000;

    max_number_of_iterations = 1e6;
    states = zeros(max_number_of_iterations,1);
    actions = zeros(max_number_of_iterations,1);

    states(1) = 1;
    actions(1) = 1;

    Q=zeros(number_of_states,2);

    for k=2:max_number_of_iterations

        alpha = A/(B+k);
        actions(k) = randi(number_of_actions);
        states(k) = find(mnrnd(1,P(states(k-1),:,actions(k))));

        r = R(states(k-1), states(k), actions(k));
        q = max(Q(states(k),:));
        Q(states(k-1), actions(k)) = (1-alpha)*Q(states(k-1),actions(k))...
                                        + alpha*(r+lambda*q);

    end

    for i=1:number_of_states
        [argvalue,argmax] = max(Q(i,:));
        Jn(i) = max(Q(i,:));
        mu(i) = argmax;
    end
    
    iter
    mu
    Jn;
    Q;

    [GC,GR] = groupcounts(states);
    
    vfrl(:,iter) = Jn;
end

for i=1:number_of_states
    avg = mean(vfrl(i,:));
    standdev = std(vfrl(i,:));
    x = linspace(avg-6*standdev,avg+6*standdev,100);
    y = normpdf(x,avg,standdev);
    figure(i)
    plot(x,y, 'k')
    hold on
    grid on
    xline(avg,'k')
    xline(vfqf(i), 'r')
end