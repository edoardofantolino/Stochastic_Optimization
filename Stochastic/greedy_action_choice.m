close all
clear
clc

iter=200;
Bk=0.5;
k=2:iter;
y=zeros(iter,1);
y=1-Bk./k;
ymy=1-Bk./log(k);

plot(k,y)
grid on
hold on
plot(k,ymy)