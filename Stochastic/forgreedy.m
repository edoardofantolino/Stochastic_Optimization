close all
clear
clc

Bk = 0.5;
k = 200;

cnt_1 = 0;
cnt_0 = 0;
for i=0:1000
    p = binornd(1, 1-Bk/k);
    if p == 1
        cnt_1 = cnt_1 + 1;
    else
        cnt_0 = cnt_0 + 1;
    end
end

cnt_1
cnt_0