clear all
clc
figure
A = load('fort.16');
x = A(:,1);
y = A(:,2);

plot(x,y)

set(gca,'Xlim',[-1.5 1.5],'Ylim',[0 1])