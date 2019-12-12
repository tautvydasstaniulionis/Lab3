clear all
x = 0.1:1/22:1;

%d = (1 + 0.6*sin(2*pi*x/0.7)) + 0.3*sin(2*pi*x)/2;
y1 = (1 + 0.6 * sin (2 * pi * x / 0.7)) + 0.3 * sin (2 * pi * x) / 2;
plot(x,y1)
%------------------------------------------------------------------
Pos = []/10;
Hgt = [];
Wdt = []/100;
Gauss = [];

for n = 1:length(Pos)
    Gauss(n,:) =  Hgt(n)*exp(-((x - Pos(n))/Wdt(n)).^2);
end

PeakSig = sum(Gauss)+y1;
pks = findpeaks(y1);
btm = islocalmin(y1);
pks=pks/2;
[pks,locs,widths,proms] = findpeaks(PeakSig,x);
widths;

% r1 = 0.32 - 0.19;
% r2 = 0.87 - 0.71;
% % r1=pks(1)/2;
% % r2=pks(2)/2;
% c1=widths(1);
% c2=widths(2);


% c1=pks(1)-0.373/2;
% c2=pks(2)-0.373/2;
% r1=widths(1);
% r2=widths(2);


c1 = 0.19;
c2 = 0.87;
r1 = 0.17;
r2 = 0.15;



[p,n1] = findpeaks(y1);
for i = 1 : length(n1)
    c(i)= x(n1(i));
end

%------------------------------------------------------------------
plot(x,y1)

% 
% F1 = exp(-(x(i)-c(1))^2/(2*r1^2));
% F2 = exp(-(x(i)-c(2))^2/(2*r2^2));

e_avg=0;

w0 = randn(1);
w1 = randn(1);
w2 = randn(1);
e=0;
n=0.1;%mokymo zingsnis
y2 = zeros(1,20);
ee=0;
for a = 1:1000

%while e(i) ~= 0

    for i = 1:20

F1 = exp(-(x(i)-c(1))^2/(2*r1^2));
F2 = exp(-(x(i)-c(2))^2/(2*r2^2));
y2(i) = (F1*w1) + (F2*w2) + w0;
e(i)=y1(i)-y2(i);

%update parameters using current inputs ant current error
w0 = w0 + n*e(i)*1;
w1 = w1 + n*e(i)*F1;
w2 = w2 + n*e(i)*F2;

     end
%     ee=abs(e)
    
% e_avg=sum(e)/length(e) 
end

for i=1:20
F1 = exp(-(x(i)-c(1))^2/(2*r1^2));
F2 = exp(-(x(i)-c(2))^2/(2*r2^2));
y2(i) = (F1*w1) + (F2*w2) + w0;
e(i)=y1(i)-y2(i);


% w0 = w0 + n*e(i)*1;
% w1 = w1 + n*e(i)*F1;
% w2 = w2 + n*e(i)*F2;

end


%     
% 
% 
 plot(x,y1,x,y2)

