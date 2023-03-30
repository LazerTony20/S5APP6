% S5 APP6
% Procédural 1
clc
close all
clear
clc

%% Probleme no 2
disp('Probleme 2')
xn = [1 2 3 4 5];
yn = [11.4 6.1 3.6 3.3 2.1];
% X = 1/x
% Y = y
% m = beta
% b = alpha
% Xn = [1 0.5 1/3 1/4 1/5];
Xn = 1./xn;
Yn = yn;

sumYn = sum(Yn);
N = sum(Xn.^0);
sumXn1 = sum(Xn.^1);
sumXn2 = sum(Xn.^2);
sumYnXn = sum(Xn.*Yn);

Ymat = [sumYn;
        sumYnXn];

A = [N sumXn1;
    sumXn1 sumXn2];
Ainv = inv(A);

Params = Ainv*Ymat;
b = Params(1);
m = Params(2);
disp(['b = ', num2str(b)])
disp(['m = ', num2str(m)])

Ymoy = (1/N).*sum(Yn);
% (A)
disp('(A)')
disp(' ')
% (B)
disp('(B)')
disp(' ')
% (C)
disp('(C)')
ynChapeau = b + m./xn;
YnChapeau = m.*Xn + b;

R2 = (sum((YnChapeau-Ymoy).^2))./(sum((yn-Ymoy).^2));
RMSY = sqrt((1/N).*sum((YnChapeau-Yn).^2));
RMSy = sqrt((1/N).*sum((ynChapeau-yn).^2));

disp(['R2 = ', num2str(R2)])
disp(['RMSY = ', num2str(RMSY)])
disp(['RMSy = ', num2str(RMSy)])

%% Problème 3 : Approximation par méthode « boîte grise » et racines de fonction non linéaire
% La pression aérodynamique sur un avion se déplaçant à une vitesse v








