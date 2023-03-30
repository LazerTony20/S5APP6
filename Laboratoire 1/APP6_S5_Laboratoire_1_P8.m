% S5 APP6
% Laboratoire 1 - Problème 8
% Racines d’une équation non linéaire avec Newton-Raphson
%
% Faire le problème E23 (p. 50) du Chapitre 5 de Méthodes numériques (Notes_JdeL_ver10_rev1) avec une
% tolérance de 1.0e-08. Revoir l’exercice E17 du Guide MATLAB pour la construction de boucle « while ».
% Dans ce problème, on remarque que dans le cas de racines multiples, le point de départ des itérations
% détermine vers quelles racines la solution va converger.
%
clc
close all
clear
clc

showGraphics = 1;

% Données du problème 8
% Avec la méthode de NR, trouver les racines de la fonction f(x) = x^3 − 6x^2 + 7x + 2 en démarrant les
% itérations aux valeurs initiales suivantes 
% (1) x = 0.85
% (2) x = 1.10
% (3) x = 0.709
% (4) x = 1.0
x_start = [0.85 1.10 0.709 1.0];

Tol = 1.0e-08;
fx = [1 -6 7 2];
dfx = [3 -12 7];
f = @(x) x.^3-6.*x.^2+7.*x+2;
d = @(x) 3.*x.^2-12.*x+7;
xlims = [-1 5];
ylims = [-15 15];
TooMuch = 500;

% x = 0.85
x = x_start(1);
F = f(x);
D = d(x);
iteration = 0;
while and(abs(F) > Tol,iteration <= TooMuch)
    x = x - F/D;
    F = f(x);
    D = d(x);
    iteration = iteration + 1;
    f_iterations(iteration) = F;
    x_iterations(iteration) = x;
end

if showGraphics == 1
    figure('Name',['Racine avec NR (tol = ', num2str(Tol), ') avec x = ', num2str(x_start(1))])
    hold on
    fplot(f,xlims,'Color',[.7 .7 .7])
    plot(x_start(1),f(x_start(1)),'ro')
    plot(x_iterations,f_iterations,'gp')
    plot(x,F,'bp')
    legend('f(x) = x^3 -6x^2 + 7x + 2', ['depart a ', num2str(x_start(1))],[num2str(iteration),' itérations'], ['Solution finale : ', num2str(x)],'Location','NorthWest')
    title(['Racine avec NR (tol = ', num2str(Tol), ')'])
    xlabel('x')
    ylabel('f(x)')
    axis([xlims ylims])
    grid on
    hold off
end

% x = 1.10
x = x_start(2);
F = f(x);
D = d(x);
iteration = 0;
while and(abs(F) > Tol,iteration <= TooMuch)
    x = x - F/D;
    F = f(x);
    D = d(x);
    iteration = iteration + 1;
    f_iterations(iteration) = F;
    x_iterations(iteration) = x;
end
if showGraphics == 1
    figure('Name',['Racine avec NR (tol = ', num2str(Tol), ') avec x = ', num2str(x_start(2))])
    hold on
    fplot(f,xlims,'Color',[.7 .7 .7])
    plot(x_start(2),f(x_start(2)),'ro')
    plot(x_iterations,f_iterations,'gp')
    plot(x,F,'bp')
    legend('f(x) = x^3 -6x^2 + 7x + 2', ['depart a ', num2str(x_start(2))],[num2str(iteration),' itérations'], ['Solution finale : ', num2str(x)],'Location','NorthWest')
    title(['Racine avec NR (tol = ', num2str(Tol), ')'])
    xlabel('x')
    ylabel('f(x)')
    axis([xlims ylims])
    grid on
    hold off
end


% x = 0.709
x = x_start(3);
F = f(x);
D = d(x);
iteration = 0;
while and(abs(F) > Tol,iteration <= TooMuch)
    x = x - F/D;
    F = f(x);
    D = d(x);
    iteration = iteration + 1;
    f_iterations(iteration) = F;
    x_iterations(iteration) = x;
end

if showGraphics == 1
    figure('Name',['Racine avec NR (tol = ', num2str(Tol), ') avec x = ', num2str(x_start(3))])
    hold on
    fplot(f,xlims,'Color',[.7 .7 .7])
    plot(x_start(3),f(x_start(3)),'ro')
    plot(x_iterations,f_iterations,'gp')
    plot(x,F,'bp')
    legend('f(x) = x^3 -6x^2 + 7x + 2', ['depart a ', num2str(x_start(3))],[num2str(iteration),' itérations'], ['Solution finale : ', num2str(x)],'Location','NorthWest')
    title(['Racine avec NR (tol = ', num2str(Tol), ')'])
    xlabel('x')
    ylabel('f(x)')
    axis([xlims ylims])
    grid on
    hold off
end

% x = 1.00
x = x_start(4);
F = f(x);
D = d(x);
iteration = 0;
while and(abs(F) > Tol,iteration <= TooMuch)
    x = x - F/D;
    F = f(x);
    D = d(x);
    iteration = iteration + 1;
    f_iterations(iteration) = F;
    x_iterations(iteration) = x;
end

if showGraphics == 1
    figure('Name',['Racine avec NR (tol = ', num2str(Tol), ') avec x = ', num2str(x_start(4))])
    hold on
    fplot(f,xlims,'Color',[.7 .7 .7])
    plot(x_start(4),f(x_start(4)),'ro')
    plot(x_iterations,f_iterations,'gp')
    plot(x,F,'bp')
    legend('f(x) = x^3 -6x^2 + 7x + 2', ['depart a ', num2str(x_start(4))],[num2str(iteration),' itérations'], ['Solution finale : ', num2str(x)],'Location','NorthWest')
    title(['Racine avec NR (tol = ', num2str(Tol), ')'])
    xlabel('x')
    ylabel('f(x)')
    axis([xlims ylims])
    grid on
    hold off
end

% Réponses :
% o x = 0.85 converge vers la racine à x = +4.2361 après 6 itérations
% o x = 1.10 converge vers la racine à x = +2.0000 après 4 itérations
% o x = 0.709 converge vers la racine à x = - 0.2361 après 32 itérations
% o x = 1.00 ne converge pas (arrêtée après 501 itérations): oscillations entre x = 1 et x = 3.