% S5 APP6
% Laboratoire 1 - Problème 9
% Approximation par méthode « boîte grise » et racines de fonction non linéaire
%
% On reprend le Problème 3 sur MATLAB. La pression aérodynamique sur un avion se déplaçant à une
% vitesse v est donnée par P = 1/2 p v^2 où p est la densité de l’air qui varie de façon exponentielle avec
% l’altitude ℎ : p = p0 e^(-h/hs).
%
clc
close all
clear
clc

% Données du problème 9

P = [262.615 195.318 174.949 155.430 150.345 153.575 188.786 221.194 242.943 280.956 332.294]';
v = [20.000 20.628 22.553 25.894 30.862 37.769 47.048 59.284 75.244 95.931 122.646]';
h = [0 100 200 300 400 500 600 700 800 900 1000]';



% (a) La série de mesures {Pi, vi, ℎi} est fournie ci-dessous. Utiliser les équations développées au Problème 3
% pour identifier les valeurs des paramètres p0 et ℎs par la méthode des moindres carrés qui donne la
% meilleure approximation des mesures. Donner le coefficient R2 et l’erreur RMS de l’approximation
% linéaire. Donner l’erreur RMS dans la pression dynamique.  
  



  
  
% (b) En supposant que la vitesse v soit une fonction de l’altitude ℎ de la forme v = 20cosh[ℎ/400] et
% que p = 1.25e^(-h/300), utiliser la fonction non linéaire F(ℎ) = 0 développée au Problème 3 et calculer
% sur MATLAB avec la méthode de Newton-Raphson les altitudes ℎ auxquelles la pression aérodynamique
% est égale à une constante P0 = 200 N/m^2. Utiliser la dérivée de la fonction F(ℎ) par rapport à l’altitude
% ℎ développée au Problème 3.  
  










  
  