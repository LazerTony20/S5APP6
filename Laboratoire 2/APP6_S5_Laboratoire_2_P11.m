% S5 APP6
% Laboratoire 2 - Problème 11
% Ordre de la méthode d’intégration et convergence
%
% Faire l’exercice E26 (p. 50) du Chapitre 6 de Méthodes numériques (Notes_JdeL_ver10_rev1). Cet exercice
% illustre la façon de déterminer à quel pas d’intégration un intégrateur numérique devient stable.
%
% E26 Cinq simulations numériques d’un problème dynamique ont été effectuées avec 5 pas d’intégration t
% différents. Les erreurs globales e correspondantes ont été mesurées :
%    e         t
% 0.08443    0.050
% 0.02603    0.040
% 0.01048    0.030
% 0.00319    0.020
% 0.00040    0.010
% À partir de ces données, déduire l’ordre de la méthode d’intégration numérique. Déduire aussi à partir de quel
% pas d'intégration la solution numérique est convergente.
%
clc
close all
clear
clc

showGraphics = 1;
showInitialData = 0;

% Données du problème 11
err =       [0.08443 0.02603 0.01048 0.00319 0.00040];
Delta_t =   [0.050   0.040   0.030   0.020   0.010];
longeur = length(Delta_t);

if showGraphics == 1
    if showInitialData == 1
        figure('Name','Données initiales')
        hold on
        plot(Delta_t,err)
        xlabel('\Deltat')
        ylabel('Erreur Globale')
        grid on
        hold off
    end
end

for m = 1:1:longeur
    p(m) = log(err(m)./err(end))./log(Delta_t(m)./Delta_t(end));
end

if showGraphics == 1
    figure('Name','Données P')
    hold on
    plot(Delta_t,p)
    xlabel('\Deltat')
    ylabel('P')
    title('Rapport d’erreur globale')
    grid on
    hold off
end















