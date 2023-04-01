% S5 APP6
% Problematique
% Anthony Royer et Jérémy Goulet
% ROYA2019 - GOUJ
%
clc
close all
clear
clc
%
% Variables de contrôle
showGraphics = 1;           % Afficher des figures
showInitialGraphics = 1;    % Afficher les données initiales
showTerminalOutput = 1;     % Afficher les sorties dans le terminal
%
% Données de la problématique
APP6_S5_Prob_Constantes
load('Accelero_Data_from_NASA')

if and(showGraphics == 1, showInitialGraphics == 1)
    figure('Name','Données initiales')
    hold on
    plot(t,acc_mes,'r')
    xlabel('Temps (t)')
    ylabel('Accélération mesurée (m/s^2)')
    title('Données de la NASA')
    grid on
    hold off
end






