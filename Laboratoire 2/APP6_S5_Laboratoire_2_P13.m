% S5 APP6
% Laboratoire 2 - Problème 13
% Calcul de la trajectoire de la capsule Apollo avec ODE45
%
% Faire l’exercice E28 (p. 52) du Chapitre 6 de Méthodes numériques (Notes_JdeL_ver10_rev1). Cet exercice
% utilise l’intégrateur ODE45 de MATLAB pour démontrer l’effet de la précision d’intégration sur le pas
% d’intégration et sur la qualité des résultats. Il sera utile d’utiliser « help ODE45 » sur MATLAB. Le
% « driver » de la partie (a) du Problème 12 ci-dessus donne un exemple d’appel à ODE45.
%
% E28 Développer le programme MATLAB tel que décrit ci-dessous et l’utiliser pour permettre d’atteindre les objectifs suivants :
% o Traduire un problème aux équations différentielles ordinaires en une fonction MATLAB
% o Initier à l'utilisation d'un intégrateur numérique à pas variable de type Runge-Kutta.
% o Vérifier l'effet de la tolérance d'intégration sur la solution numérique et sur le nombre de pas requis.
%
%  Utiliser l'intégrateur ODE45 sur MATLAB. Faire help ode45 sur MATLAB pour avoir l’information.
%  [T,Y] = ODE45('F',TSPAN,Y0,OPTIONS) avec TSPAN = [T0 TFINAL]
%  Intègre le système d'équations dy/dt = F(t,y) de T0 à TFINAL avec conditions initiales Y0.
%  'F' est le nom du fichier qui calcule les EDO d'ordre 1 à intégrer.
%  La fonction F(T,Y) doit retourner un vecteur colonne.
%  Chaque rangée de la solution Y correspond au temps dans la colonne T.
%  Sans l'argument "OPTIONS", les paramètres d'intégration nominaux sont utilisés.
%  Pour introduire de nouveaux paramètres, il faut les passer par la variable OPTIONS, qui est elle-même créée par la fonction ODESET (voir help odeset).
%  Pour changer la tolérance d'intégration relative ("RelTol") de 1.0e-03 (valeur par défaut) à 1.0e-06, on utilise: 
%    options = odeset('RelTol',1e-6);.
%
clc
close all
clear
clc

showGraphics = 1;

% Données du problème 13
currentFile = mfilename('fullpath');
[pathstr,~,~] = fileparts(currentFile);
addpath([pathstr,'\APP6e_PROBLEME_no_13']);

Driver_Apollo





