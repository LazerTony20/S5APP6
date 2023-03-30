% S5 APP6
% Laboratoire 1 - Problème 6
% Intégration numérique avec méthodes du trapèze et de Simpson
%
% Dans le chapitre 4 de Méthodes numériques (Notes_JdeL_ver10_rev1) :
% Faire E17 (p. 41) : intégration numérique avec la méthode du trapèze
% Faire E18 (p. 43) : intégration numérique avec la méthode de Simpson.
% Comparer la précision de chaque méthode et l'effet de la résolution des données.
%
clc
close all
clear
clc

showGraphics = 0;

% Données du problème 6 (exercices E17, E18)
x = [-2 -1.5 -1.0 -0.5  0.0 0.5  1.0 1.5   2.0  2.5   3.0  3.5   4.0]';
F = [14  8.75 5.0  2.75 2.0 2.75 5.0 8.75 14.0 20.75 29.0 38.75 50.00]';

% Affichage des données initiales
if showGraphics == 1
    figure('Name','Données initiales')
    hold on
    plot(x,F)
    xlabel('Distance (m)')
    ylabel('Force (N)')
    legend('F(x)','Location','NorthWest')
    title('Force appliquée sur un véhicule selon sa position')
    grid on
    hold off
end
%% Intégration numérique avec la méthode du trapèze
disp('==========Intégration numérique avec la méthode du trapèze==========')
% E17 - La force appliquée F sur un véhicule en fonction de sa position x a été mesurée à des intervalles de 0.5 m.
% (a) Avec la méthode trapézoïdale, calculer le travail effectué par cette force sur le véhicule.
disp('(a) Avec la méthode trapézoïdale, calculer le travail effectué par cette force sur le véhicule.')
ha = 0.5;        % Intervalle
Nm = length(x);

InteTrapeze = (F(1) + F(Nm) + 2*sum(F(2:Nm-1)))*ha/2;

delx = ha;
fpa = (F(2) - F(1))/ delx;
fpb = (F(end) - F(end-1))/delx;
ErreurTrapeze = ((ha^2)/12)*(fpb - fpa);

disp(['Le travail effectué par la force sur le véhicule est de : ', num2str(InteTrapeze), ' Joules'])
disp(['L`erreur d`intégration est approximativement : ', num2str(ErreurTrapeze), ' Joules'])
disp(' ')

% (b) Faire le même calcul si l’intervalle entre les mesures avait été de 1.0 m plutôt que de 0.5 m. 
% Dans les deux cas, calculer approximativement l’erreur d’intégration en calculant approximativement 
% les dérivées premières f′(a) et f′(b) aux extrémités de l’intervalle.
disp('(b) Faire le même calcul si l’intervalle entre les mesures avait été de 1.0 m plutôt que de 0.5 m.')
hb = 1;        % Intervalle
som = 0;
for m = 3:(hb/ha):Nm-2
    som = som + F(m);
end

InteTrapeze2 = (F(1) + F(Nm) + 2*som)*hb/2;

delx2 = hb;
fpa2 = (F(3) - F(1))/ delx2;
fpb2 = (F(end) - F(end-2))/delx2;
ErreurTrapeze2 = ((hb^2)/12)*(fpb2 - fpa2);

disp(['Le travail effectué par la force sur le véhicule est de : ', num2str(InteTrapeze2), ' Joules'])
disp(['L`erreur d`intégration est approximativement : ', num2str(ErreurTrapeze2), ' Joules'])
disp(' ')

% Note : En intégrant la fonction analytique qui a généré ces données, on obtient la réponse exacte de 84 Joules. Pour
% l’intervalle de 0.5, l’erreur prédite est 0.6875 et l’erreur réelle est 0.7500. Pour l’intervalle de 1.0, l’erreur prédite est
% 2.50 et l’erreur réelle est 3.00.

%% Intégration numérique avec la méthode de Simpson.
disp('==========Intégration numérique avec la méthode de Simpson==========')
h3 = ha;

InteSimp1 = (F(1) + F(Nm) + 4*sum(F(2:2:Nm -1)) + 2*sum(F(3:2:Nm-1)))*h3/3;

delx3 = h3;
fpppa = (F(4) - 3*F(3) + 3*F(2) - F(1))/delx3^3;
fpppb = (F(end) - 3*F(end-1) + 3*F(end-2) - F(end-3) )/delx3^3;
ErreurSimpson1 = ((h3^4)/180)*(fpppb - fpppa);

disp(['Le travail effectué par la force sur le véhicule est de : ', num2str(InteSimp1), ' Joules'])
disp(['L`erreur d`intégration est approximativement : ', num2str(ErreurSimpson1), ' Joules'])
disp(' ')

%% Comparer la précision de chaque méthode et l'effet de la résolution des données.  