% Fichier de constantes de la problématique

%% ANNEXE A : Notre mandat
% Paramètres
m_capsule = 50;             % En Kg
J_capsule = 1.5;            % Kg-m^2
R_Mars = 3397e03;           % En m
u_Mars = 42830e09;          % En m^3/s^2
S_aero_capsule = 0.80;      % En m^2
d_aero_capsule = 0.05;      % En m
C_do = 1.20;                % Coeff de trainée
C_Lalf = 0.80;              % Coeff de portance
C_Malf = -0.07;             % Coeff de couple
C_Mq = -0.05;               % Coeff amortissement
C_Mdelta = 0.10;            % Coeeff de volet aéro

% Conditions initiales à l’entrée atmosphérique
v_ini = 6100;               % En m/s
y_ini = deg2rad(-20.5);     % En degrés (converti en rad/s)
h_ini = 120000;             % En m
s_ini = deg2rad(0.0);       % En degrés (converti en rad/s)
Delta_ini = deg2rad(-80);   % En degrés (converti en rad/s)
q_ini = 0.0;                % En degrés/s

% Conditions finales désirées
v_fin = [250 300];          % En m/s (deux valeurs car deux options)
h_fin =  10000;             % En m


%% ANNEXE B : Résultats de l’expérience de la NASA
% Conditions initiales à l’entrée atmosphérique
v_ini_NASA = 6100;              % En m/s
y_ini_NASA = deg2rad(-90);      % En degrés (converti en rad/s)
h_ini_NASA = 120000;            % En m
s_ini_NASA = 0.0;               % En deg
Delta_ini_NASA = deg2rad(-90);  % En degrés (converti en rad/s)

% Bruit de mesure de l’accéléromètre
sigma_n_NASA = 0.035;
% Paramètre balistique
B_NASA = (C_do*S_aero_capsule)/m_capsule;


