% S5 APP6
% Problematique
% Anthony Royer et Jérémy Goulet
% ROYA2019 - GOUJ2711
%
clc
close all
clear
clc
opengl software
%
% Variables de contrôle
showGraphics = 1;           % Afficher des figures
showInitialGraphics = 0;    % Afficher les données initiales
showTerminalOutput = 1;     % Afficher les sorties dans le terminal
%
% Données de la problématique
APP6_S5_Prob_Constantes
load('Accelero_Data_from_NASA')

if and(showGraphics == 1, showInitialGraphics == 1)
    figure('Name','Données initiales')
    hold on
    plot(t,-acc_mes,'bo')
    xlabel('Temps (s)')
    ylabel('Accélération (Daero/masse) en m/s^2')
    title('Mesures accélérométriques de la NASA avec \gamma = -90 deg')
    grid on
    hold off
end

%% Identification par intégration numérique et par lissage
% Intégration numérique avec la méthode du trapèze
if showTerminalOutput == 1
    disp('==================================================================================================')
    disp('======================Identification par intégration numérique et par lissage=====================')
    disp('==================================================================================================')
end

ha = t(2)-t(1);        % Intervalle
Nm = length(acc_mes);

for n = 1:Nm
    vit_mes(n,1) = v_ini - (acc_mes(1) + acc_mes(n)+2*sum(acc_mes(2:n-1)))*ha/2;
end
% Calcul de l'erreur
fpa_vit = (acc_mes(2)-acc_mes(1))/ha;
fpb_vit = (acc_mes(end)-acc_mes(end-1))/ha;
err_vit_trap = (ha^2/12)*(fpb_vit-fpa_vit);

if showTerminalOutput == 1
    disp('==========Intégration numérique des données accélérométriques avec la méthode du trapèze==========')
    disp(['L`erreur d`intégration de la vitesse : ', num2str(err_vit_trap), ' m/s'])
    disp(' ')
end

if showGraphics == 1
    figure('Name','Données trapèze (Vitesse)')
    hold on
    plot(t,vit_mes,'ro')
    xlabel('Temps (s)')
    ylabel('Vitesse (Daero/masse) en m/s')
    title('Mesures de vitesse de la NASA avec \gamma = -90 deg avec trapèze')
    grid on
    hold off
end

% Intégration numérique avec la méthode de Simpson
alt_mes(1,1) = h_ini;
for n = 3:2:Nm
    alt_mes(((n-1)/2)+1,1) = h_ini - ((vit_mes(1) + vit_mes(n) + 4*sum(vit_mes(2:2:n-1)) + 2*sum(vit_mes(3:2:n-1)))*ha/3);
    t_simp(((n-1)./2)+1,1) = t(n);
end

% Calcul de l'erreur
fpppa_simp = (vit_mes(4) - 3*vit_mes(3) + 3*vit_mes(2) - vit_mes(1))./ha^3;
fpppb_simp = (vit_mes(end) - 3*vit_mes(end-1) + 3*vit_mes(end-2) - vit_mes(end-3))./ha^3;
err_alt_simp = ((ha^4)/180)*(fpppb_simp - fpppa_simp);

if showTerminalOutput == 1
    disp('==========Intégration numérique de la vitesse avec la méthode de Simpson==========')
    disp(['L`erreur d`intégration de l`altitude : ', num2str(err_alt_simp), ' m'])
    disp(' ')
end

if showGraphics == 1
    figure('Name','Données Simpson (Altitude)')
    hold on
    plot(t_simp,alt_mes,'go')
    xlabel('Temps (s)')
    ylabel('Altiude (Daero/masse) en m')
    title('Mesures d`altitude de la NASA avec \gamma = -90 deg avec Simpson')
    grid on
    hold off
end

% Équations dynamiques ainsi que Gravité et aérodynamique
r = R_Mars + alt_mes; 
gr = u_Mars./r.^2;
D_aero = (acc_mes(1:2:end).*m_capsule);
P_dyn = D_aero./(S_aero_capsule.*C_do);
D_aero_bruite = -(acc_mes(1:2:end) + gr).*m_capsule;
P_dyn_bruite = D_aero_bruite./(S_aero_capsule.*C_do);

% Approximation à deux paramètres
% forme : Y = mX + b
% Depuis ρ = ρ0*exp(-h/hs)
% Avec Pdyn = (ρ*v^2)/2
%   => ρ = (2*Pdyn)/(v^2)
% On obtient ln(ρ) = ln(ρ0) - h/hs
% où α = ρ0, β = -1/hs et x = h
% où X = h, Y = ln(ρ), m = β = -1/hs, b = ln(ρ0) 
% 

if showTerminalOutput == 1
    disp('==========Équation des Y et X transformés pour l’approximation linéaire à deux paramètres==========')
    disp('X = h')
    disp('Y = ln(ρ)')
    disp(' ')
    disp('==========Coefficients de l’approximation à 2 paramètres==========')
    disp('m = β = -1/hs')
    disp('b = ln(α) = ln(ρ0)')
    disp(' ')
end

% Lissage en "Boîte grise"
y_bruite = (2*P_dyn)./(vit_mes(1:2:end).^2);
Yn = log(y_bruite);
Xn = alt_mes;
N_mat = length(Yn);
sum_Yn = sum(Yn);
sum_YnXn = sum(Xn.*Yn);
sum_Xn = sum(Xn);
sum_Xn_2 = sum(Xn.^2);
Yn_moy = (1/N_mat).*sum_Yn;

Mat_sums = [N_mat sum_Xn;
            sum_Xn sum_Xn_2];
Y_mat = [sum_Yn;
        sum_YnXn];
Params = pinv(Mat_sums)*Y_mat;

b = Params(1);
m = Params(2);
hs = -1/m;
roh0 = exp(b);
roh_approx = roh0.*exp(-alt_mes./hs);
P_approx = (1/2)*(roh_approx).*(vit_mes(1:2:end).^2);
D_aero_approx = P_approx.*S_aero_capsule.*C_do;
acc_approx = -D_aero_approx./m_capsule;

Yn_chapeau = m.*Xn + b;
R2 = (sum((Yn_chapeau - Yn_moy).^2))./(sum((Yn - Yn_moy).^2));
RMSY = sqrt((1/N_mat).*sum((Yn_chapeau - Yn).^2));
Err_abs = mean((-acc_mes(1:2:end)-acc_approx).^2);
Err_rel = mean(((-acc_mes(1:2:end)-acc_approx)./(-acc_mes(1:2:end))).^2);
RMS_acc_abs = sqrt(Err_abs);
RMS_acc_rel = sqrt(Err_rel);

if showTerminalOutput == 1
    disp(['m = ', num2str(m)])
    disp(['b = ', num2str(b)])
    disp(['R2 = ', num2str(R2)])
    disp(['RMS = ', num2str(RMSY)])
    disp(' ')
    disp('==========Densité de référence à la surface de la planète ρ0 et facteur d’échelle de la densité hs==========')
    disp(['ρ0 = ', num2str(roh0)])
    disp(['hs = ', num2str(hs)])
    disp(' ')
    disp('==========Erreur RMS absolue dans les accélérations (en m/s2)==========')
    disp([num2str(RMS_acc_abs), ' m/s^2'])
    disp(' ')
    disp('==========Erreur RMS relative dans les accélérations==========')
    disp(num2str(RMS_acc_rel))
    disp(' ')
end

if showGraphics == 1
    figure('Name','Accélération approximée par lissage et mesurée')
    hold on
    plot(t,-acc_mes,'bo')
    plot(t_simp,acc_approx,'r')
    xlabel('Temps (s)')
    ylabel('Accélération (Daero/masse) en m/s^2')
    legend('Accélération Mesurée','Accélération Approximée','Location','SouthWest')
    title('Accélération approximée par lissage et mesurée')
    grid on
    hold off
end

%% Loi de guidage : validation de la RAA
if showTerminalOutput == 1
    disp('==================================================================================================')
    disp('===============================Loi de guidage : validation de la RAA==============================')
    disp('==================================================================================================')
    disp(' ')
    disp('Génération des graphiques...')
    disp(' ')
end
% On utilise la RAA sans la gravité ici parce qu’un accéléromètre ne mesure pas la force de gravité
roh_ini = roh0*exp(-h_ini/hs);
vit_RAA = v_ini*exp((1/2)*B_NASA*hs*((roh_approx - roh_ini)/sin(gamma_ini_NASA)));

P_dyn_ini = (1/2).*(roh_approx).*vit_RAA.^2; 
D_aero_ini = P_dyn_ini.*S_aero_capsule.*C_do;
acc_RAA = -D_aero_ini./m_capsule; % À valider la démarche et la réponse du graphique

if showGraphics == 1
    figure('Name','V(h) calculée avec la RAA VS obtenue par intégration')
    hold on
    plot(alt_mes,vit_mes(1:2:end),'bo')
    plot(alt_mes,vit_RAA,'r')
    xlabel('hauteur m')
    ylabel('Vitesse v(h) m/s')
    legend('Vitesse Approximée (depuis NASA)', 'Vitesse avec la RAA','Location','SouthEast')
    title('V(h) calculée avec la RAA VS obtenue par intégration')
    grid on
    hold off
    
    figure('Name','Accélération calculée avec la RAA superposée sur les mesures accéléro de la NASA')
    hold on
    plot(alt_mes,-acc_mes(1:2:end),'bo')
    plot(alt_mes, acc_RAA,'r')
    xlabel('hauteur m')
    ylabel('Accélération m/s^2')
    legend('Accélération mesurée (NASA)', 'Accélération avec la RAA','Location','SouthEast')
    title('Accélération RAA superposée sur les mesures accéléro de la NASA')
    grid on
    hold off
end

%% Loi de guidage : limites structurelles
if showTerminalOutput == 1
    disp('==================================================================================================')
    disp('==============================Loi de guidage : limites structurelles==============================')
    disp('==================================================================================================')
end

% Newton-Raphson (Commun)
r_ini = R_Mars + h_ini;
r_fin = R_Mars + h_fin;
roh_fin = roh0*exp(-h_fin/hs);
roh_all = roh0*exp(-alt_mes./hs);
% Newton-Raphson (250m/s)
h_start_min_250 = 18000; % Point de départ choisi
[incr_min_250, Gamma_ref_250, h_min_250, v_min_250, ~, D_aero_min_250, ~, ~, ~,] = APP6_S5_Newton_Raphson(h_start_min_250,v_fin(1),roh0,roh_ini,roh_fin,hs,r_ini,r_fin);
h_start_max_250 = 53000; %alt_ini;
[incr_max_250, ~, h_max_250, v_max_250, ~, D_aero_max_250, ~, ~, ~,] = APP6_S5_Newton_Raphson(h_start_max_250,v_fin(1),roh0,roh_ini,roh_fin,hs,r_ini,r_fin);

v_moy_250 = (v_min_250 + v_max_250)/2;
delta_t_250 = (h_min_250-h_max_250)/(v_moy_250*sin(Gamma_ref_250));
vit_f_250 = v_ini*exp(1/2*B_NASA*hs*((roh_all-roh_ini)/sin(Gamma_ref_250)));
P_f_250 = 1/2*(roh_all).*vit_f_250.^2;

% Newton-Raphson (300m/s)
h_start_min_300 = 18000; % Point de départ choisi
[incr_min_300, Gamma_ref_300, h_min_300, v_min_300, ~, D_aero_min_300, ~, ~, ~,] = APP6_S5_Newton_Raphson(h_start_min_300,v_fin(2),roh0,roh_ini,roh_fin,hs,r_ini,r_fin);
h_start_max_300 = 53000; % Point de départ choisi
[incr_max_300, ~, h_max_300, v_max_300, ~, D_aero_max_300, ~, ~, ~,] = APP6_S5_Newton_Raphson(h_start_max_300,v_fin(2),roh0,roh_ini,roh_fin,hs,r_ini,r_fin);

v_moy_300 = (v_min_300 + v_max_300)/2;
delta_t_300 = (h_min_300-h_max_300)/(v_moy_300*sin(Gamma_ref_300));
vit_f_300 = v_ini*exp(1/2*B_NASA*hs*((roh_all-roh_ini)/sin(Gamma_ref_300)));
P_f_300 = 1/2*(roh_all).*vit_f_300.^2;

if showTerminalOutput == 1
    disp('=====VALEURS POUR V_FIN = 250m/s=====')
    disp(['γref = ', num2str(rad2deg(Gamma_ref_250)), ' degrés'])
    disp(['h_min = ', num2str(h_min_250), ' m'])
    disp(['v_min = ', num2str(v_min_250), ' m/s'])
    disp(['h_depart = ', num2str(h_start_min_250)])
    disp(['#Itérations = ', num2str(incr_min_250)])
    disp(['h_max = ', num2str(h_max_250), ' m'])
    disp(['v_max = ', num2str(v_max_250), ' m/s'])
    disp(['h_depart = ', num2str(h_start_max_250)])
    disp(['#Itérations = ', num2str(incr_max_250)])
    disp(['P_dyn_max : ', num2str(P_f_250(16)), ' N/m^2'])
    disp(['Δt^ lim = ', num2str(delta_t_250), ' s'])
    disp(' ')
    disp('=====VALEURS POUR V_FIN = 300m/s=====')
    disp(['γref = ', num2str(rad2deg(Gamma_ref_300)), ' degrés'])
    disp(['h_min = ', num2str(h_min_300), ' m'])
    disp(['v_min = ', num2str(v_min_300), ' m/s'])
    disp(['h_depart = ', num2str(h_start_min_300)])
    disp(['#Itérations = ', num2str(incr_min_300)])
    disp(['h_max = ', num2str(h_max_300), ' m'])
    disp(['v_max = ', num2str(v_max_300), ' m/s'])
    disp(['h_depart = ', num2str(h_start_max_300)])
    disp(['#Itérations = ', num2str(incr_max_300)])
    disp(['P_dyn_max : ', num2str(P_f_300(16)), ' N/m^2'])
    disp(['Δt^ lim = ', num2str(delta_t_300), ' s'])
    disp(' ')
end

%% Conception d’asservissements : loi de commande en translation
if showTerminalOutput == 1
    disp('==================================================================================================')
    disp('===================Conception d’asservissements : loi de commande en translation==================')
    disp('==================================================================================================')
end




%% Conception d’asservissements : loi de commande en rotation
if showTerminalOutput == 1
    disp('==================================================================================================')
    disp('====================Conception d’asservissements : loi de commande en rotation====================')
    disp('==================================================================================================')
end




%% Validation par simulation numérique
if showTerminalOutput == 1
    disp('==================================================================================================')
    disp('===============================Validation par simulation numérique================================')
    disp('==================================================================================================')
end








