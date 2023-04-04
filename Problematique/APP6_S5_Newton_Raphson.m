% S5 APP6
% Problematique - Newton-Raphson
% Anthony Royer et Jérémy Goulet
% ROYA2019 - GOUJ2711
% Usage
% [incr, Gamma_ref, h_sim, v_eval, P_eval, D_aero_eval, D_aero_eval, f_iterations, h_iterations] = APP6_S5_Newton_Raphson(h_start,v_fin,roh0,roh_ini,roh_fin,hs,r_ini,r_fin);
%
function [nr_out1, nr_out2, nr_out3, nr_out4, nr_out5, nr_out6, nr_out7, nr_out8, nr_out9] = APP6_S5_Newton_Raphson(h_start,vit_fin,roh0,roh_ini,roh_fin,hs,r_ini,r_fin)
% Données de la problématique
APP6_S5_Prob_Constantes

% Données du Newton-Raphson
incr = 0;
h_sim = h_start;
Delta_V_aero = vit_fin - sqrt(v_ini^2 + 2.*u_Mars.*(1./r_fin - 1./r_ini));
Gamma_ref = asin((1./2).*B_NASA.*hs.*((roh_fin - roh_ini)/log(1 + Delta_V_aero./v_ini)));

% F(x)
roh_eval = roh0*exp(-h_sim/hs);
v_eval = v_ini*exp(1/2*B_NASA*hs*((roh_eval-roh_ini)/sin(Gamma_ref)));
P_eval = 1/2*(roh_eval)*v_eval^2; 
D_aero_eval = P_eval*S_aero_capsule*C_do-D_aero_lim;
% D(x)
% dP_max = (roh0*vit_ini^2)/2.*exp((B.*hs.*roh-roh_ini)/sin(gamma_ref)).*exp(-alt_mes./hs).*((-B.*roh0.*exp(-alt_mes./hs))/sin(gamma_ref)-1/hs);
dD_aero_eval = S_aero_capsule*C_do*((roh0*v_ini^2)/2*exp((B_NASA*hs*roh_eval-roh_ini)/sin(Gamma_ref))*exp(-h_sim/hs)*((-B_NASA*roh_eval)/sin(Gamma_ref)-1/hs));

while and(abs(D_aero_eval) > tol_nr,incr < 500)
    % x = x - F/D;
    h_sim = h_sim-D_aero_eval/dD_aero_eval;
    roh_eval = roh0*exp(-h_sim/hs);
    % F = f(x);
    v_eval = v_ini*exp(1/2*B_NASA*hs*((roh_eval-roh_ini)/sin(Gamma_ref)));
    P_eval = 1/2*(roh_eval)*v_eval^2;
    D_aero_eval = P_eval*S_aero_capsule*C_do-D_aero_lim;
    % D = d(x);
    dD_aero_eval = S_aero_capsule*C_do*((roh0*v_ini^2)/2*exp((B_NASA*hs*roh_eval-roh_ini)/sin(Gamma_ref))*exp(-h_sim/hs)*((-B_NASA*roh_eval)/sin(Gamma_ref)-1/hs));
    incr = incr+1;
    f_iterations(incr) = D_aero_eval;
    h_iterations(incr) = h_sim;
end

nr_out1 = incr;       % Iterations
nr_out2 = Gamma_ref;  % Gamma_ref (EN RADIANS)
nr_out3 = h_sim;      % h_des (h_min ou h_fin selon l'appel de la fonction)
nr_out4 = v_eval;     % v_des (v_min ou v_max selon l'appel de la fonction)
nr_out5 = P_eval;     % P_dyn_des (min ou max selon l'appel de la fonction)
nr_out6 = D_aero_eval;% D_aero_eval (min ou max selon l'appel de la fonction)
nr_out7 = dD_aero_eval;% dD_aero_eval (min ou max selon l'appel de la fonction)
nr_out8 = f_iterations;
nr_out9 = h_iterations;
