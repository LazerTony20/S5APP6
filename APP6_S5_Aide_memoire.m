% S5 APP6
% Aide-Mémoire
% 
% 
% 

%% 3.2 Utilisation des dérivées dans le calcul d’erreur d’intégration numérique 

% === Premiere derivée======

% Première dérivée au point a (différence avant) :
% f'(a)
% fpa = (y(2)-y(1))/ delx;

% Première dérivée au point b (différence arrière) :
% f'(b)
% fpb = (y(end)-y(end-1))/delx;

% === Deuxieme derivée======

% Deuxième dérivée au point a (différence avant)
% f''(a)
% fppa = (y(3) - 2*y(2)+y(1))/delx^2;

% Deuxième dérivée au point b (différence arrière) :
% f''(b)
% fppb = (y(end) - 2*y(end-1)+y(end-2))/ delx^2;

% === Troisieme derivée======
% Pour la 3e dérivée, on utilise la double dérivée appliquée à la simple dérivée.
%
% Troisième dérivée au point a (différence avant) :
% f'''(a)
% fpppa = ( y(4) − 3*y(3) + 3*y(2) − y(1) )/delx^3;

% Troisième dérivée au point b (différence arrière):
% f'''(b)
% fpppb = ( y(end) - 3*y(end-1) + 3*y(end-2) - y(end-3) )/delx^3;


%% Méthode de trapèze
%
% Intégration numérique
% F = (y(1) + y(Nm) + 2*sum(y(2:Nm-1)) )*h/2;
%
% Intégration numérique (pas par pas)
% F(1) = 0 % Initialisation de l’intégrale
% Nm = N + 1;
% for n = 2:Nm
%     F(n) = (y(1) + y(n) + 2*sum(y(2:n-1)) )*h/2;
% end

% Erreur 
% (EN CAS DE PAS PAR PAS, IL FAUT LA RECALCULER DANS LA BOUCLE)
% 
% Erreur = ((h^2)/12)*(f'(b) - f'(a))

%% Méthode de Simpson.
%
% Intégration numérique
% F = (y(1) + y(Nm) + 4*sum(y(2:2:Nm -1)) + 2*sum(y(3:2:Nm-1)))*h/3;    %
% DOIT AVOIR UN NOMBRE IMPAIR DE VALEURS
%
% Intégration continue
% Fs(1) = 0;        % Initialisation de l'intégrale
% xs(1) = x(1);
% for n = 3:2:Nm
%     Fs(((n-1)./2)+1) = (y(1) + y(n) + 4*sum(y(2:2:n-1)) + 2*sum(y(3:2:n-1)))*h/3;
%     xs(((n-1)./2)+1) = x(n);
% end

% Erreur
%
% Erreur = ((h^4)/180)*(fpppb - fpppa)
%
%



