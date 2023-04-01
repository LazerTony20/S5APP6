% Fichier créé par l'étudiant
% S5 APP6 - H2023
function f = apollo(t,z)

u = 1.0/(82.45);
u_des = 1 - u;

x = z(1);
y = z(2);
vx = z(3);
vy = z(4);

r1 = sqrt((x + u).^2 + y.^2);
r2 = sqrt((x - u_des).^2 + y.^2);

f(1) = vx;
f(2) = vy;
f(3) = 2*vy + x - u_des*(x + u)/(r1.^3) - u*(x - u_des)/(r2.^3);
f(4) = -2*vx + y - (u_des*y)/(r1.^3) - (u*y)/(r2.^3);

  f = f(:);