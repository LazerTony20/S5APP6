% Fichier créé par l'étudiant
function f = eqn1b(t,y)

  f(1) = 998*y(1) + 1998*y(2);
  f(2) = -999*y(1) - 1999*y(2);
  
  f = f(:);