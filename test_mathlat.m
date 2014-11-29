matrix = [1.5 1.764; 3.523 0.2];
digits(4);
s = sym(matrix,'d');
v = vpa(s,5); 
latex(v)