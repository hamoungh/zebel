N=10;
N=N+1; 
D=0.2; 
Z=16; 
cvx_begin 
    variables  R(1,N) X(1,N) Q(1,N) 
    R(2:N) == D.*(1 + Q(1:N-1));    
% this equality constraint is not a?ne
   % X(1:N) == (1:N)./(Z + R(1:N));   
         X(1:N) .* (Z + R(1:N)) == (1:N);  
 %        X(1:N) Z + X(1:N) R(1:N) - (1:N) ==0;  
    Q(2:N) == X(2:N) .* R(2:N);   
    Q(1)==0; 
cvx_end

cvx_begin 
    variables  y R0 X0 Q0  R1 X1 Q1 
    Q0==0; 
    R0 == D.*(1 + Q0);
    y == [X0 R0];
%    X0 .* (Z + R0) == 0;
%       X1 .* (Z + R1) == 1;
%       ZX  + X R - 1 ==0;
%       ax1 + x1 x2 -1 = 0
  %     log(exp([1 0]*y+log(Z)) + exp([1 1]*y)) ==0   
%          
%          log(exp(
    Q1 == X1 .* R1;   
   
cvx_end








