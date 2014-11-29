classdef test_qm_solve < handle

    properties
    end
    
    methods
        
        
        
        function X=encode_multiclass_X(o, r , x , q) 
            X = vertcat(reshape(r,[prod(size(r)),1]),... % prod(size(r)
                reshape(x,[prod(size(x)),1]),...
                reshape(q,[prod(size(q)),1]));
        end
        
        function [r , x , q] = decode_multiclass_X(o ,X,C, K,N)
               r = reshape(X(1: prod([C K N]) ), [C K N]);  
             x =  reshape(X(prod([C K N])+1:prod([C K N])+prod([C N])),[C N]); 
             q = reshape(X(prod([C K N])+prod([C N])+1 : prod([C K N])+prod([C N])+prod([C K N])), [C K N]);  
        end 
        
        function F = multiclass_closed_model(o,X,d,z,C,K,N)
            [r , x , q] =  o.decode_multiclass_X(X , C, K , N);
   %         q = [zeros(K,1) q];

%             for each class c we write C equations, it r(k , [N1 N2 ... N_C]) can be derived in C ways 
%             0 = r(k , [N1 N2]) - d(2,k)  .* (1+ Q(k , [N1 N2-1])) 
%             0 = r(k , [N1 N2]) - d(1,k)  .* (1+ Q(k , [N1-1 N2])) 
%             ... 

            nd = ndims(N); 
            oneCust=eye(nd); 
            tt=ones(size(N,2),1) * N; 
            tt(eye(size(tt))~=0)=1; 
            for c=1:nd 
                index = repmat({':'},1,nd); 
                indexPlusOne = index; 
                indexPlusOne{c} = 2:N(c);  
                 indexMinusOne = index; 
                indexMinusOne{c} = 1:N(c)-1;  
                %   eq_r =  r(k, 2:N1, 1:N2 ) - d(k,:)*ones(1,N) .* (1+q(k , 1:N1-1, 1:N2 ))  
                % d(k,c)  % [1 nd] 
              eq_r =  r(c, :, indexPlusOne{:}) -  repmat(d(c,:), [1 1 N-oneCust(c,:)])  .* (1+q(c, : , indexMinusOne{:} ))    ;        
              eq_x =  permute(x(c,index{:}),[2 3 1]) - ( repmat((1:N(c))',tt(c,:)) ./(z(c)+ permute(sum(r(c, : ,index{:}), 2) , [2+1:2+nd,1,2]) )) ; 
              % x(c,N1,N2) r(c,k,N1,N2)  
              % ex: > repmat(reshape(.5*ones(3,2),[3 1 2]), [1 3 1]) .* reshape(1:18,[3 3 2])
              eq_q = q(c,:, index{:})-  repmat(x(c, index{:}),[1 K]) .* r(:, index{:}); 
              F =  vertcat(reshape(eq1,[K*N,1]),...
                                   reshape(eq2,[N,1]),...
                                   reshape(eq3,[K*N,1]));
            end
        end

        function [x_ r_] = qm_multiclass_fsolve_(o,d,z,N)
            K=size(d,1); 
            C=size(d,2);  
            Nc=num2cell(N+1);  % to store the 0 related metrics  
            r=zeros(C,K,Nc{:});     % d * ones(1,N)
            x=zeros(C,Nc{:});     % ((1:N) ./ (z+sum(r(:,1:N),1)))
            q=zeros(C,K,Nc{:});         %      q(:,1)=[0 0]';
            
            X0 = o.encode_multiclass_X(r , x , q);
            [X,Fval,exitflag] = fsolve(@(X)o.multiclass_closed_model(X,d,z,C,K,N),X0);
            [r , x , q] = o.decode_multiclass_X(X,K,N+1);
            x_=x(N);
            r_=sum(r(:,N));
        end
        
       % --------------------------------------------------  
        function X=encode_X(o, r , x , q)
            K=size(q,1); 
            N=size(q,2); 
            X = vertcat(reshape(r,[K*N,1]),...
                reshape(x,[N,1]),...
                reshape(q,[K*N,1]));
        end
        
        function [r , x , q] = decode_X(o ,X,K,N)
               r = reshape(X(1:K*N),[K,N]);  
             x =  reshape(X(K*N+1:K*N+N),[1,N]); 
             q = reshape(X(K*N+N+1 : K*N+N+K*N),[K,N]);  
        end 
        
        function F = singleclass_closed_model(o,X,d,z,K,N)       
              [r , x , q] =  o.decode_X(X , K , N); 
               q = [zeros(K,1) q]; 
              eq1 =  r(:,1:N)-d*ones(1,N) .*(1+q(:,1:N));
              eq2 =  x(1:N)-((1:N) ./ (z+sum(r(:,1:N),1)));
              eq3 = q(:, 2:N+1)-ones(2,1)*x(1:N) .* r(:,1:N);
               F =  vertcat(reshape(eq1,[K*N,1]),...
                                    reshape(eq2,[N,1]),...
                                    reshape(eq3,[K*N,1]));  
        end
        
        function [x_ r_] = qm_fsolve_(o,d,z,N)
            K=size(d,1); 
            r=zeros(K,N);     % d * ones(1,N) 
            x=zeros(1,N);     % ((1:N) ./ (z+sum(r(:,1:N),1)))  
            q=zeros(K,N);         %      q(:,1)=[0 0]';
            
            X0 = o.encode_X(r , x , q);        
            [X,Fval,exitflag] = fsolve(@(X)o.singleclass_closed_model(X,d,z,K,N),X0); 
            [r , x , q] = o.decode_X(X,K,N); 
            x_=x(N);
            r_=sum(r(:,N));
        end       
       % --------------------------------------------------    
       
        function [x_ r_]= qm_solve_(o,d,z,N) 
            K=size(d,1);
            r=zeros(K,N);
            x=zeros(1,N);
            q=zeros(K,N+1);
            q(:,1)=[0 0]';
            nn=0:9;
            for n=1:N
                r(:,n)=d.*(1+q(:,n));
                x(n)=n/(z+sum(r(:,n),1));
                q(:,n+1)=x(n) * r(:,n);
            end
            
            x_=x(N);
            r_=sum(r(:,N));
        end
        
         % --------------------------------------------------    
        function test_qm_fsolve_(o)
            K=2;
            d=[1  .1  % c k 
                  1  1];   
            z=[40 60]; 
            N=[3 4];    
            o.qm_multiclass_fsolve_(d,z,N); 
            
             d=[1  
                 .1];   
             N=3;  
             z=40;  
            [x_ r_]= o.qm_fsolve_(d,z,N) 
            [x_ r_]= o.qm_solve_(d,z,N)         
          end
        
        % o=test_qm_solve()
        % [DFDY,FAC] = numjac( @(t,y)o.test_qm_solve1(t,y) ,[5],[1 1 40]',[0.2352 0.2963  0.2963 1.2597 1.2597]',1e-6*ones(3,1), .5*ones(3,1),0)
        function vec=test_qm_solve1(o,y, N)
            dd=y(1:2,1);
            z=y(3);
            [x r]=o.test_qm_solve_(dd,z, N);  
            vec=[x; r];
        end
        
        function main(o)  
%          Z = x1 r1
%          X = d1 d2 z1
            % generate 10 numbers with same N but changing z
            T=60; 
             d=[1
                   1];   
             N=10; 
                z=40+20*sin((1:T)/10);  
            X_orig = [d*ones(1,T) ; z];
            o=test_qm_solve(); 
            Z=[]; 
            for i=1:T  
                 [x_ r_]= o.test_qm_solve_(d,z(i),N); 
                Z=[Z [x_;r_]];    
            end
            X_orig
     %     Z
             
          z=40; 
          [x , r] = o.test_qm_solve_(d,z,N) 
          % [DFDY,FAC] = NUMJAC('F',T,Y,FTY,THRESH,FAC,VECTORIZED) 
          % numerically computes the Jacobian of function F(T,Y), 
          % here my function is a handle with (t,y) parameter signature which passes N hiddenly and discards the t 
           f_handle = @(t,XX)o.test_qm_solve1(XX,N);        
          [A,FAC] = numjac(f_handle  , [5], [d ; z] , [x ; r]  ,1e-6*ones(3,1), .5*ones(3,1),0);  
          A
          
           Pi_= .1 * eye(3);   
           Q= .5 * eye(3);   
           R(1,1)= .001;   R(2,2)=.001;
           % to denote relationship between thoughput  and response time   
           % R(1,2)=-0.04; R(2,1)=-0.04; 
           % denotes the a priori estimate of x0, that is, the estimate of x0 with knowledge of no measurements.
           X0Bar = [1 1 40]'; 
%           X(2+1,T)  
         
          cvx_begin  quiet 
                cvx_precision best
                variables    X(3,T+1) Y(2,T)  v(2,T)   w(3,T) 
                X>0; 
                Y - [x ; r] *ones(1,T) == A* (X(:,1:T) - [d ; z]*ones(1,T)) 
                 v == Y - Z ;
                 w == X(:,2:T+1)-X(:,1:T)
 
%                minimize sum(sum_square(v,2))  
               
            minimize ( (X(:,1)-X0Bar)' * Pi_ ^ -1 * (X(:,1)-X0Bar)+ ...
                              reshape( (w' * Q^-1)' , [1,3*T]) *   reshape( w, [1,3*T])' +...
                              reshape( (v' * R^-1)' , [1,2*T]) *   reshape(v, [1,2*T])' )

            cvx_end
            [X(:,1:T)' X_orig'] 
            [Y' Z'] 
            
   %         [y(1,:)' Z(1,1:T)']
        end
        
    end % methods
    
end % class
