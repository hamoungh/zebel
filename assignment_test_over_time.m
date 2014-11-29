num_host=3;
num_service=4;
num_classes = 3; 
RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(10)));
 % rng(sum(100*clock),'v4')
 c=rand(num_host,num_service)*10; 
%c=round(rand(num_host,1)*10); 
%c=round(ones(num_host,1)*10); 

% c_r = reshape(c(:,:),[num_host*num_service,1]); 
w=workload(); 
nsteps = 1200;
T=1200; 
tmp= w.get_workload('data/day42_per_min.txt', 1, 23*60)';
 N=[];
 N(1,:) = [tmp(1:nsteps)];
 % lambda(nsteps+1:end)=0;
 N(2,:) =  [sin((1:nsteps)/300)*200+300];
 N(3,:) = [sin(((1:nsteps)+400)/300)*100+400];
 plot(N')
 RT_sla=[10 2 7]'; 
 Z=[15 15 15]'; 
 
 % I assume perfect knowlege of workload 
 % in computing 'next step desired' throughput 
  f_sla = [zeros(3,1) N./((RT_sla + Z ) * ones(1,T))];   
 
% f_sla=[10 2 7]'*ones(1,T) ;
cap=[60 90 120]';

r1 = [10000 10000 10000];

% taking container (or task) as application   

% total direct and indirect requests for a service 
% described in 'requests per user response'
% calculated for each user class, and used in every equality constraint
% involving that class
% example: 
% Users1 to DB2Serv it is:  
% Yuser1DB2Serv = 1 x 0.3 x (0.7 + 0.3 x 1.4) = 0.336 
% found by following all paths from the user class to the service.  
% in system with 3 user classes and 4 services/applications each class
% might use a combination of services, here we assume all users use all service equaly  
% d= ones(num_service,num_classes,T)/num_service;   
d= rand(num_service,num_classes); 


% d(1,1)=1; 

% x=round(rand(m,n))
 
%Y_cs be the total direct and indirect mean
 % requests to entry s for one request from user class c
Y=[]

cvx_begin  
    % f_serv is the requested throughput of the user class 
    % d is flow ratio parameters which convert the class flow f_c, in 
    % units of user requests/sec, to demand flows ?sc for services
     
     % gamma is service rate in terms of cpu cycles per sec
          %  outputed from services to user classes 
    variables  alpha(num_host , num_service,T) ...   % control inputs                     
                    gamma(num_service , num_classes,T) ... 
                    f_c(num_classes,T+1)...  %  output variable in MPC with one step delay based on input    
                    consumption(num_host,T) ... % cost related variables 
                    consumption_cost(num_host,T) ...
                    consumption_cost_tot(1,T) 
   
      alpha>=0
      f_c(:,1) == zeros(3,1);
     % assumotion is that all hardware level service rates is going to 
     % become service throughputs      
     % sum(alpha,1)' == sum(gamma,2) 
      permute(sum(alpha,1),[2,3,1])  == permute(sum(gamma, 2),[1,3,2])
    % reshape(sum(alpha,1),[2,3,1])  == reshape(sum(gamma, 2),[1,3,2])
    
     % For seach single user request by class c, a demand of d_sc CPU-sec
     % is required for service s, giving this flow proportionality: _sc = d_sc f_c
     % this formula assumes that for every request-response from a client,
     % consumes service rate of the servers proportionally 
    % or the service rate at the host level is distributed among the user
    % classes 
     % according to menase this is not right, if there is queuing 
     % service rate of a SERVER, demand of a FLOW on server 
     % anyways this d is something that has to be tracked 
     %  associated with 
     % for each service s, class c uses its harware flow proportional to its
     % demand 
     % gamma(:,:,1) == d .* (ones(4,1)*f_c(:,:,1)') % here im talking about each class node 
     % here im talking about each class node  
     % perfect information assumption assumes we know input or workload (N) ahead
     % and proper output (f_sla) can be 
     % describing how gamma is going to affect next step's f_c 
     % pretty much like 'X(:,2:T+1) == A*X(:,1:T)+B*U' but we are stateless
     % 

     gamma(:,:,1:T) == ...
         repmat(d,[1,1,T]) .* reshape(...
                (ones(num_service,1)* reshape(f_c(:,2:T+1),[1,num_classes*T])), ...
                [num_service,num_classes, T])
%        gamma==[...
%        0.0919   16.3750   14.9768
%        11.9246    6.6220    3.9246
%        18.9190   20.8222   10.3059
%         14.2605   16.2121   18.6092];
             
     % f_c >= f_sla
      f_c == f_sla;
     % minimize sum(sum(c.*alpha)) 
     
     consumption == permute(sum(alpha,2),[1,3,2])
     consumption_cost == permute(sum( repmat(c,[1 1 T]) .*alpha,2),[1,3,2])
     consumption -  cap*ones(1,T) <= 0 ;
      % for linear case use pos instead of square_pos
     
      
     %  consumption_host == c_r*ones(1,T) 
      
      % first term is over host and over time 
      % 
     % minimize(sum( sum(c*ones(1,T) .* consumption,1)+ r1 * pos( f_sla*ones(1,T)-f_c) ,2))
    % not that cost is affine to the resources and flows (for now) so I
    % could extract it but the sla penalty is not (because of pos function) , so I have to put it in
    % objective 
     consumption_cost_tot ==  sum(consumption_cost,1); 

     % minimize( sum(consumption_cost_tot + r1 * pos( f_sla-f_c) ,2))
     minimize( sum(consumption_cost_tot ,2))
     
     % in fact should have been 
     % minimize( X(:,2:T+1) - target_rt ) 
     % minimize( RT_c(:,2:T+1) - RT_sla )  
 cvx_end 
 % alpha  
 
% %  gamma
% %  f_c
% %  f_sla
% %  c
%   consumption_cost 
% % consumption
% consumption_cost_tot
% permute(sum(gamma, 2),[1,3,2])
%  
