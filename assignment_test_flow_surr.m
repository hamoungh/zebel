 classdef assignment_test_flow_surr < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        num_host
        num_container
        num_service
        num_classes  % 3 4 3
    % num_inst_types
       
        beta_deployed
        call_graph
        c 
        cost  
        d
        
        N
        Z
        RT_sla
        
        %f_sla=[10]';%[10 2 7]';
        % cap becomes the multiplicity of the host later on
        cap
        speed_ratio
      %  r1 
      
      % outputs after solving the model 
         placement 
         pl_delta
    end 
    
    methods
     % Take a as state 
     function f_sla = calculate_throuput_objective(o, N, Z, RT_sla, T) 
                 % I assume perfect knowlege of workload 
                 % in computing 'next step desired' throughput 
                 % f_sla .* (RT_sla + Z) ==  N ;                   
                  f_sla = [N./((RT_sla + Z ) * ones(1,T))];   
                  % f_c >= f_sla
                  
                  %  f_c(:,2:T+1) == f_sla;  % f_sla*ones(1,T);
     end

         function [u, alphaV,betaV,gammaV, f_c, mu_host ]=solve_nfm_over_time(o, f_sla, T, alphaV_0, cap) 
              d=o.d';  
               r1 = 1000*ones(1,o.num_classes); 
            cvx_begin  quiet 
                cvx_precision best
                variables  alphaV(o.num_host , o.num_container,T+1) ...   % control inputs  
                                betaV(o.num_container,o.num_service, T+1) ...
                                gammaV(o.num_service , o.num_classes,T+1) ... 
                                f_c(o.num_classes,T+1)...  %  output variable in MPC with one step delay based on input                      
                                consumption_cost(o.num_host,T+1) ... % cost related variables 
                                consumption_cost_tot(1,T+1)...
                                u(o.num_host , o.num_container,T) 

                  %-------------- constraints ----------------------
                  alphaV>=0
                  % sum(alphaV,2)<=o.cap
                 permute(sum(alphaV(:,:, 1:T+1),2),[1,3,2]) -  cap  <= 0 ;  %  o.cap*ones(1,T+1) <= 0 ;  
                  betaV>=0
                  
                  %--------------- system dynamics ---------------
                  alphaV(:,:,1) == alphaV_0 
                  alphaV(:,:,2:T+1) == alphaV(:,:,1:T) + u(:,:,1:T) 
                   
%                  %------------- objective calculation ---------------
%                  % I assume perfect knowlege of workload 
%                  % in computing 'next step desired' throughput 
%                  % f_sla .* (RT_sla + Z) ==  N ;                   
%                   f_sla == [N./((RT_sla + Z ) * ones(1,T))];   
%                   % f_c >= f_sla
%                   
%                   %  f_c(:,2:T+1) == f_sla;  % f_sla*ones(1,T);
                  
                  %------------------- system output -----------------
                  % sum(alphaV,1)' == sum(betaV ,2) 
                  permute(sum(alphaV,1),[2,3,1])  == permute(sum(betaV, 2),[1,3,2]) 
                  
                  % betaV .* o.beta_deployed == betaV;
                  betaV .* repmat(o.beta_deployed,[1,1,T+1]) == betaV;
             
                  % sum(betaV,1)' == sum(gammaV,2)
                  permute(sum(betaV,1),[2,3,1])  == permute(sum(gammaV, 2),[1,3,2]) 
                   
                 % here im talking about each class node  
                 % perfect information assumption assumes we know input or workload (N) ahead
                 % and proper output (f_sla) can be 
                 % describing how gamma is going to affect next step's f_c 
                 % pretty much like 'X(:,2:T+1) == A*X(:,1:T)+B*U' but we are stateless
                 % 
                 % gammaV == d .* (ones(o.num_service,1)*f_c')
                 gammaV(:,:,1:T+1) == ...
                     repmat(d,[1,1,(T+1)]) .* reshape(...
                            (ones(o.num_service,1)* reshape(f_c(:,1:T+1),[1,o.num_classes*(T+1)])), ...
                            [o.num_service,o.num_classes, (T+1)])
             %   f_c(:,1) == zeros(o.num_classes,1);
                 
                %--------------------- cost --------------------
                % minimize sum(sum(o.c.*alphaV)) +  o.r1 * pos( f_sla - f_net)
                % for linear case use pos instead of square_pos
                % not that cost is affine to the resources and flows (for now) so I
                % could extract it but the sla penalty is not (because of pos function) , so I have to put it in
                % objective   
                consumption_cost == permute(sum( repmat(o.c,[1 1 T+1]) .*alphaV,2),[1,3,2])
                consumption_cost_tot ==  sum(consumption_cost,1); 
                % minimize( sum(consumption_cost_tot ,2)) 
                
                %----------------- objective ------------------ make that
                %500 a 0 if you want no cost 
                minimize( sum(... % time
                                consumption_cost_tot(:,2:T+1) + ... % infrastructure cost
                                 r1 * pos( f_sla-f_c(:,2:T+1))   +  ... % sla violation cost
                                 0 * sum(sum(sum(abs(u(:,:,:))))) ... % change of control cost summed over host and container 
                                 ,2)) 
                % this is the compact notation used in the thesis: 
                %  minimize( sum(sum(permute(sum( repmat(c,[1 1 T]) .*alphaV,2),[1,3,2]),1) ,2))

                 % in fact should have been 
                 % minimize( X(:,2:T+1) - target_rt ) 
                 % minimize( RT_c(:,2:T+1) - RT_sla )  
             cvx_end 
             u= round((10^4).*u)/(10^4);
             alphaV = round((10^4).*alphaV)/(10^4);       
           %  consumption_cost 
            mu_host = permute(sum(alphaV(:,:, 1:T+1),2),[1,3,2]) ; 
        %     u
         %       f_sla
          %     f_c(:,2:T+1)
        %      RT_c = pos((N./f_c(:,2:T+1)) -  Z*ones(1,T))
%             sum_square(u)
%              sum(sum_square(u)) 
             % alphaV
             % pl_delta = (alphaV(:,:,2:T)>0)-(alphaV(:,:,1:T-1)>0) is not convex
             % thus we substitute it with square( alphaV(:,:,2:T) -  alphaV(:,:,1:T-1))
         end % solve_nfm_over_time
   
         function [U_scaling, VarthetaOut, Vartheta, cap, alphaV,betaV,gammaV, f_c, f_sla]=...
                 solve_nfm_over_time_pub(o, f_sla, T, alphaV_0, Vartheta_0) 
             provisioning_time=1;
             hold_duration=2;
                d=o.d';  
               r1 = 1000*ones(1,o.num_classes); 
              stations=provisioning_time+hold_duration;
             A=[zeros(1,stations-1)
                 eye(stations-1)];
             A=[A zeros(stations,1)];
             Ain=zeros(stations,1);
             Ain(1,1)=1;
             Aout = [zeros(1,provisioning_time) ones(1,hold_duration)];      
             num_inst_types = prices_amz().get_number_of_types;  
             
             capScale = 10;  
            cuPerInst = prices_amz().get_compute_units(); 
        %    cuPerInst=cuPerInst(2: o.num_inst_types+1); 
            pricePerInst = prices_amz().get_hourly_price();
         %   pricePerInst=pricePerInst(2: o.num_inst_types+1); 
             
           
%              Vartheta0=zeros(stations,1);
%              Vartheta0(provisioning_time+1,1) = 5;
             
             r=1;
       %      vartheta = varthetaPerInst * prices_amz().get_compute_units(); 
            cvx_begin  quiet
                cvx_precision best
                variables  alphaV(num_inst_types , o.num_container,T+1) ...   % control inputs  
                                betaV(o.num_container,o.num_service, T+1) ...
                                gammaV(o.num_service , o.num_classes,T+1) ... 
                                f_c(o.num_classes,T+1)...  %  output variable in MPC with one step delay based on input                       
                                consumption_cost_tot(1)...
                                u(num_inst_types , o.num_container,T)...
                                Vartheta(stations, num_inst_types, T+1)... 
                                VarthetaOut( num_inst_types, T+1)... 
                                U_scaling(num_inst_types,T)...
                                cap(num_inst_types ,T+1)

                  %-------------- constraints ----------------------
                  alphaV>=0
                  % U_scaling==ones(o.num_inst_types,T);  
                  U_scaling >=  0; 
                  % sum(alphaV,2)<=o.cap
                  VarthetaOut == reshape(Aout * reshape(Vartheta,[stations,num_inst_types*(T+1)])...   
                                                                        ,[num_inst_types,T+1]); 
                  cap ==  capScale * cuPerInst*ones(1,T+1)  .*  VarthetaOut 
                  permute(sum(alphaV,2),[1,3,2]) - cap <= 0 ;  % o.cap*ones(1,T+1)  % Aout * Vartheta  
                % permute(sum(alphaV,2),[1,3,2]) -   sum(o.cap*ones(1,T+1))  <= 0 ;  %
                  betaV>=0
                % Vartheta>0; 
                  %--------------- system dynamics ---------------
                   Vartheta(1:provisioning_time, : ,1) == zeros(provisioning_time ,  num_inst_types);       
                   Vartheta(provisioning_time+1,:,1) == Vartheta_0;  
                   Vartheta(provisioning_time+2:provisioning_time+hold_duration,:,1) ==...
                       zeros(hold_duration-1 ,  num_inst_types) ;  
               %    Vartheta(2,1,2) == 4;  
              
                   Vartheta(:,:,2:T+1) == reshape(A* reshape(Vartheta(:,:,1:T),[stations,num_inst_types*T])...
                                                        ,[stations,num_inst_types,T])+...
                                                     reshape(Ain*reshape(U_scaling(:,1:T),[1,num_inst_types*T]),...
                                                         [stations, num_inst_types, T]); 
         
                
         
                  alphaV(:,:,1) == alphaV_0 
                  alphaV(:,:,2:T+1) == alphaV(:,:,1:T) + u(:,:,1:T) 
                   
%                  %------------- objective calculation ---------------
%                  % I assume perfect knowlege of workload 
%                  % in computing 'next step desired' throughput 
%                  % f_sla .* (RT_sla + Z) ==  N ;                   
%                   f_sla == [N./((RT_sla + Z ) * ones(1,T))];   
%                   % f_c >= f_sla
%                   
%                   %  f_c(:,2:T+1) == f_sla;  % f_sla*ones(1,T);
                  
                  %------------------- system output -----------------
                  % sum(alphaV,1)' == sum(betaV ,2) 
                  permute(sum(alphaV,1),[2,3,1])  == permute(sum(betaV, 2),[1,3,2]) 
                  
                  % betaV .* o.beta_deployed == betaV;
                  betaV .* repmat(o.beta_deployed,[1,1,T+1]) == betaV;
             
                  % sum(betaV,1)' == sum(gammaV,2)
                  permute(sum(betaV,1),[2,3,1])  == permute(sum(gammaV, 2),[1,3,2]) 
                   
                 % here im talking about each class node  
                 % perfect information assumption assumes we know input or workload (N) ahead
                 % and proper output (f_sla) can be 
                 % describing how gamma is going to affect next step's f_c 
                 % pretty much like 'X(:,2:T+1) == A*X(:,1:T)+B*U' but we are stateless
                 % 
                 % gammaV == d .* (ones(o.num_service,1)*f_c')
                 gammaV(:,:,1:T+1) == ...
                     repmat(d,[1,1,(T+1)]) .* reshape(...
                            (ones(o.num_service,1)* reshape(f_c(:,1:T+1),[1,o.num_classes*(T+1)])), ...
                            [o.num_service,o.num_classes, (T+1)])
             %   f_c(:,1) == zeros(o.num_classes,1);
                 
                %--------------------- cost --------------------
                % minimize sum(sum(o.c.*alphaV)) +  o.r1 * pos( f_sla - f_net)
                % for linear case use pos instead of square_pos
                % not that cost is affine to the resources and flows (for now) so I
                % could extract it but the sla penalty is not (because of pos function) , so I have to put it in
                % objective   
      %          consumption_cost_tot ==  r*sum( pricePerInst' * U_scaling(:,1:T)); % permute(sum( repmat(o.c,[1 1 T+1]) .*alphaV,2),[1,3,2])
                % sum(consumption_cost,1); 
                % minimize( sum(consumption_cost_tot ,2)) 
                
                %----------------- objective ------------------
                minimize(  sum(... % time
                               r* pricePerInst' * U_scaling(:,1:T)+ ...  % infrastructure cost
                               r1 * pos( f_sla-f_c(:,2:T+1))   +  ... % sla violation cost
                                 permute( 0 * sum(sum(abs(u(:,:,:)))) , [1 3 2]) ... % change of control cost summed over host and container 
                                 ,2)) 
                % this is the compact notation used in the thesis: 
                %  minimize( sum(sum(permute(sum( repmat(c,[1 1 T]) .*alphaV,2),[1,3,2]),1) ,2))

                 % in fact should have been 
                 % minimize( X(:,2:T+1) - target_rt ) 
                 % minimize( RT_c(:,2:T+1) - RT_sla )  
             cvx_end 
             u= round((10^4).*u)/(10^4);
             alphaV = round((10^4).*alphaV)/(10^4);       
        %     u    
       % alphaV  
       % U_scaling 
        %      f_sla
%                f_c
%                Vartheta(:,:)
%                Vartheta
        %      RT_c = pos((N./f_c(:,2:T+1)) -  Z*ones(1,T))
%             sum_square(u)
%              sum(sum_square(u)) 
                % alphaV
             % pl_delta = (alphaV(:,:,2:T)>0)-(alphaV(:,:,1:T-1)>0) is not convex
             % thus we substitute it with square( alphaV(:,:,2:T) -  alphaV(:,:,1:T-1))
             
   
             
%              subplot(3,1,1); stairs(Uopt); title('Uopt');
%              subplot(3,1,2); plot(Xopt); title('Q opt');
%              subplot(3,1,3); plot(lambda_smooth); title('lambda');
%              figure; plot(([zeros(1,provisioning_time) ones(1,hold_duration)]*Vartheta)');
%              subplot(4,1,1); plot(Xopt); title('Xopt')
%              subplot(4,1,2); plot(Uopt'); title('Uopt')
%              subplot(4,1,3); plot(Xall); title('Xall')
%              subplot(4,1,4); plot(Uall'); title('Uall')
%              
%              subplot(3,1,3); plot(lambda(2:end))
         end
         
          function [VarthetaOut, Vartheta, cap, alphaV,betaV,gammaV, f_c, f_sla]=solve_nfm_over_time_known_u(o, U_scaling, N, Z, RT_sla, T, alphaV_0 ,Vartheta_0) 
             provisioning_time=1;
             hold_duration=2;
                d=o.d';  
               r1 = 1000*ones(1,o.num_classes); 
              stations=provisioning_time+hold_duration;
             A=[zeros(1,stations-1)
                 eye(stations-1)];
             A=[A zeros(stations,1)];
             Ain=zeros(stations,1);
             Ain(1,1)=1;
             Aout = [zeros(1,provisioning_time) ones(1,hold_duration)];      
             num_inst_types = prices_amz().get_number_of_types;  
             capScale = 10;  
            cuPerInst = prices_amz().get_compute_units(); 
            pricePerInst = prices_amz().get_hourly_price();

             
             r=1;
       %      vartheta = varthetaPerInst * prices_amz().get_compute_units(); 
            cvx_begin  quiet
                cvx_precision best
                variables  alphaV(num_inst_types , o.num_container,T+1) ...   % control inputs  
                                betaV(o.num_container,o.num_service, T+1) ...
                                gammaV(o.num_service , o.num_classes,T+1) ... 
                                f_c(o.num_classes,T+1)...  %  output variable in MPC with one step delay based on input    
                                f_sla(o.num_classes,T) ...                          
                                u(num_inst_types , o.num_container,T)...
                                Vartheta(stations, num_inst_types, T+1)...                                 
                                cap(num_inst_types ,T+1)

                  %-------------- constraints ----------------------
                  alphaV>=0
                  VarthetaOut = reshape(Aout * reshape(Vartheta,[stations,num_inst_types*(T+1)])...   
                                                                        ,[num_inst_types,T+1]); 
                  cap ==  capScale * cuPerInst*ones(1,T+1)  .*  VarthetaOut 
                  permute(sum(alphaV,2),[1,3,2]) - cap <= 0 ;  % o.cap*ones(1,T+1)  % Aout * Vartheta  
                  betaV>=0
                  %--------------- system dynamics ---------------
                   Vartheta(1:provisioning_time, : ,1) == zeros(provisioning_time ,  num_inst_types);       
                   Vartheta(provisioning_time+1,:,1) == Vartheta_0;  
                   Vartheta(provisioning_time+2:provisioning_time+hold_duration,:,1) ==...
                       zeros(hold_duration-1 ,  num_inst_types) ;   
                   
                   Vartheta(:,:,2:T+1) == reshape(A* reshape(Vartheta(:,:,1:T),[stations,num_inst_types*T])...
                                                        ,[stations,num_inst_types,T])+...
                                                     reshape(Ain*reshape(U_scaling(:,1:T),[1,num_inst_types*T]),...
                                                         [stations, num_inst_types, T]);  
         
                  alphaV(:,:,1) == alphaV_0 
                  alphaV(:,:,2:T+1) == alphaV(:,:,1:T) + u(:,:,1:T) 
                   
                 %------------- objective calculation ---------------                
                  f_sla == [N./((RT_sla + Z ) * ones(1,T))];   
                  
                  %------------------- system output -----------------
                  permute(sum(alphaV,1),[2,3,1])  == permute(sum(betaV, 2),[1,3,2]) 
                  betaV .* repmat(o.beta_deployed,[1,1,T+1]) == betaV;
                  permute(sum(betaV,1),[2,3,1])  == permute(sum(gammaV, 2),[1,3,2]) 
                  gammaV(:,:,1:T+1) == ...
                     repmat(d,[1,1,(T+1)]) .* reshape(...
                            (ones(o.num_service,1)* reshape(f_c(:,1:T+1),[1,o.num_classes*(T+1)])), ...
                            [o.num_service,o.num_classes, (T+1)])
                 
                %----------------- objective ------------------
                minimize( sum(... % time                                
                                 r1 * pos( f_sla-f_c(:,2:T+1))   +  ... % sla violation cost
                                 0 * sum(sum(sum(abs(u(:,:,:))))) ... % change of control cost summed over host and container 
                                 ,2)) 

             cvx_end 
             u= round((10^4).*u)/(10^4);
             alphaV = round((10^4).*alphaV)/(10^4);       

          end
         
          function VarthetaOut=calculate_VarthetaOut_by_U_scale(o, U_scaling, T,Vartheta_0)
                  provisioning_time=1;
             hold_duration=2;
              stations=provisioning_time+hold_duration;
             A=[zeros(1,stations-1)
                 eye(stations-1)];
             A=[A zeros(stations,1)];
             Ain=zeros(stations,1);
             Ain(1,1)=1;
             Aout = [zeros(1,provisioning_time) ones(1,hold_duration)];      
             capScale = 10;  
            cuPerInst = prices_amz().get_compute_units(); 
          num_inst_types = prices_amz().get_number_of_types;  

             
             r=1;
       %      vartheta = varthetaPerInst * prices_amz().get_compute_units(); 
            cvx_begin  quiet
                cvx_precision best
                variables  Vartheta(stations, num_inst_types, T+1)...
                                VarthetaOut(num_inst_types,T+1)

                  %-------------- constraints ----------------------      
                  VarthetaOut == reshape(Aout * reshape(Vartheta,[stations,num_inst_types*(T+1)])...   
                                                                        ,[num_inst_types,T+1]); 
                                                                    
                  %--------------- system dynamics ---------------
                   Vartheta(1:provisioning_time, : ,1) == zeros(provisioning_time ,  num_inst_types);       
                   Vartheta(provisioning_time+1,:,1) == Vartheta_0;  
                   Vartheta(provisioning_time+2:provisioning_time+hold_duration,:,1) ==...
                       zeros(hold_duration-1 ,  num_inst_types) ;   
                   
                   Vartheta(:,:,2:T+1) == reshape(A* reshape(Vartheta(:,:,1:T),[stations,num_inst_types*T])...
                                                        ,[stations,num_inst_types,T])+...
                                                     reshape(Ain*reshape(U_scaling(:,1:T),[1,num_inst_types*T]),...
                                                         [stations, num_inst_types, T]);  
             cvx_end     
          end
          
        function [alphaV,betaV,gammaV]=solve_nfm(o, N, Z, RT_sla) 
            d=o.d';
             r1 = 1000*ones(1,o.num_classes); 
            cvx_begin quiet
            % f_serv is the requested throughput of the user class
            % d is flow ratio parameters which convert the class flow f_c, in
            % units of user requests/sec, to demand flows ?sc for services
            
            % gamma is service rate in terms of cpu cycles per sec
            %  outputed from services to user classes
            variables  alphaV(o.num_host,o.num_container) ...
                betaV(o.num_container,o.num_service)...
                gammaV(o.num_service,o.num_classes) ...
                f_c(o.num_classes,1) ...
                f_sla(o.num_classes,1) ...
                f_net(o.num_classes,1)...
                RT_c(o.num_classes,1)
            
            alphaV>=0
            sum(alphaV,2)<=o.cap .* o.speed_ratio
            betaV>=0
            
            % this comes from equation1 of the network
            f_sla .* (RT_sla + Z) ==  N ;  % here each element of RT is \sum_k RT_c,k
         %   f_c >= f_sla
             f_net  ==  f_c;
            
            % assumotion is that all hardware level service rates is going to
            % become service throughputs
            sum(alphaV,1)' == sum(betaV ,2) %--
            betaV .* o.beta_deployed == betaV;
            % betaV(o.beta_deployed==0) == 0; 
            sum(betaV,1)' == sum(gammaV,2)
            
            % if service has a flow from container only if its deployed on it
           
            % For each single user request by class c, a demand of d_sc CPU-sec
            % is required for service s, giving this flow proportionality: _sc = d_sc f_c
            % this formula assumes that for every request-response from a client,
            % consumes service rate of the servers proportionally
            % or the service rate at the host level is distributed among the user
            % classes
             gammaV == d .* (ones(o.num_service,1)*f_c') %--
           
            
            % ones(1,num_host) sum(sum(c.*x)) ones(num_service,1)
            % for linear case use pos instead of square_pos 
            minimize sum( sum(o.c .* alphaV,2)) +  r1 * pos( f_sla - f_net) 
            %minimize sum(sum( alpha)) +  r1 * pos( f_sla - f_c) 
            cvx_end
            % RT_c=N./f_net - Z;
        end % solve_nfm
        
        function [f_,rt_]=solve_lqm(o,d,N,Z)
            try
                model=OpModel();
                
                % 2 PMs different capacities
                % 1 container on each
                % 1 class of user
                % variable workload
                
                model.nodes= [OpNode('ClientH','client',1,1)];
                for i=1:o.num_host
                    model.nodes=[model.nodes,  OpNode(sprintf('H%d',i),'server', o.cap(i) ,1)];
                end
                
                model.services=[OpService('ClientS','ClientT')];
                for i=1:o.num_service
                    model.services=[model.services,  OpService(sprintf('S%d',i),...
                        sprintf('T%d',i))];
                end
                
                % for now i put each service on a container
                model.containers=[OpContainer('ClientT', 1000, 'ClientH', 'false' )];
                for i=1:o.num_service
                    model.containers=[model.containers,   OpContainer(sprintf('T%d',i), 1000, sprintf('H%d',i), 'true' )];
                end
                model.clusters=[OpCluster('ClientCluster',[model.containers(1)]),...
                    OpCluster('TaskCluster',[model.containers(2:end)])];
                
                % this is going to be more complex for different examples,
                % later on, maybe comes from a matrix
                % but now its just equal calls from the client to all services
                for i=1:o.num_service
                    for j=1:num_classes
                        % OpCall(caller, callee, invocations, CPUDemand, DiskDemand)
                        model.calls=[model.calls...
                            OpCall('ClientS', sprintf('S%d',i), 1, d(o.num_service,o.num_classes), 0)];
                    end
                end
                model.scenarios=OpScenario('select','ClientS',model.calls);
                
                scWorkloads=containers.Map({'select'},{OpClosedWorkload(N,Z)});
                model.workload = OpWorkload(N,scWorkloads);
                
                nodes = ['ClientH '];
                for i=1:o.num_host
                    nodes=[nodes sprintf('H%d ',i)];
                end
                model.networks = OpNetwork(nodes);
            catch err
                disp(err);
            end
       
            
            %+++++++++++++++++
            model.solve();
            %         util = [];
            %         for i=3:length(opM.nodes)
            %             util=[util opM.nodes(i).cpuUtilization];
            %         end
            %         util_ = mean(util);
            f_=model.scenarios(1).throughput
            rt_=model.scenarios(1).responseTime
            
        end % solve_lqm
       
        
        function  Y=calculate_contact_level(o,call) % write using optimization  
              % total direct and indirect requests for a service
                % described in 'requests per user response'
                % Yuser1DB2Serv = 1 x 0.3 x (0.7 + 0.3 x 1.4) = 0.336
                % found by following all paths from the user class to the service.
                % calculating the contact level between class and service
                % nodes
                con=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
                for a=1:o.num_service
                    con_old=con;
                    con= con+ call^a;
                    if (con_old==con) 
                        break; 
                    end;
                end  
                Y=con(1:o.num_classes, o.num_classes+1:o.num_classes+o.num_service); 
        end
        
        function setup_hosts(o) 
              o.num_host=6;
              o.cap=[16 6 6 16 6 6]'; %[6 9 12]';
              o.speed_ratio=[1 1.2 0.9 1.1 0.8 1.2]'; 
               RandStream.setDefaultStream ...
                (RandStream('mt19937ar','seed',sum(10)));
                % rng(sum(100*clock),'v4')
                o.c=rand(o.num_host,o.num_container)+10;
                o.c;
                o.cost=round(rand(o.num_host,1)*10)+1;
                % c=round(rand(num_host,1)*10);
        end
            
        function setup_services(o)   
                o.num_service=14;
                o.num_container=10; 
                o.num_classes = 2;  % 3 4 3
                
                beta=zeros(o.num_container,o.num_service);
                beta(1,1)=1; 
                beta(2,2)=1; beta(2,3)=1; 
                beta(3,4)=1; 
                beta(4,5)=1; beta(4,6)=1; 
                beta(5,7)=1; 
                beta(6,8)=1; 
                beta(7,9)=1; 
                beta(8,10)=1; 
                beta(9,11)=1; beta(9,12)=1; beta(9,13)=1; 
                beta(10,14)=1; 
                o.beta_deployed=beta; 
                
                call=zeros(o.num_service+o.num_classes,o.num_service+o.num_classes);
                call(1,3)=1; call(1,4)=2;  call(1,5)=1;
                call(2,5)=2; call(2,6)=1; 
                call(3,7)=0.8; 
                call(4,7)=1; 
                call(5,8)=2.5; call(5,9)=1;
                call(6,8)=2; call(6,10)=1; call(6,15)=0.5; 
                call(7,11)=1; 
                call(8,11)=1; call(8,12)=2; 
                call(9,13)=0.5; call(9,14)=1; 
                call(11,16)=1.5;     
                o.call_graph = call; 
                Y=calculate_contact_level(o,call);
                
                
                D=...
                [0.005, 0.004, 0.002, 0.008, ...
                 0.001, 0.002, 0.005, 0.008, ...
                 0.006, 0.008, 0.004, 0.005, 0.006,...
                 0.005]; 
                 
                o.d= Y .* (ones(o.num_classes,1) * D); 
        
           
        end % setup_services
        
          function test_nfm_static(o)      
                o.setup_hosts(); 
               o.setup_services();
                N=[250 100]';
                Z=[1 1]';
                RT_sla=[0.146  0.267]';      
                [alphaV,betaV,gammaV]= solve_nfm(o, N, Z, RT_sla)
               % solve_lqm(o,d,N,Z,)
          end
        
         function test_nfm_mpc(o)   
             o.setup_hosts(); 
             nsteps = 1;
             T=1; 
             o.setup_services();
              N=[];
%              N(1,:) =  [sin((1:nsteps)/300)*200+300];
%              N(2,:) = [sin(((1:nsteps)+400)/300)*100+400];    %    plot(N')
              N=[250 100]' * ones(1,T); 
             Z=[1 1]';
             RT_sla=[0.146  0.267]'; 
            %   o.r1 = 1000*ones(1,o.num_classes); 
               [alphaV,betaV,gammaV] = solve_nfm_over_time(o, N, Z, RT_sla, T) 
         end
         
         % o= assignment_test_flow_surr()
         % o.test_nfm_mpc_static_compare
            function   test_nfm_mpc_static_compare(o)      
              nsteps = 1;
              T=3; 
              o.setup_hosts(); 
              o.setup_services();
              N=[];
              N=[250 100 200
                  100   250 50]; 
              Z=[1 1]';
              RT_sla=[0.146  0.267]'; 
             alphaV_0 = zeros(o.num_host , o.num_container); 
             
              f_sla = calculate_throuput_objective(o, N, Z, RT_sla, T); 

              [u, alphaV,betaV,gammaV,~,~] = solve_nfm_over_time(o,f_sla, T, alphaV_0, o.cap*ones(1,T+1));
            
             
            % alphaV
            placement = (alphaV(:,:,:)~=0);
            o.placement = placement;
            % placement
             pl_delta = placement(:,:,2:T+1) -  placement(:,:,1:T)
             o.pl_delta = pl_delta; 
             fprintf('initial_deployment=%d \t additions=%d \t removes=%d \t changes=%d \n',...
                  sum(sum(sum(pl_delta(:,:,1)))),...
                  sum(sum(sum(pl_delta>0))),...
                 sum(sum(sum(pl_delta<0))),...
                 sum(sum(sum(pl_delta~=0)))); 
      
            % additions 
             [i j k]=ind2sub(size(pl_delta),find(pl_delta(:,:,2:end)>0));
             additions=[i j k+1]
             % removes 
             [i j k]=ind2sub(size(pl_delta),find(pl_delta(:,:,2:end)<0)); 
             removals=[i j k+1]
             o.print_delta_sequece(pl_delta)
             o.print_mat( pl_delta(:,:,1) );    
            end
            
            function  t=get_container_of(o,s) 
                  [t,~] = find(o.beta_deployed(:,s)>0); 
            end
            
             function draw_task_deployment(o, varargin) 
                  filename = 'test_graph1.dot'; 
                 
                fid = fopen(filename, 'w');
                fprintf(fid, 'digraph G {\n');

                  for t=1:o.num_container  
                    [~,s]=find(o.beta_deployed(t,:)>0);
                     fprintf(fid, 'T%d [shape=record, label="{{',t); 
                     for s_i=1:size(s,2); 
                          fprintf(fid,'<s%d>s%d', s(s_i),s(s_i));  
                          if s_i ~=  size(s,2)
                              fprintf(fid,'|');
                          end 
                     end
                     fprintf(fid, '}}"];\n');
                     fprintf(fid, 'T%d -> T%d [label = "T%d" color = transparent]; \n',t,t,t );  
                  end
                  
               for i=1:o.num_host
                    fprintf(fid, 'H%d;\n', i  );     
                    for j=1:o.num_container 
                        if o.pl_delta(i,j,1)>0
                             fprintf(fid, 'T%d -> H%d; \n',j,i );  
                        end
                    end             
                end 
                  
                 fprintf(fid, '}\n');
             end 
             
            % o= assignment_test_flow_surr()
            % o.draw_service_callgraph() 
            % varargin{1} is set to true whenever you want show the
            % deployment on the hosts in a specific timestamp.  the
            % timestamps stored in varargin{2}
            function draw_service_callgraph(o, varargin) 
                %   'filename'  -  if omitted, writes to 'tmp.dot'
                o.setup_hosts(); 
               o.setup_services();
                N=[250 100 200
                  100   250 50]; 
               Z=[1 1]';
               RT_sla=[0.146  0.267]'; 
              
                i=1; 
                 filename = 'test_graph1.dot'; 
                 arctxt = '->';  
                 labeltxt = '[label="%s"]';
                 
                fid = fopen(filename, 'w');
                fprintf(fid, 'digraph G {\n');
           %     fprintf(fid, 'node [style=rounded];\n splines=false;');  
                fprintf(fid, ' splines=false;');
              
                
                for c=1:o.num_classes
                      % fprintf(fid, 'C%d [shape=record, label="{{<c1>Z=..., N=...}| C%d}"];\n',c,c); 
                       fprintf(fid, 'C%d [shape=record, label="{{<c1>Z=%d, N=%s}}"];\n',c ,Z(c),vect2str(N(c,:), 'formatstring', '%d') );   
                end
                  
                for t=1:o.num_container  
                    [~,s]=find(o.beta_deployed(t,:)>0);
                     fprintf(fid, 'T%d [shape=record, label="{{',t); 
                     for s_i=1:size(s,2); 
                          fprintf(fid,'<s%d>s%d', s(s_i),s(s_i));  
                          if s_i ~=  size(s,2)
                              fprintf(fid,'|');
                          end 
                     end
                     fprintf(fid, '}}"];\n');
                     fprintf(fid, 'T%d -> T%d [label = "T%d" color = transparent]; \n',t,t,t );  

                   %  fprintf(fid, '}| T%d}"];\n',t);
                end
                
                for i=1:o.num_classes
                    for j=1:o.num_service
                         inv = o.call_graph(i,j+o.num_classes); 
                         if  inv>0
                             if varargin{1}==true 
                                 fprintf(fid,'C%d:c1 ->T%d:s%d; \n',i, o.get_container_of(j), j); 
                             else 
                                 fprintf(fid,'C%d:c1 ->T%d:s%d[label="(%.1f)"]; \n',i, o.get_container_of(j), j, inv); 
                             end
                         end
                    end
                end
                
                 for i=1:o.num_service
                    for j=1:o.num_service
                         inv = o.call_graph(i+o.num_classes,j+o.num_classes); 
                         if  inv>0
                             if varargin{1}==true 
                                  fprintf(fid,'T%d:s%d ->T%d:s%d; \n ', o.get_container_of(i), i, o.get_container_of(j), j ); 
                             else
                                fprintf(fid,'T%d:s%d ->T%d:s%d [label="(%.1f)"]; \n ', o.get_container_of(i), i, o.get_container_of(j), j , inv); 
                             end
                         end
                    end
                 end
                 
                  timestamp=varargin{2};
                 if varargin{1}==true 
                     for i=1:o.num_host
                         fprintf(fid, 'H%d;\n', i  );
                         for j=1:o.num_container
                             % this is the case when the container is just
                             % added to host and should be drawn in red
                             % link
                             %
                             if (timestamp>1 && o.pl_delta(i,j,timestamp-1)>0)
                                 fprintf(fid, 'T%d -> H%d [style=dashed, dir=none,color=blue,penwidth=4]; \n',j ,i );
                             elseif  (timestamp>1 && o.pl_delta(i,j,timestamp-1)<0)
                                 fprintf(fid, 'T%d -> H%d [style=dashed, dir=none,color=red,penwidth=4]; \n',j ,i );
                             elseif  (o.placement(i,j, timestamp )>0)
                                 fprintf(fid, 'T%d -> H%d [style=dashed, dir=none]; \n',j,i );
                             else
                                 %do nothing
                             end
                         end
                     end
                 end %if 

%                      if timestamp>1
%                          [i j]=find(o.pl_delta(:,:,timestamp-1)~=0);
%                          c=o.pl_delta(sub2ind(size(o.pl_delta),  i, j, ones(size(i,1),1)*2 ));
%                          for k=1:size(i,1)
%                              if (o.pl_delta(i(k),j(k),timestamp-1)==1 )
%                                  % 'added to';
%                                  fprintf(fid, 'T%d -> H%d [style=dashed, dir=none,color=red,penwidth=4]; \n',j(k) ,i(k) );
%                              else
%                                  % 'removed from';
%                                  % fprintf(fid, 'T%d -> H%d [style=dashed, dir=none,color=blue,penwidth=4]; \n',j(k) ,i(k) );
%                              end
%                          end
%                      end
%                      
                     
               
  
                fprintf(fid, '}\n');
  
       
            end
            
            function test_nfm_over_time_pub(o)      
              nsteps = 1;
              T=3; 
              o.setup_services();
              N=[];
              N=[250 100 200
                  100   250 50]; 
              Z=[1 1]';
              RT_sla=[0.146  0.267]';  
%              o.num_inst_types=2; 
            num_inst_types = prices_amz().get_number_of_types; 
             % [u, alphaV,betaV,gammaV] = solve_nfm_over_time(o, N, Z, RT_sla, T, alphaV_0);
          
             Vartheta_0 = [5 0 0];
             cuPerInst = prices_amz().get_compute_units(); 
             f_sla_fixed = calculate_throuput_objective(o, N, Z, RT_sla, T); 
             f_sla_fixed
             
             f_sla = f_sla_fixed;   
             f_c_old=f_sla_fixed; 
             while true
            %      f_sla
                  alphaV_0 = zeros(num_inst_types , o.num_container);
                 [U_scaling, ~ ,  ~, ~, ~,~,~ , f_c] = solve_nfm_over_time_pub(o, f_sla, T, alphaV_0, Vartheta_0);
               %  f_c
                 %              % ceil
                 %              [VarthetaOut_, Vartheta_  , cap_, alphaV_,betaV_,gammaV_ , f_c, f_sla]=solve_nfm_over_time_known_u(o, floor(U_scaling), N, Z, RT_sla, T, alphaV_0, Vartheta_0);
                 %              f_c
                 %              f_sla
                 %              VarthetaOut_
                 VarthetaOut_=calculate_VarthetaOut_by_U_scale(o, floor(U_scaling),T,Vartheta_0);   
                 floor(U_scaling)  
                 num_machines=max(max(VarthetaOut_));
                 
                 cap1 = zeros(num_machines,num_inst_types, T+1 );   % T+1
                 for t=1:T+1
                     cap1(1,:,t)=VarthetaOut_(:,t);
                     for m=1:num_machines-1
                         cap1(m+1,:,t)=max(cap1(m,:,t)-1,0);
                     end
                     cap1(:,:,t) =  (cap1(:,:,t) > 0) * diag(cuPerInst);
                 end
                 cap1 = reshape(cap1 , [num_machines*num_inst_types T+1]);
              %    cap1
                 o.num_host =  num_machines*num_inst_types;
                 alphaV_0 = zeros(o.num_host , o.num_container);
                 capScale = 10;   
                 o.cap= cap1*capScale;
                 o.speed_ratio=[1 1.2 0.9 1.1 0.8 1.2]';
                 RandStream.setDefaultStream ...
                     (RandStream('mt19937ar','seed',sum(10)));
                 o.c=rand(o.num_host,o.num_container)+10;
                 o.c;
                 o.cost=round(rand(o.num_host,1)*10)+1;
                 [u, ~,~,~ , f_c, mu_host]=solve_nfm_over_time(o, f_sla, T, alphaV_0, o.cap);
                 %      alphaV
                   f_c
%                  mu_host
                if abs(f_c(:,2:T+1) - f_c_old)<abs(0.01 * f_c(:,2:T+1)) 
                    break;
                end
     
                 f_sla = f_sla + pos(f_sla_fixed - f_c(:,2:T+1));       
                 f_c_old = f_c(:,2:T+1);   
                 % c=round(rand(num_host,1)*10);
                 % cc = cumsum(VarthetaOut_(:,t))'
                 % uu=xor(vertcat(zeros(1,5),(cc.^-1)' *(1:5)<=1) , vertcat((cc.^-1)' *(1:5)<=1,zeros(1,5)))
             end % for loop
             f_sla 
             f_sla_fixed 
                 
            end
        
        function main1(o)
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

            d= rand(num_service,num_classes); 
        end  % main1

        function print_delta_sequece(o,pl_delta)
            s=size(pl_delta);
            fprintf('\\begin{eqnarray*}\n\\begin{split}\n'); 
            for t=2:s(3)
                 fprintf('\\text{t=%d:} \\\\\n',t);        
                [i j]=find(pl_delta(:,:,t)~=0); 
                c=pl_delta(sub2ind(size(pl_delta),  i, j, ones(size(i,1),1)*2 )); 
                for k=1:size(i,1) 
                    if (pl_delta(i(k),j(k),t)==1 )
                        text='added to';
                    else
                        text='removed from';
                    end     
                    fprintf('& %d   &  \\text{%s }  &  %d\\\\ \n',j(k) ,text,i(k))
                end
            end
            fprintf('\\end{split}\n\\end{eqnarray*}\n');  
            
%             \begin{eqnarray*} 
%             \begin{split}
%             
%             \text{t=2:} \\
%             & 8   &  \text{added to }  & 2\\
%             & 8  &   \text{added to } & 4 \\
%             
%             \text{t=3:}  \\
%              &  3   &  \text{remove from } & 1 \\
%             &   7  &    \text{added to } & 4 \\
%              &   8  &    \text{remove from } & 4 \\
%             \end{split}
%             \end{eqnarray*}  

        end
        
        function print_mat( obj, mat )           
            digits(4);
            s = sym(mat,'d');
            v = vpa(s,5); 
            latex(v)
        end
        
     
        
     
        
    end %methods
    end
