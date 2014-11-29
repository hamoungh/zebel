num_host=2;
num_service=2;
num_classes = 1;  % 3 4 3 
RandStream.setDefaultStream ...
     (RandStream('mt19937ar','seed',sum(10)));
 % rng(sum(100*clock),'v4')
 c=round(rand(num_host,num_service)*10)
% c=round(rand(num_host,1)*10); 

N=100; 
Z=8; 
RT_sla=[2]'; 
% f_sla =  N./(RT_sla + Z ) ;   
  
%f_sla=[10]';%[10 2 7]';
% cap becomes the multiplicity of the host later on 
cap=[6 10]'; %[6 9 12]';        

r1 = 1000*ones(1,num_classes);

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
d= ones(num_service,num_classes); %/num_service;   
%d(1,1)=1; 

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
    variables  alpha(num_host,num_service) ...          
                        gamma(num_service,num_classes) ...
                        f_c(num_classes,1) ...
                        f_sla(num_classes,1) ...
                        f_net(num_classes,1)...
                        RT_c(num_classes,1)
                    
     alpha>=0
     
     % this comes from equation1 of the networl
     f_sla .* (RT_sla + Z )==  N ;  % here each element of RT is \sum_k RT_c,k
     
     % assumotion is that all hardware level service rates is going to 
     % become service throughputs  
     % sum(alpha,1)' == sum(gamma,2) %--
     % For each single user request by class c, a demand of d_sc CPU-sec
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
    %  gamma == d .* (ones(num_service,1)*f_c') %--
      sum(alpha,1)' == sum(d .* (ones(num_service,1)*f_c'),2)
      
     % f_c >= f_sla
     % minimize sum(sum(c.*alpha)) 
     
     sum(alpha,2)<=cap
     
     f_net + 1 ==  f_c;
   
      % ones(1,num_host) sum(sum(c.*x)) ones(num_service,1)
      % for linear case use pos instead of square_pos
      minimize sum(sum(c.*alpha)) +  r1 * pos( f_sla - f_net) 
     %minimize sum(sum( alpha)) +  r1 * pos( f_sla - f_c) 
    
 cvx_end 
 RT_c=N./f_net - Z;
%  alpha
%  gamma
%  d
 f_net
 RT_c
 %c
% c*ones(1,num_service)


try
%             opM=OpModel();
%             CloudSimulate().buildtest2(opM,instance_quantity_to_fire, demand, obj.workload(obj.t));            
            %+++++++++++++++++
         
            model=OpModel();            
        
        % 2 PMs different capacities
        % 1 container on each  
        % 1 class of user
        % variable workload                         

            model.nodes= [OpNode('ClientH','client',1,1)];
            for i=1:num_host
                model.nodes=[model.nodes,  OpNode(sprintf('H%d',i),'server', cap(i) ,1)]; 
            end

            model.services=[OpService('ClientS','ClientT')];    
            for i=1:num_service
                model.services=[model.services,  OpService(sprintf('S%d',i),...
                                                                        sprintf('T%d',i))]; 
            end

            % for now i put each service on a container 
            model.containers=[OpContainer('ClientT', 1000, 'ClientH', 'false' )];            
            for i=1:num_service
                model.containers=[model.containers,   OpContainer(sprintf('T%d',i), 1000, sprintf('H%d',i), 'true' )];  
            end           
            model.clusters=[OpCluster('ClientCluster',[model.containers(1)]),...
                            OpCluster('TaskCluster',[model.containers(2:end)])]; 
                        
            % this is going to be more complex for different examples,
            % later on, maybe comes from a matrix
            % but now its just equal calls from the client to all services 
            for i=1:num_service
                for j=1:num_classes 
                    % OpCall(caller, callee, invocations, CPUDemand, DiskDemand)
                    model.calls=[model.calls... 
                        OpCall('ClientS', sprintf('S%d',i), 1, d(num_service,num_classes), 0)]; 
                end
            end            
            model.scenarios=OpScenario('select','ClientS',model.calls); 
            
             scWorkloads=containers.Map({'select'},{OpClosedWorkload(N,Z)});            
             model.workload = OpWorkload(N,scWorkloads);         
                          
             nodes = ['ClientH ']; 
             for i=1:num_host
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
             model.scenarios(1).throughput
             model.scenarios(1).responseTime
             
       %     k=1
%catch err
%end  
        
