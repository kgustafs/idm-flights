function [transfer_rates,system_state,time] = tauleapfractional(rates, end_time, input_counts, transfer_prob, transfer_int, tophubs, minpop, startfly_time, cycle_period, cycle_phase, cycle_x, alpha_0, beta_0, googleflu, flufuncarray, ffcarrayid)

% framework for a tau-leaping simulation
% from Bayati, JCP, 2013
% SSA adapted from Burrage and Carletti

global alpha beta syst_type

%transport = '1dfractional';
transport = '2dfractional';
%transport = 'airports';
% initial state is X0
num_nodes = size(input_counts,2);
% intended to determine the number of reaction channels
num_channels = size(input_counts,1);

total_pop_in = sum(input_counts(1,:),2);

rand('state', sum(100*clock));

current_time = 0;
alpha = 1; beta = 1;

Dalpha = 1.5;
center = 1;
alphafrac = 1.5;

muc = 1/(2*abs(cos(alphafrac*pi/2)));

tau = 1/(12*sqrt(2));
timestep = 1;

year = 0;
infecthub = tophubs(1); % change index to year+1 to move the infected hub each flu season

infecthub = 0;

if infecthub>0
    input_counts(:,infecthub) = [0.5*input_counts(1,infecthub); 0.5*input_counts(1,infecthub)];
end

if transport == '1dfractional'
    
    hfracgrid1d = 1;
    ffmu = (Dalpha/hfracgrid1d^alphafrac)*muc;
    
    vvector = zeros(num_channels,num_nodes);
    kvector = zeros(1,num_nodes);
    % fullstep refers to the moment after the diffusion takes place
    fullstep_counts = input_counts;
    % halfstep refers to the moment after tau-leaping reactions happen
    halfstep_counts = zeros(size(input_counts));
    
    transfer_rates(:,1,:) = zeros(size(input_counts));
    
    system_state(:,1,:) = input_counts;
    
elseif transport == '2dfractional'
    
    input_counts = zeros(2,num_nodes,num_nodes); %repmat(input_counts,[1,1,1,size(input_counts,3)]);
    input_counts(:,floor(num_nodes/2),floor(num_nodes/2)) = 1000;
    
    halfstep_counts = zeros(size(input_counts));
    fullstep_counts = input_counts;
    
    hfracgrid2da = 1;
    ffmu2da = (Dalpha/hfracgrid2da^alphafrac)*muc;
    
    hfracgrid2db = sqrt(2);
    ffmu2db = (Dalpha/hfracgrid2db^alphafrac)*muc;
    
    vvector = zeros(size(input_counts)); % this is 2d in space
    kvector = zeros(1,num_nodes);
    
    transfer_rates(:,1,:,:) = zeros(size(input_counts));
    
    system_state(:,1,:,:) = input_counts;
    
end

if rem(num_nodes,2) == 1
    lBigG = (num_nodes+1)/2
else
    lBigG = num_nodes/2;
end

all_a0 = zeros(1,num_nodes);
all_propensity = zeros(num_channels,num_nodes);

while current_time < end_time
    vvector = zeros(num_channels,num_nodes);
    if cycle_period>0
        alpha = alpha_0.*((1 - cycle_x) + cycle_x*sin(2*pi*current_time/cycle_period + cycle_phase));
        beta = beta_0;
        % loop over all nodes
        for idnode = 1:num_nodes
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts);
            % define basis for SSA interval
            all_a0(idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    elseif googleflu == 1
        % loop over all nodes
        for idnode = 1:num_nodes
            if current_time > 0
                node = ffcarrayid(idnode);
                timenow = ceil(current_time*52);
                ffnow = flufuncarray(timenow,node);
                maxffnow = max(flufuncarray(:,node));
                alpha = alpha_0.*ffnow./maxffnow; % normalize rate: why??
                beta = beta_0;
                %            beta  = 2.5*beta_0.*((1 - .3) + 0.3*sin(2*pi*[1:1:538]/52 + 104));
                %            beta = beta_0.*(1 + 2.*sin(2*pi*current_time/52/2+78).^12);
            else
                alpha = 0;
            end
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts);
            % define basis for SSA interval
            all_a0(idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    elseif strcmp(syst_type,'exponential')
        for idnode = 1:num_nodes
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts);
            % define basis for SSA interval
            all_a0(idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    elseif strcmp(syst_type,'SIS')
        for idnode = 1:num_nodes
            counts = fullstep_counts(:,idnode);
            % recompute propensity vector
            all_propensity(:,idnode) = propensity(counts);
            % define basis for SSA interval
            all_a0(idnode) = sum(all_propensity(:,idnode));
            %         if all_a0(idnode)==0
            %             disp('there is no probability of any reaction')
            %             idnode
            %             break
            %         end
        end
    else
        all_a0 = zeros(1,num_nodes);
        all_propensity = zeros(num_channels,num_nodes);
    end
    
    % the SSA timestep - without the random value
    % here we need to implement the selection procedure from Cao et al,
    % which will be problematic without a further condition to consider
    % the timestep from every node, and take the minimum
    % tau leaping
    % this is a temporary and strict solution for finding a timestep that is
    % acceptable for all cities' reaction systems
    tau = min(tau,1/max(all_a0));
    tau = min(1,tau);
    %    tau = max(tau,0.01);
    
    if transport =='1dfractional'
        
        for idnode = 1:num_nodes
            
            % this only computes infection for hubs, defined in find_probs
            counts = fullstep_counts(:,idnode);
            ktau = zeros(num_channels,1);
            % scale the number of reactions during the timestep by a Poisson
            % distributed random variate
            for idx = 1:num_channels
                % Line 6 in Bayati Sec:IV.B
                poismean = tau*all_propensity(idx,idnode);
                %ktau(idx,:) = random('poiss',poismean);
                ktau(idx,:) = poissrnd(poismean);
                % Line 8 in Bayati Sec:IV.B
                counts = counts + ktau(idx).*rates(:,idx);
                if counts(idx) < 0
                    % crude way of preventing negativity - Cao et al describe a
                    % better way
                    counts(idx) = 0;
                end
            end
            halfstep_counts(:,idnode) = counts;
        end
        
        vvector = zeros(num_channels,num_nodes);
        for idx = 1:num_channels
            BigG2 = frac_kernel(tau,ffmu,hfracgrid1d,num_nodes,alphafrac,center,Dalpha);
            for lambda=1:num_nodes
                BigG2rc = BigG2(num_nodes-lambda+1:end-lambda+1);
                BigG2rc = BigG2rc./sum(BigG2rc); % renormalization: fair?
                kvector = mnrnd(halfstep_counts(idx,lambda),BigG2rc);
                % mistake in Bayati paper?kvector(lambda) = kvector(lambda) - halfstep_counts(idx,1,lambda);
                vvector(idx,lambda) = vvector(idx,lambda) - sum(kvector(1:end~=lambda));
                vvector(idx,1:end~=lambda) = vvector(idx,1:end~=lambda) + kvector(1:end~=lambda);
            end
        end
        fullstep_counts = halfstep_counts + vvector;
        %     % attempt to renormalize
        %     fullstep_counts
        %     total_pop_in
        %     sum(fullstep_counts,2)
        %     fullstep_counts = 10.*fullstep_counts;
        % step forward in time
        current_time = current_time + tau;
        
        if mod(timestep,100)==0
            current_time
            tau
        end
        
        % the index for the system trajectory
        timestep = timestep + 1;
        % keep track of the time, which is only used for diagnostics in the
        % end, since the value of end_time in ssa_flights will determine the
        % vector of time points
        time(timestep) = current_time;
        % record the system state
        %     if current_time>year+1
        %         year = year+1
        %         infecthub = tophubs(1);
        %         fullstep_counts = input_counts;
        %         halfstep_counts = zeros(num_channels,num_nodes);
        %         fullstep_counts(1,infecthub) = 0.5*input_counts(1,1,infecthub);
        %         fullstep_counts(1,infecthub) = 0.5*input_counts(1,1,infecthub);
        %     end
        system_state(:,timestep,:) = fullstep_counts;
        transfer_rates(:,timestep,:) = vvector;
        
        %if transport=='airports'
        
        % % if transfer is largely negative, it will drive the population
        % % negative, but this must be prevented
        %     for idnode=1:num_nodes
        %         if fullstep_counts(1,idnode)<minpop;
        %             fullstep_counts(1,idnode)=minpop;
        %         end
        %         % however, the number of infected should be allowed to go to zero
        %         % if the epidemic dies away
        %         if fullstep_counts(2,idnode)<minpop;
        %             fullstep_counts(2,idnode)=0;
        %         end
        %     end
        %     min(fullstep_counts(1,:));
        
    elseif transport=='2dfractional'
        
        % INSERT REACTIONS HERE
        
        halfstep_counts = fullstep_counts; % placeholder skipping reactions
        
        vvector = zeros(num_channels,num_nodes);
        
        for nodepointx=1:num_nodes
            for idx = 1:num_channels
                        if path<3
                            shortpath = 1;
                            BigG2 = frac_kernel(tau,ffmu2da,hfracgrid2da,num_nodes,alphafrac,center,Dalpha);
                            if path<2
                                XXX = nodepointx;
                                distlB = nodepointx;
                                distrB = num_nodes/2-nodepointx;
                                lineput = squeeze(halfstep_counts(idx,XXX,:));
                                
            for nodepointy=1:num_nodes
                
                XXX = zeros(num_nodes,1);
                YYY = zeros(num_nodes,1);
                
                for path = 1:2
                    
                    % for each node point (x0,y0), compute dispersion along four lines
                    %  x = x0, all y
                    %  y = y0, all x
                    %  xi = yi
                    %  xi = -yi
                    
                            else
                                YYY = nodepointy;
                                distlB = nodepointy;
                                distrB = num_nodes - nodepointy;
                                lineput = squeeze(halfstep_counts(idx,:,YYY));
                            end
                        else
                            longpath = 1;
                            BigG2 = frac_kernel(tau,ffmu2da,hfracgrid2da,num_nodes,alphafrac,center,Dalpha);
                            if path>3
                                distlB = 0;
                                distrB = 0;
                            else
                                distlB = 0;
                                distrB = 0;
                            end
                        end
                        
                        % then compute the fractional kernel for the particular
                        % line of dispersion, which is a slice through
                        % halfstep_counts that is defined in a new line
                        
                        % this loop can probably remain the same...vvector will span
                        % the entire 2d grid, and should know about movements due to
                        % computing the dispersion at other nodes, which is probably as
                        % good as we can do for now - we must choose whether to do all
                        % the dispersions independently or sequentially. In any case,
                        % vvector will be updated in this lambda loop for a new index,
                        % which is specified in the new nodepoint loop (which probably
                        % needs to be a double loop.
                        
                        for lambda=1:num_nodes
                            BigG2rc = BigG2(num_nodes-lambda+1:end-lambda+1);
                            BigG2rc = BigG2rc./sum(BigG2rc); % renormalization: fair play?
                            kvector = mnrnd(lineput(lambda),BigG2rc)
                            %                             if jumpsout
                            %                                 keep in same node
                            %                             end
                            % add boundary condition: rejection if hopping off grid
                            % mistake in Bayati paper?kvector(lambda) = kvector(lambda) - lineput(lambda);
                            vvector(idx,lambda) = vvector(idx,lambda) - sum(kvector(1:end~=lambda))
                            vvector(idx,1:end~=lambda) = vvector(idx,1:end~=lambda) + kvector(1:end~=lambda)
                        end
                    end
                    if path<2
                        fullstep_counts(:,XXX,:) = halfstep_counts(:,XXX,:) + reshape(vvector,size(halfstep_counts(:,XXX,:)));
                    else
                        fullstep_counts(:,:,YYY) = halfstep_counts(:,:,YYY) + reshape(vvector,size(halfstep_counts(:,:,YYY)));
                    end
                end
            end
        end
    end
    
    
    
    %     % attempt to renormalize
    %     fullstep_counts
    %     total_pop_in
    %     sum(fullstep_counts,2)
    %     fullstep_counts = 10.*fullstep_counts;
    % step forward in time
    current_time = current_time + tau;
    
    if mod(timestep,100)==0
        current_time
        tau
    end
    
    % the index for the system trajectory
    timestep = timestep + 1
    % keep track of the time, which is only used for diagnostics in the
    % end, since the value of end_time in ssa_flights will determine the
    % vector of time points
    time(timestep) = current_time;
    % record the system state
    %     if current_time>year+1
    %         year = year+1
    %         infecthub = tophubs(1);
    %         fullstep_counts = input_counts;
    %         halfstep_counts = zeros(num_channels,num_nodes);
    %         fullstep_counts(1,infecthub) = 0.5*input_counts(1,1,infecthub);
    %         fullstep_counts(1,infecthub) = 0.5*input_counts(1,1,infecthub);
    %     end
    system_state(:,timestep,:,:) = fullstep_counts;
    %transfer_rates(:,timestep,:,:) = vvector;
    
    
end

end
%size(t)
%size(X)
