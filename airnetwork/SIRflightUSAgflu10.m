% simple reaction simulation
% from Bayati JCP 2013

% 1.) a = (a1,a2,...)^T : probability of reactions
% 2.) deltat_r = Tt : sample from exponential
% 3.) l = Qq ~ a_j(t)/a_0(t) : sample point-wise distribution of a
% 4.) n_bar(t + deltat_r) = n_bar(t) + nu_l

% multinomial stochastic algorithm with Lie-Trotter splitting and tau-leap

% 1.) a = (a1,a2,...)^T : the probablility of reactions
% 2.) deltat_r = Tt : sample from exponential or other type of reaction
% 3.) deltat_d = Dd : time-step restriction
% 4.) deltat = min(deltat_r,deltat_d)
% for all i do
% 5.) k_i ~ Pp(a_i*deltat) : number of executions for each channel i
% end for
% 6.) n_bar(t + deltat/2) = n_bar(t) + sum_l[k_l*nu_l] : state vector,
%                                                        number of particles
% 7.) nu_d = 0 : diffusion transitions
% for all lambda do : loop over cells, spatial discretization
% 8.) k1 ~ Bb(n_barlambda(t + deltat/2),Gg_lambdam1_lambda(deltat) :
%           binomial distribution Bb(L,p) of L trials, probability p
% 9.) k2 ~ Bb(n_barlambda(t + deltat/2) - k1,
%               (Gg_lambdap1_lambda(deltat))/(1-Gg_lambdam1_lambda(deltat))
% 10.) nu_lambda_d = nu_lambda_d - (k1 + k2) : movement out of lambda
% 11.) nu_lambdam1_d = nu_lambdam1_d + k1 : movement into neighbor
% 12.) nu_lambdap1_d = nu_lambdap1_d + k2 : movement into other neighbor
% end for
% 13.) n_bar(t + deltat) = n_bar(t + deltat/2) + nu_d : new state vector

clear all
clear reset

% index for species
global species1 species2 syst_type

% input parameters

syst_type = 'SIS';
algorithm = 'tauleap';
country = 'USA';

% use the open source fftf filtering program to get top 10 Fourier coeffs
load ../usagflu.mat
top10usagflu = fftf([1:1:538],usagflu,[],10);
% readjust the signal to account for floor
top10usagflufloor = top10usagflu - 0.3;
slip = find(top10usagflufloor<0);
top10usagflufloor(slip) = 0.;
top10usagflufloor = top10usagflufloor + 0.1;

species1 = 1; % susceptible
species2 = 2; % infected

R_0 = 6;  % reproductive number magnitude
beta_0 = 4; % recovery rate
alpha_0 = R_0*beta_0; % infection rate per unit infected

% allow alpha to cycle for seasonality
cycle_period = 0; % units of years
cycle_x = 0;
cycle_phase = 0;

% for using Google flu data or similar
gflu = 1; % on is 1

% only a single city ?
singlecity = 0;
% for a fake simulation ?
fakey=0;

% rate matrix (stoichiometry)
if strcmp(syst_type,'SIR')==1
    rates = [-1 0; 1 -1];
elseif strcmp(syst_type,'SIS')==1
    % a symmetric, sometimes oscillatory system
    rates = [-1 1; 1 -1];
end
    
[num_species,num_channels] = size(rates);

% probability matrix of transfer from city to city (cell to cell)
% how to determine the probability of remaining?
% v_lambda; v_lambda-1; v_lambda+1
% these test cases are set up so that the row indicates output travelers
% from that city to the nth city

% this factor is scaling the transfer matrix to account for the rescaling
% of the population of the airport service regions
% see tauleapdiffusion.m
trfactor = 0/10000000;

if fakey==1
    transfer_prob1 = [0 0; 2 2; 2 2; 2 2];
    transfer_prob2 = [0 0; 0 0; 0 0; 0 0];
    transfer_prob3 = [0 0; 0 0; 0 0; 0 0];
    transfer_prob4 = [4 0; 0 0; 0 0; 0 0];
    
    transfer_prob = cat(3,transfer_prob1',transfer_prob2',transfer_prob3',transfer_prob4').*trfactor;
else
    if strcmp(country,'USA') == 1
        load ../transferUSA.mat
        transfer_prob = symtransferUSA.*trfactor;
    elseif strcmp(country,'USA') == 1
        load ../transferUSA.mat
        transfer_prob = symtransferUSA.*trfactor;
    end
end


load ../usadata.mat usaorigin

hub_size = 1;
% infected airport city ids
infecthub = 0; %default
infecthub = tophubsUSA(1);
infecthub = 397;

% count_unique is an add-on function
[uniqueIATA,nn] = count_unique(usaorigin);
hubsUSA = find(nn>hub_size);

num_nodes = size(nn,1);

clear usaorigin

% initial state baseline
init_inf_percent = 1; % measure infected percentage-wise
initial_states = zeros(num_species,1,num_nodes);
% the total number of individuals using the hub is scaled according to the 
% connectivity of the hub
total_indiv = 1000; %round(1000.*nn./mean(nn(tophubsUSA)))+100; % 100 sets a minimum
init_inf = round(total_indiv*init_inf_percent/100); % set the number of infected
initial_states(species1,1,:) = total_indiv-init_inf; 
initial_states(species2,1,:) = init_inf; % assign number of infected
% a special infected hub can be selected
if infecthub>0
    initial_states(:,:,infecthub) = [0.5*initial_states(1,1,infecthub); 0.5*initial_states(1,1,infecthub)];
end

if singlecity==1
    transfer_prob = transfer_prob(:,:,infecthub);
end

% requiring a floor on population to avoid divide by zero
minpop = 1;
% time at which to start the transport
startfly_time = 0;
% total time in years
end_time = 5;
% interval between transfers due to air travel
transfer_int = 0.0192;

% the initial state of the SIR model for all cities
%initial_states = cat(3,X_initial,[100; 0],[100; 0],[100; 0]);
%initial_states = repmat(X_initial,[1,1,num_nodes]);
% infect half the population of the infected airport 
%initial_states(:,:,infect1) = [ceil(total_indiv/2); floor(total_indiv/2)];

% make the main hubs have 10x the population
% if singlecity==0
%     for kk = 1:size(tophubsUSA,1)
%         initial_states(:,:,tophubsUSA(kk)) = [total_indiv*10; init_inf*10];
%     end
% end

if strcmp(algorithm,'ssa')==1
    error('invalid algorithm')
elseif strcmp(algorithm,'tauleap')==1
    [transfer_rates,final_state,time_vector] = ...
        tauleapdiffusion10(rates, end_time, initial_states, transfer_prob, transfer_int, hubsUSA, minpop, startfly_time, cycle_period, cycle_phase, cycle_x, alpha_0, beta_0, gflu, top10usagflufloor);
else
    error('invalid algorithm')
end

