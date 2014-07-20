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

global species1 species2 syst_type

% input parameters

syst_type = 'noreax';
algorithm = 'tauleap';
country = 'USA';

% use the open source fftf filtering program to get top 10 Fourier coeffs
% for the google flu data for the entire country

species1 = 1; % susceptible
species2 = 2; % infected

R_0 = 5;  % reproductive number magnitude
beta_0 = 6; % recovery rate
alpha_0 = R_0*beta_0; % infection rate per unit infected

% allow alpha to cycle for seasonality
cycle_period = 0; % units of years
cycle_x = 0.5;
cycle_phase = 0.1;

% for using Google flu data
googleflu = 0; % on is 1

% only a single city ?
singlecity = 0;
% for a fake simulation ?
fakey=0;

rates = zeros(2,2);
% rate matrix (stoichiometry)
if strcmp(syst_type,'SIR')==1
    rates = [-1 0; 1 -1];
elseif strcmp(syst_type,'SIS')==1
    % a symmetric, sometimes oscillatory system
    rates = [-1 1; 1 -1];
elseif strcmp(syst_type,'exponential')
    rates = [1 0; 1 0];
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
trfactor = 1/10000000;

if fakey==1
    transfer_prob1 = [0 0; 2 2; 2 2; 2 2];
    transfer_prob2 = [0 0; 0 0; 0 0; 0 0];
    transfer_prob3 = [0 0; 0 0; 0 0; 0 0];
    transfer_prob4 = [4 0; 0 0; 0 0; 0 0];
    
    transfer_prob = cat(3,transfer_prob1',transfer_prob2',transfer_prob3',transfer_prob4').*trfactor;
else
    if strcmp(country,'USA') == 1
        load ../transferUSAnodes.mat
        load ../gflu51states.mat
        transfer_prob = symtransferUSA.*trfactor;
        % this is a relic...
    elseif strcmp(country,'USA') == 1
        load ../transferUSAnodes.mat
        transfer_prob = symtransferUSA.*trfactor;
    end
end

% temporary num_nodes = size(nodecodesUSA,1);
num_nodes = 60;

% initial state baseline
init_inf_percent = 1; % measure infected percentage-wise
initial_states = zeros(num_species,num_nodes);
% the total number of individuals using the hub is scaled according to the 
% connectivity of the hub
total_indiv = 5000; %round(1000.*nn./mean(nn(tophubsUSA)))+100; % 100 sets a minimum
init_inf = round(total_indiv*init_inf_percent/100); % set the number of infected
% temporary initial_states(species1,1,:) = total_indiv-init_inf; 
% temporary initial_states(species2,1,:) = init_inf; % assign number of infected
% a special infected hub can be selected

initial_states(:,:) = 0; % can't be zero for Gillespie to work
initial_states(:,ceil(num_nodes/2)) = 1000;

if singlecity==1
    transfer_prob = transfer_prob(:,:,infecthub);
end

% requiring a floor on population to avoid divide by zero
minpop = 1;
% time at which to start the transport
startfly_time = 0;
% total time in years
end_time = 1;
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
        tauleapfractional_wking(rates, end_time, initial_states, transfer_prob, transfer_int, tophubsUSA, minpop, startfly_time, cycle_period, cycle_phase, cycle_x, alpha_0, beta_0, googleflu, gflu51states, idState);
else
    error('invalid algorithm')
end
