function a = propensity(input_counts)

% returns propensities for present state
% species, from Schmelzer, Burrage, Carletti

global species1 species2 syst_type beta alpha

total_indiv = sum(input_counts);
    
if strcmp(syst_type,'SIR')
           
    c = [alpha/total_indiv beta];
    
    input_counts = input_counts';
    
    a(1) = c(1)*input_counts(species1)*input_counts(species2);
    a(2) = c(2)*input_counts(species2);
    
    a=[a(1) a(2)]';
    
elseif strcmp(syst_type,'SIS')
    % apparently the same propensity vector can be used for SIS as for SIR
    % just need to change the stochiometry in SIR_flights.m
    
    c = [alpha/total_indiv beta];
    
    input_counts = input_counts';
    
    a(1) = c(1)*input_counts(species1)*input_counts(species2);
    a(2) = c(2)*input_counts(species2);
    
    a=[a(1) a(2)]';
elseif strcmp(syst_type,'exponential')
        
    c = [alpha/total_indiv]; % where alpha is rho for growth rate
    
    input_counts = input_counts';
    
    a(1) = c(1)*input_counts(species1);
    a(2) = c(1)*input_counts(species2);
    
    a=[a(1) a(2)]';
else
    error('invalid system type')
end