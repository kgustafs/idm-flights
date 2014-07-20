% generate a Levy walk with stepsize probability r^-mu and t = r^-nu
% inner loop on each particle 1) step through a time grid with step size
% delta_t, 2) find the current time cur_t=i*delta_t, 3) while the Levy step
% time lev_t (initialized to zero) is less than cur_t 4) generate step_r and
% step_t, 5) set lev_t = lev_t + step_t, 6) init_r = final_r, 7) final_r =
% init_r + step_r, 8) end while, 9) linear interpolate between init_r and
% final_r to find step size during each delta_t

type = input('distribution type? ');

in_par = input('input parameters? (1 or 0) ');

delta_t = input('delta_t = ');
wrfrq = input('write frequency = ');

if in_par==1
    if strcmp(type,'powerlaw')==1
        mu = input('mu = ')
        kpareto = 1/(mu-1);
        thpareto = input('lower cutoff for power law = ')
        sigpareto = kpareto*thpareto;
        nu = input('nu = ')
        plscale = input('power law step scale = ')
    elseif strcmp(type,'normal')==1
        mean = input('mean = ');
        sigma = input('sigma = ');
    end
    Nt = input('Nt = ')
    Np = input('Np = ')
else
    if strcmp(type,'powerlaw')==1
        plscale = 1
        mu = 3
        kpareto = 1/(mu-1);
        thpareto = 0.01;
        sigpareto = kpareto*thpareto;
        nu = 1
    elseif strcmp(type,'normal')==1
        mean = 0
        sigma = 1
        nu = delta_t/100
    end
    Nt = 100
    Np = 10
end

pos_r = zeros(Np,Nt/wrfrq);

% time loop over the grid of steps

for j_p = 1:Np
    if mod(j_p,1)==0
        j_p
    end
    init_t = 0;
    final_t=0;
    final_r = 0;
    init_r = 0;
    
    step_t = 0;
    step_r = 0;
    cur_t = 0;
    
    pos_r_temp = 0;
    inc_step_r = 0;
    
    for i_t = 2:Nt
        
        cur_t = (i_t-1)*delta_t;
        while final_t <= cur_t
            if strcmp(type,'normal')==1
                step_r = random('norm',mean,sigma);
                step_t = random('exp',1000);
            elseif strcmp(type,'powerlaw')==1
                %step_r = (randht(1,'powerlaw',mu)-1).*sign(random('uniform',-1,1,1));
                step_r = (random('gp',kpareto,sigpareto,thpareto,1,1)./plscale).*sign(random('uniform',-1,1,1));
                step_t = abs(step_r)^(1/nu);
            end            
            init_t = final_t;
            final_t = init_t + step_t;
            init_r = final_r;
            final_r = init_r + step_r;
        end
        
        %inc_step_r = final_r*delta_t/lev_t;
        inc_step_r = step_r*(cur_t-init_t)/step_t;
        
        pos_r_temp=init_r+inc_step_r;
        
        
        if mod(i_t,wrfrq)==0
            pos_r(j_p,i_t/wrfrq) = pos_r_temp;
        end
    end
end



time = [0:delta_t*wrfrq:(Nt-1)*delta_t];


dpos_r= zeros(size(pos_r));
for k=1:size(pos_r,2)
    dpos_r(:,k) = pos_r(:,k)-pos_r(:,1);
end
stddpos_r= std(dpos_r,0,1);

save mu2p685walk_cut0p01.mat






