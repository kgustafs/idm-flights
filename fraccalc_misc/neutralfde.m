%Neutral fractional diffusion script

alpha=1; % alpha = beta
theta_extr=min(alpha,2-alpha);
% input all of the parameters
theta=input('theta = (-1<0<1)    ')
%theta=-0.4;
eta=input('eta = (~14) ')
%eta=5;
shift=input('shift = (-1.5)     ')
norm=input('norm = (3.5)     ')
% define the domain in which to find the solution
Np=10000;
xc=5;
xx_p=linspace(-0.001,xc,Np);
% compute the solution from Mainardi et al Eq. 4.38
frac_p=(sin(pi*(alpha-theta)/2)/pi)*xx_p.^(alpha-1)./(1+2*cos(pi*(alpha-theta)/2)*xx_p.^alpha+xx_p.^(2*alpha));
xx_n=linspace(0,-xc, Np);
theta=-theta; % give theta opposite sign but use positive half of domain
frac_n=(sin(pi*(alpha-theta)/2)/pi)*xx_p.^(alpha-1)./(1+2*cos(pi*(alpha-theta)/2)*xx_p.^alpha+xx_p.^(2*alpha));
% concatenate domain and plot figure
xx=[xx_n xx_p];
frac=[frac_n frac_p];
set(0,'DefaultTextInterpreter','latex')
figure; semilogy(xx./eta-shift,norm.*eta.*frac,'b.');
title(['$\theta$ =' num2str(theta)],'FontSize',20);