% ivp_fit_012908.m

clear
%
% Fractional diffusion parameters 
% neutral asymmetric case with tail on left
%
alpha=input('alpha =   ');
abs_theta_max=min(alpha,2-alpha);
if alpha <1
    theta= abs_theta_max;
else
    theta= -abs_theta_max;
end
%
% Initial value parameters
%
epsilon_ic=input('epsilon  =    ');
time=936;   % This corresponds to the first frame of the data
chi=input('chi  (=0.175))=    ');
tau=time*(chi^(1/alpha))
%tau=input('tau =    ');
%
% Rescaled varaiables for integral
%
epsilon=epsilon_ic/tau;     
Ng=2000; % half of total # of grid points
xmin=-50; xmax=2; cero=10^-3;
delta_x=(xmax-xmin)/(2*Ng);
i_epsilon=floor(0.5*epsilon/delta_x);
%
% Evaluation of integrand
%
% x>0 
%
xp=linspace(cero,xmax,Ng);
s=sin(pi*(alpha-theta)/2);
c=cos(pi*(alpha-theta)/2);
Lp=(1/pi)*(xp.^(alpha-1)*s./(1+2*c*xp.^alpha+xp.^(2*alpha)));
%
% x<0 
%
xn=linspace(xmin,-cero,Ng);
s=sin(pi*(alpha+theta)/2);
c=cos(pi*(alpha+theta)/2);
Ln=(1/pi)*((-xn).^(alpha-1)*s./(1+2*c*(-xn).^alpha+(-xn).^(2*alpha)));

x=[xn,xp];
L=[Ln,Lp];
%
% Integral computation using a aimple Riemann sum
%
ik=0;
for ii=1+i_epsilon : 2*Ng-i_epsilon,
    L_ic(ii-i_epsilon)=(1/epsilon)*delta_x*sum(L(ii-i_epsilon:ii+i_epsilon));
    x_ic(ii-i_epsilon)=x(ii);
end
%
% Normalize distribution 
%
Ln_ic=L_ic/(sum(L_ic)*(x_ic(2)-x_ic(1)));
%
% Plot in physical varaibles
%
figure(10)
semilogy(tau*x_ic,Ln_ic/tau,'r')
hold on
axis([-800 100 10^-5 10^-2])

figure(20)
loglog(abs(tau*x_ic),Ln_ic/tau,'r')
hold on
%axis([-800 100 10^-5 10^-2])

%
% Plot fit
%
load PDFonly
pdf_x_t=PDFrossbyfltr3r_Ax; % Nx=334; Nt=101
x_bin=bincntrossbyfltr3r_Ax;
[Nx Nt]=size(pdf_x_t);
time=linspace(0,(Nt-1)*52,Nt);  % 52=scale factor to go to physical time
iframes=[19, 29, 39, 49];
    
figure(10);
semilogy(x_bin,pdf_x_t(:,iframes(1)))
hold on
xlabel('x'); ylabel('P(x,t)')
title(['Raw pdf at t=',num2str(time(iframes))])

figure(20);
loglog(abs(x_bin),pdf_x_t(:,iframes(1)))
xlabel('x'); ylabel('P(x,t)')
title(['Raw pdf at t=',num2str(time(iframes))])
hold on

NN=size(iframes);
for ip=1:NN(2)-1,
    figure(10)
    semilogy(x_bin,pdf_x_t(:,iframes(ip+1)))
    figure(20);
    loglog(abs(x_bin),pdf_x_t(:,iframes(ip+1)))
end
figure(10)
axis([-2000 100 10^-5 10^-2])
figure(20)
axis([10 10^3 10^-4 10^-2])

figure(30)
subplot(1,2,1)
semilogy(x_bin,pdf_x_t(:,iframes(1)))
hold on
semilogy(tau*x_ic,Ln_ic/tau,'r')
axis([-800 100 10^-4 10^-2])
xlabel('x'); ylabel('P(x,t_0)')

subplot(1,2,2)
loglog(abs(x_bin),pdf_x_t(:,iframes(1)))
hold on
loglog(abs(tau*x_ic),Ln_ic/tau,'r')
xlabel('x'); ylabel('P(x,t_0)')
axis([10 1000 10^-4 0.01])
%
% Rescale maximum to have the same normalization
%
cc=max(pdf_x_t(:,iframes(1))/max(Ln_ic/tau));
figure(40)
subplot(1,2,1)
semilogy(x_bin,pdf_x_t(:,iframes(1)))
hold on
semilogy(tau*x_ic,cc*Ln_ic/tau,'r')
axis([-800 100 10^-4 10^-2])
xlabel('x'); ylabel('P(x,t_0)')
title(['\alpha=\beta= ', num2str(alpha)]);

subplot(1,2,2)
loglog(abs(x_bin),pdf_x_t(:,iframes(1)))
hold on
loglog(abs(tau*x_ic),cc*Ln_ic/tau,'r')
xlabel('x'); ylabel('P(x,t_0)')
axis([10 1000 10^-4 0.01])
title(['\chi= ', num2str(chi), '        \epsilon= ', num2str(epsilon)]);
