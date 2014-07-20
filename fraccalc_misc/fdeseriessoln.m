%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% solution for nu=0.6 tanh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear alpha beta theta Np xmin xmax chi t tau x eta G Gs Gl K Ks Kl LL Picn

alpha      = input('alpha = ');
beta       = input('beta = ');
theta_extr = min(alpha,2-alpha)
theta      = input('theta = ');
Np         = 3000
xmin = input('xmin = '); % close to the end of the tail
xmax = input('xmax = '); % arbitrary
chi        = input('chi = ');
t          = input('t = ');
tau        = (chi^(1/beta)*t)^(-beta/alpha)
tau        = input('tau = ');


x_p = linspace(xmax/Np,xmax,Np);
x_n = linspace(xmin,xmin/Np,Np);
x = [x_n x_p];

eta_p = x_p.*tau;
eta_n = x_n.*tau;
eta = [eta_n eta_p];

if alpha==beta

    epsilon = input('epsilon = ');
    delta = (epsilon/2)./tau;
    %delta = (epsilon/2).*((chi^(1/beta).*t).^(-beta/alpha));
    LL_n = KN(-x_n,alpha,beta,-theta); % fde solution for x<0
    LL_p = KN(x_p,alpha,beta,theta); % x>0
    LL = [LL_n,LL_p];
    figure; semilogy(x,LL,'k'); grid; axis([xmin xmax min(LL) max(LL)]);
    title(['$\theta =$' num2str(theta) '$ \alpha = $' num2str(alpha) ...
        '$ \epsilon = $' num2str(epsilon) '$ \chi = $' num2str(chi)],'FontSize',20);
    if epsilon > 0
        Pic_n=zeros(1,Np); %initialize
        Pic_p=zeros(1,Np); %"
        for k = 1:size(t,2)
            k
            Ptemp2_n = []; % reset for new time point
            Ptemp2_p = [];
            for j = 1:Np
                etapdelta_p = (x_p(j)+epsilon/2)/tau(k); % set integration limits
                etamdelta_p = (x_p(j)-epsilon/2)/tau(k);
                etapdelta_n = (-x_n(j)+epsilon/2)/tau(k);
                etamdelta_n = (-x_n(j)-epsilon/2)/tau(k);
                I_n = quad(@(z)KN(z,alpha,beta,-theta),etamdelta_n,etapdelta_n); %integrate
                Ptemp_n = 1/(2*delta(k))/tau*I_n; % compute P(x,t)
                Ptemp2_n = [Ptemp2_n Ptemp_n]; % add to array
                I_p = quad(@(z)KN(z,alpha,beta,theta),etamdelta_p,etapdelta_p);
                Ptemp_p = 1/(2*delta(k))/tau*I_p;
                Ptemp2_p = [Ptemp2_p Ptemp_p];
            end
            Pic_n = [Pic_n; Ptemp2_n];
            Pic_p = [Pic_p; Ptemp2_p];
        end
        Pic = [Pic_n Pic_p];
        Picn=Pic(2,:)/(sum(Pic(2,:))*(x(2)-x(1))); %normalize to self
        %cc=max(PDFrfltr3r_Ax(:,936/52)/max(Picn/tau));
        figure; semilogy(x,Picn,'r'); %plot P(x,t)
        title(['$\theta =$' num2str(theta) '$ \alpha = $' num2str(alpha) ...
            '$ \epsilon = $' num2str(epsilon) '$ \chi = $' num2str(chi)],'FontSize',20);
        hold on; semilogy(bincntr_Ax,PDFrfltr3r_Ax(:,936/52),'k'); %plot data
        xlabel('$\delta x$'); ylabel('$P(\delta x)$');
    end

elseif alpha > beta

    preKs_p1 = 0;
    preKs_n1 = 0;
    preKs_p2 = 0;
    preKs_n2 = 0;

    kmax = 125

    % convergent representation for alpha > beta
    
    for k=0:kmax

        preKs_p1 = preKs_p1 + (gamma(1 - alpha*k)/gamma(1 - ...
            beta*k))*sin(k*pi*(theta - alpha)*0.5).*(-eta_p.^alpha).^k;

        preKs_n1 = preKs_n1 + (gamma(1 - alpha*k)/gamma(1 - ...
            beta*k))*sin(k*pi*((-theta) - alpha)*0.5).*(-(-eta_n).^alpha).^k;

        preKs_p2 = preKs_p2 + (gamma(1 - k/alpha)*gamma(1 + ...
            k/alpha))/(factorial(k)*gamma(1 - beta*k/alpha))*sin(k*pi*(theta - alpha)/(2*alpha)).*(-eta_p).^k;

        preKs_n2 = preKs_n2 + (gamma(1 - k/alpha)*gamma(1 + ...
            k/alpha))/(factorial(k)*gamma(1 - beta*k/alpha))*sin(k*pi*((-theta) - alpha)/(2*alpha)).*(-(-eta_n)).^k;

    end

    Ks_p = (preKs_p1 + preKs_p2)./(pi.*eta_p);
    Gs_p = tau.*Ks_p;
    Ks_n = (preKs_n1 + preKs_n2)./(pi.*(-eta_n));
    Gs_n = tau.*Ks_n;

    Ks = [Ks_n Ks_p];
    Gs = [Gs_n Gs_p];

    % large eta solution

    preKl_p = 0;
    preKl_n = 0;

    nmax =  20

    % asymptotic representation for alpha > beta
    
    for n=0:nmax

        preKl_p = preKl_p + (gamma(1 + alpha*n)/gamma(1 + ...
            beta*n))*sin(n*pi*0.5*(theta - alpha)).*(-eta_p.^(-alpha)).^n;

        preKl_n = preKl_n + (gamma(1 + alpha*n)/gamma(1 + ...
            beta*n))*sin(n*pi*0.5*((-theta) - alpha)).*(-(-eta_n).^(-alpha)).^n;

    end

    Kl_p = preKl_p./(pi.*eta_p);
    Gl_p = tau.*Kl_p;
    Kl_n = preKl_n./(pi.*(-eta_n));
    Gl_n = tau.*Kl_n;

    Kl = [Kl_n Kl_p];
    Gl = [Gl_n Gl_p];

    set(0,'DefaultTextInterpreter','latex');
    figure; semilogy(x,Gs,'k'); grid; axis([xmin xmax 0.001 1]);
    hold on; semilogy(x,Gl,'b');
    title(['$\theta =$ ' num2str(theta) '$ \alpha =$ ' num2str(alpha) ...
        '$ \beta =$ ' num2str(beta)],'FontSize',20);
    xlabel('$x$'), ylabel('$G$');

elseif alpha < beta
    
    preKs_p1 = 0;
    preKs_n1 = 0;
    preKs_p2 = 0;
    preKs_n2 = 0;

    kmax = 20

    % asymptotic representation for alpha < beta
    
    for k=0:kmax

        preKs_p1 = preKs_p1 + (gamma(1 - alpha*k)/gamma(1 - ...
            beta*k))*sin(k*pi*(theta - alpha)*0.5).*(-eta_p.^alpha).^k;

        preKs_n1 = preKs_n1 + (gamma(1 - alpha*k)/gamma(1 - ...
            beta*k))*sin(k*pi*((-theta) - alpha)*0.5).*(-(-eta_n).^alpha).^k;

        preKs_p2 = preKs_p2 + (gamma(1 - k/alpha)*gamma(1 + ...
            k/alpha))/(factorial(k)*gamma(1 - beta*k/alpha))*sin(k*pi*(theta - alpha)/(2*alpha)).*(-eta_p).^k;

        preKs_n2 = preKs_n2 + (gamma(1 - k/alpha)*gamma(1 + ...
            k/alpha))/(factorial(k)*gamma(1 - beta*k/alpha))*sin(k*pi*((-theta) - alpha)/(2*alpha)).*(-(-eta_n)).^k;

    end

    Ks_p = (preKs_p1 + preKs_p2)./(pi.*eta_p);
    Gs_p = tau.*Ks_p;
    Ks_n = (preKs_n1 + preKs_n2)./(pi.*(-eta_n));
    Gs_n = tau.*Ks_n;

    Ks = [Ks_n Ks_p];
    Gs = [Gs_n Gs_p];
    
    
    preKl_p = 0;
    preKl_n = 0;

    nmax =  125

    % convergent representation for alpha < beta
    
    for n=0:nmax

        preKl_p = preKl_p + (gamma(1 + alpha*n)/gamma(1 + ...
            beta*n))*sin(n*pi*0.5*(theta - alpha)).*(-eta_p.^(-alpha)).^n;

        preKl_n = preKl_n + (gamma(1 + alpha*n)/gamma(1 + ...
            beta*n))*sin(n*pi*0.5*((-theta) - alpha)).*(-(-eta_n).^(-alpha)).^n;

    end

    Kl_p = preKl_p./(pi.*eta_p);
    Gl_p = tau.*Kl_p;
    Kl_n = preKl_n./(pi.*(-eta_n));
    Gl_n = tau.*Kl_n;

    Kl = [Kl_n Kl_p];
    Gl = [Gl_n Gl_p];

    set(0,'DefaultTextInterpreter','latex');
    figure; semilogy(x,Gl,'b'); grid; axis([xmin xmax 0.001 1]);
    hold on; semilogy(x,Gs,'k');
    title(['$\theta =$ ' num2str(theta) '$ \alpha =$ ' num2str(alpha) ...
        '$ \beta =$ ' num2str(beta)],'FontSize',20);
    xlabel('$x$'), ylabel('$G$');

end
