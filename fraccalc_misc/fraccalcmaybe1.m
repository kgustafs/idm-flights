alpha = input('alpha = ');
beta = alpha; % neutral FDE solution
theta = input('theta = ');
xmin = -1000; % close to the end of the tail
xmax = 2; % arbitrary
Np = 2000; % arbitrary, but good enough
xx_p = linspace(xmax/Np,xmax,Np);
xx_n = linspace(xmin,-xmax/Np,Np);
xx = [xx_n xx_p];
t = 936; %arbitrary but reasonable
chi = input('chi = ');
epsilon = input('epsilon = ');
tau = chi^(1/beta).*t; % scaling from x to eta
delta = (epsilon/2)./tau;
%delta = (epsilon/2).*((chi^(1/beta).*t).^(-beta/alpha));
LL_n = KN(-xx_n,alpha,beta,-theta); % fde solution for x<0
LL_p = KN(xx_p,alpha,beta,theta); % x>0
LL = [LL_n,LL_p]; 
Pic_n=zeros(1,Np); %initialize
Pic_p=zeros(1,Np); %"
for k = 1:size(t,2)
    k
    Ptemp2_n = []; % reset for new time point
    Ptemp2_p = [];
    for j = 1:Np
        etapdelta_p = (xx_p(j)+epsilon/2)/tau(k); % set integration limits
        etamdelta_p = (xx_p(j)-epsilon/2)/tau(k);
        etapdelta_n = (-xx_n(j)+epsilon/2)/tau(k);
        etamdelta_n = (-xx_n(j)-epsilon/2)/tau(k);
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
Picn=Pic(2,:)/(sum(Pic(2,:))*(xx(2)-xx(1))); %normalize to self
%cc=max(PDFrfltr3r_Ax(:,936/52)/max(Picn/tau));
figure; semilogy(xx,Picn,'r'); %plot P(x,t)
title(['$\theta =$' num2str(theta) '$ \alpha = $' num2str(alpha) ...
    '$ \epsilon = $' num2str(epsilon) '$ \chi = $' num2str(chi)],'FontSize',20);
hold on; semilogy(bincntr_Ax,PDFrfltr3r_Ax(:,936/52),'k'); %plot data
xlabel('$\delta x$'); ylabel('$P(\delta x)$');