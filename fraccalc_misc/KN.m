function y = KN(z,alpha,beta,theta)
s = sin(pi*(alpha-theta)/2);
c = cos(pi*(alpha-theta)/2);
y = (1/pi)*(z.^(alpha-1)*s./(1+2*c*z.^alpha+z.^(2*alpha)));