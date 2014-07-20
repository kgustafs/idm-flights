function frac_kern = frac_kernel(tau,ffmu,h,num_nodes,alphafrac,center,Dalpha)
%FRAC_KERNEL
%   Detailed explanation goes here

fmu = ffmu*tau;
BigG = zeros(num_nodes,1);
frac_kern = zeros(num_nodes,1);
for jjj = 1:num_nodes
    if alphafrac>0 && alphafrac<1
        if fmu > 0.5
            tau = abs(cos(pi*alphafrac/2))*h^alphafrac/Dalpha;
        end
        if center == jjj
            BigG(jjj) = 1 - 2.*fmu;
        else
            BigG(jjj) = fmu*abs(gamma(alphafrac+1)/(gamma(jjj+1)*gamma(alphafrac-jjj+1)));
        end
    elseif alphafrac==1
        if Dalpha*tau > 1
            tau = 0.999;
        end
        BigG(jjj) = (acot(Dalpha*tau/(abs(center-jjj)*h+h/2)) - ...
            acot(Dalpha*tau/(abs(center-jjj)*h-h/2)))/pi;
    elseif alphafrac>1 && alphafrac<=2
        if fmu > 0.5/alphafrac
            tau = abs(cos(pi*alphafrac/2))*h^alphafrac/(Dalpha*alphafrac);
        end
        if center == jjj
            BigG(jjj) = 1 - 2.*alphafrac.*fmu;
        elseif abs(jjj-center) == 1
            % alpha choose 2
            BigG(jjj) = fmu*(1+gamma(alphafrac+1)/(gamma(3)*gamma(alphafrac-1)));
        else
            % alpha choose (j+1)
            BigG(jjj) = fmu*abs(gamma(alphafrac+1)/(gamma(jjj+2)*gamma(alphafrac-jjj)));
        end
    else
        alphafrac
        return
%         if rem(vectorlength,2) == 1
%             lBigG = (vectorlength+1)/2;
%         else
%             lBigG = vectorlength/2;
%         end
    end
end
%    if rem(num_nodes,2) == 1
frac_kern = [BigG(end:-1:1); BigG(2:end)];
%     else
%         BigG2 = [BigG(end:-1:1); BigG(2:end); 0];
%     end
%frac_kern = BigG2./sum(BigG2); % must be renormalized over the two-tailed kernel

end

