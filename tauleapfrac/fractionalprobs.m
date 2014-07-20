clear all

% this file contains the deparature, arrival, latitude, longitude,
% distance, seats, end dates, days per week, sorted by departure IATA code
load /Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/countries/USA/data/usadata.mat
load /Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/countries/USA/data/stateIATAUSA.mat
load /Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/countries/USA/data/statepostal.mat


min_connects = 1; % finding intersect with Code list acts as a size filter

% count_unique is an add-on function
[uniqueIATA,nn] = count_unique(usaorigin);
tophubsUSA = find(nn>15000);
nodeidUSA = find(nn>min_connects);
IATAnodes = uniqueIATA(nodeidUSA);

[nodecodesUSA, idnodes, idCode] = intersect(IATAnodes,Code);

% initialize the transfer matrix for all USA airports
transferUSA = zeros(2,size(nodecodesUSA,1),size(nodecodesUSA,1));

SRtransferratio = 1;

for ii = 1:size(nodecodesUSA,1)
    ii
    for jj = 1:size(nodecodesUSA,1)
        arrival = find(strcmp(usadestination,nodecodesUSA(jj))==1);
        departure = find(strcmp(usaorigin,nodecodesUSA(ii))==1);
        % get the index, among all Arrival and Departure Airports, sorted
        % alphabetically by departure
        connection = intersect(arrival,departure);
        transferUSA(1,ii,jj) = sum(usaseats(connection).*usadays(connection));
    end
    %totalii = sum(transferUSA(1,ii,:),3);
    %transferUSA(:,ii,:) = transferUSA(:,ii,:)./totalii;
    % transfer probability for infected, default SRtransferratio = 1
    transferUSA(2,ii,:) = transferUSA(1,ii,:).*SRtransferratio;
end

uppertri1 = triu(squeeze(transferUSA(1,:,:)));
uppertri2 = triu(squeeze(transferUSA(2,:,:)));

lowertri1 = tril(squeeze(transferUSA(1,:,:)));
lowertri2 = tril(squeeze(transferUSA(2,:,:)));

avgtri1 = (uppertri1+lowertri1)./2;
avgtri2 = (uppertri2+lowertri2)./2;

symtransferUSA(1,:,:) = avgtri1+avgtri1.';
symtransferUSA(2,:,:) = avgtri2+avgtri2.';

narc = squeeze(symtransferUSA(1,:,:));

%%
Dalpha = 1;
alphafrac = 0.5;
h = 1;
center = 1;
muc = 1/(2*abs(cos(alphafrac*pi/2)));
ffmu = (Dalpha/h^alphafrac)*muc;

tau = 1/(12*sqrt(2));

vectorlength =60;

if rem(vectorlength,2) == 1
    lBigG = (vectorlength+1)/2;
else
    lBigG = vectorlength/2;
end

nsteps = 50;
nvector = zeros(vectorlength,nsteps+1);
nvector(lBigG,1) = 1000;
%nvector(:,1) = nvector(:,1)./sum(nvector(:,1));
vvector = zeros(vectorlength,1);
kvector = zeros(vectorlength,1);

time = 0;

%nvector(:,1) = nvector(:,1)/sum(nvector(:,1));

for t=1:nsteps
    fmu = ffmu*tau;
    t
    tau
    BigG = zeros(lBigG,1);
    %BigG2 = zeros(vectorlength,1);
    vvector = zeros(vectorlength,1);
    for jjj = 1:vectorlength
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
        end
    end
    %     if rem(vectorlength,2) == 1
    % here, the kernel is being doubled in width to account for
    % shifting at each value of lambda
    BigG2 = [BigG(end:-1:1); BigG(2:end)];
    %     else
    %         BigG2 = [BigG(end:-1:1); BigG(2:end); 0];
    %     end
    %BigG2 = BigG2./sum(BigG2); % must be renormalized over the two-tailed kernel
    vvector = zeros(vectorlength,1);
    for lambda=1:vectorlength
        BigG2rc = BigG2(vectorlength-lambda+1:end-lambda+1);
        BigG2rc = BigG2rc./sum(BigG2rc); % renormalization: fair?
        %kvector = mnrnd(nvector(lambda,t),BigG2rc); % random sampling
        kvector = nvector(lambda,t).*BigG2rc; kvector = kvector'; % expectation of the kernel
        %kvector(end-lambda+1:end) = kvector(1:lambda);
        %kvector(lambda) = kvector(lambda) - nvector(lambda,t);
        vvector(lambda) = vvector(lambda) - sum(kvector(1:end~=lambda));
        vvector(1:end~=lambda) = vvector(1:end~=lambda) + kvector(1:end~=lambda)';
    end
    nvector(:,t+1) = nvector(:,t) + vvector;
    %nvector(:,t+1) = nvector(:,t+1)./sum(nvector(:,t+1));
    time = [time t*tau];
end

% BigGm = BigG(end:-1:2);
% figure; semilogy([-49:1:49],[BigGm;BigG]./sum([BigGm;BigG]));
% do multinomial simulation - divide each cell by the initial number of
% particles just after the first time step,
% that will look like the kernel for each timestep