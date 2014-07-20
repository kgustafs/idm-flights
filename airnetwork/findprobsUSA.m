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

idState = [];
for k = 1:size(idCode)
    idState(k) = find(strcmp(statepostal,State(idCode(k))));
end

% for a small sample of the largest hubs

[idatl,atl] = find(strcmp(uniqueIATA(tophubsUSA(1)),usaorigin)==1);
[idden,den] = find(strcmp(uniqueIATA(tophubsUSA(2)),usaorigin)==1); 
[idewr,ewr] = find(strcmp(uniqueIATA(tophubsUSA(3)),usaorigin)==1);
[idiah,iah] = find(strcmp(uniqueIATA(tophubsUSA(4)),usaorigin)==1);
[idlax,lax] = find(strcmp(uniqueIATA(tophubsUSA(5)),usaorigin)==1);
[idord,ord] = find(strcmp(uniqueIATA(tophubsUSA(6)),usaorigin)==1);
[idsfo,sfo] = find(strcmp(uniqueIATA(tophubsUSA(7)),usaorigin)==1);

septfois = ['idatl';'idden';'idewr';'idiah';'idlax';'idord';'idsfo'];

% for kkk = 1:size(EndDate,1)
%     if strcmp(EndDate(kkk),'2010-12-31')==0
%         usaseats(kkk)=0;
%     end
% end

% compute the transfer matrix for the six largest hubs
transferUSA7 = zeros(2,7,7);

for ii = 1:size(septfois,1)
    for jj = 1:size(septfois,1)
        transferUSA7(1,ii,jj) = sum(strcmp(usadestination(eval(septfois(ii,:))),uniqueIATA(tophubsUSA(jj))).*usaseats(eval(septfois(ii,:))).*usadays(eval(septfois(ii,:))));
        transferUSA7(2,ii,jj) = sum(strcmp(usadestination(eval(septfois(ii,:))),uniqueIATA(tophubsUSA(jj))).*usaseats(eval(septfois(ii,:))).*usadays(eval(septfois(ii,:))));
    end
    transferUSA7(:,ii,:) = transferUSA7(:,ii,:)./sum(transferUSA7(1,ii,:),3);
end

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

save /Users/gustafsonkyle/Research/EMOD/Basil1/flights_network/countries/USA/data/transferUSAnodes.mat transferUSA symtransferUSA tophubsUSA nodecodesUSA
