function [pubseq_midas groupall] = calendar_gen2(input)
%% Unpack information

K = input.K; % numger of higher frequency indicators
lagsm = input.lagsm; % number of months included
v  = input.v; % number of publications in higher frequency period
mstart = input.mstart; % Starting month to begin the nowcast cycle
pubdelay = input.pubdelay; % publication delays for all K
qmpub = input.qmpub; % number of months GDP comes out after the reference quarter
qseqpub = input.qseqpub; % publication number within the month GDP for the reference comes out.

%% Find total number of publications until GDP comes out
V = [];
%1) Find out how many of each publication types are within each point along
%the monthly cycle
[cnt_unique, unique_pub] = hist(pubdelay,unique(pubdelay));
pubtotal = [ unique_pub ; cnt_unique];
pubtotal = sortrows(pubtotal.',1,"descend").';
Vmax = v*(abs(mstart) + qseqpub + 1);
V = Vmax - (v-qseqpub+1); % adjustment for how many monthly publications are made before GDP comes out
cyclesmin = ceil(V/v); % min amount of cycles given lowest publication lag
cyclesmax = ceil(Vmax/v); % max amount of cycles given given GDP release
mismatch = input.mismatch;

% 2) Add intermittant publications (how many publications within a cycle have higher than 0 publication lag?)
        % i) Find out the publications delays for each publication within
        % cycle
        pubtypes = zeros(v,1);
        for i = 1:v
            if isempty(pubdelay(input(i).pubsec))
                pubtypes(i) = NaN;
            else
                pubtypes(i) = unique( input(1).pubdelay(input(i).pubsec));
            end
        end

        % ii) Add intermittant publications
        for i = 1:cyclesmax-cyclesmin
            pub_type_temp =  pubtotal(1,i+1:end);
        V = V + sum(sum(pubtypes(:) == pub_type_temp)) + sum(isnan(pubtypes)); % previously only nnz(pubtypes)
        end

% 3) Add number of GDP publications within time-frame
    % 1) Find months in which GDP comes out
    num_qmpubs = round(qmpub/mismatch) + round((abs(mstart)+1)/mismatch);

    V = V + num_qmpubs;

%% Add later the times when GDP lags get released
%{
% Add number of GDP publications within time-frame
    % 1) Find months in which GDP comes out
    num_qmpubs = round(qmpub/mismatch) + round((abs(mstart)+1)/mismatch);
%}

%% Create sequence for a given cycle 
pubseq = zeros(v,K);
for i = 1:v 
    pubseq(i,input(i).pubsec) = 1;
end

if size(pubseq,1)~=v
   error('Publication sequence must contain at exactly input.v number of rows.')
end

if size(pubdelay,2)~=K
    error('There needs to exactly as many publication delays be passed on as there are indicators.')
end
%% Generate full publication calender

% Create U-MIDAS availability
pubseq_midas = zeros(V,K*lagsm);
cycles = ceil(V/v);

% Loop for each publication cycle for all indicators

for j = 1:K
        pub_temp = zeros(V,lagsm);
        colstart = lagsm+mstart-pubdelay(j)+2;
        % Base availablity
        if colstart <= lagsm
            pub_temp(:,colstart:lagsm) = 1;
        end

for i = 1:cycles
    % Check if any updates are needed
     if   colstart-i>=1
for ii = 1:v
        % Add 1 if the variable comes out with the publication
        if sum(input(ii).pubsec == j)>0
            colind = colstart- i;
        else
            colind = colstart - i + 1;
        end
        if colind <= lagsm
       pub_temp(v*i + (ii-v):V,colind) = 1;
        end
end
    end 
end
pubseq_midas(:,j*lagsm-lagsm +1 : j*lagsm) = pub_temp;
end


%% Retrieve Group Structure
groupall = zeros(1,K*lagsm);

%{
% Get rid of potential empty publications
temp = [];
for i = 1:v
    if isempty(input(i).pubsec)
        inputtemp = NaN;
    else
        inputtemp = input(i).pubsec
    end
    temp = [temp inputtemp];
end

for i = 1:v-sum(isnan(temp))
    if  isempty(input(i).pubsec) ~= 1
    k_temp = unique(input(i).pubsec);
    for j = 1:length(k_temp)
        groupall(:,k_temp(j)*v : (k_temp(j)+1)*(v) -1) = i;
    end
    end
end
%}

for i = 1:K
    groupall(:,(i-1)*lagsm+1:lagsm*i) = i;
end








