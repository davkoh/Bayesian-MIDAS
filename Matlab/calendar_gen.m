function output = calendar_gen(input)
%% Unpack information

K = input.K; % numger of higher frequency indicators
lagsm = input.lagsm; % number of months included
V = input.V; % number of nowcastperiods
v  = input.v; % number of publications in higher frequency period
mstart = input.mstart; % Starting month to begin the nowcast cycle
pubdelay = input.pubdelay; % publication delays for all K

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

% Cycles (defines how many higher frequency time units pass within the whole nowcast cycle)
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

output = pubseq_midas;




