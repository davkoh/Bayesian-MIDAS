%% Change to data handling function

input.pubseq(input.pubseq==6) = 5;
input.pubseq(17) = 6; 
input.pubseq = input.pubseq(1:17);


input.K = size(y_m,2)-1; % number of higher frequency indicators
input.mismatch = mismatch; % mismatch in sampling frequency
input.mlags = monthvars; % number of months used for nowcasting
input.pubdelay = Var_delay(1:17); % vector of publication delays of dimension equal to number of higher frequency indicators.
[puball groupall] = calendar_gen(input);
Xm = Xm(:,1:end-6);
