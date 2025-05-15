%% 1) Retrieve the pseudo publication calendar
% Publication Delay
Var_delay = [Var_delay;Varq_delay];
if sur ==1,   Var_delay = [Var_delay; Vars_delay]; end
if act ==1,   Var_delay = [Var_delay; Vara_delay]; end
if lab ==1,   Var_delay = [Var_delay; Varl_delay]; end
if mort ==1,   Var_delay = [Var_delay; Varmt_delay]; end

% Publication Groups
Var_pubgroup = [Var_pubgroup;Varq_pubgroup];
if sur ==1,   Var_pubgroup = [Var_pubgroup; Vars_pubgroup]; end
if act ==1,   Var_pubgroup = [Var_pubgroup; Vara_pubgroup]; end
if lab ==1,   Var_pubgroup = [Var_pubgroup; Varl_pubgroup]; end
if mort ==1,   Var_pubgroup = [Var_pubgroup; Varmt_pubgroup]; end



%  "input" structure which contains the relevant information in order to construct a pseudo real-time calendar as in the
% paper. Please note, that in the calendar generation code below, it is assumed that the final data publication refers to the quarterly variable coming out.

input.K = size(y_m,2); % number of higher frequency indicators
input.mismatch = mismatch; % mismatch in sampling frequency
input.mlags = monthvars; % number of months used for nowcasting
input.pubdelay = Var_delay; % vector of publication delays of dimension equal to number of higher frequency indicators.
input.pubseq = Var_pubgroup; % vector of groupings that define which variables come out in which order.
[puball groupall] = calendar_gen(input);
                       





