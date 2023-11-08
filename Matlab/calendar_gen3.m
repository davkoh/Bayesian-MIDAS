function [pubseq_midas groupall] = calendar_gen3(input)
%% Unpack information

%%% Next steps:
% check which information needs to be passed on
% make sure that the selection algorithm is robust (different subset
% selections
% Make sure that the user selection is intuitive

K = input.K; % numger of higher frequency indicators
mismatch = input.mismatch;
mlags = input.mlags;
mstart = input.mstart;
pubdelay = input.pubdelay; % publication delays for all K
pubseq = input.pubseq; % publication number within the month GDP for the reference comes out.
m_end = input.mend;
qvar_cutoff = find(unique(pubseq(find(abs(pubdelay) == pubdelay(1))))==pubseq(1));
mvar_cutoff = 1+qvar_cutoff;

%% Find the cutoff for the last month

%% Number of MIDAS sampled covariates
K_midas = mlags*K;

%% Generate full publication calender

% Retrieve months in which quarterly variables are updated
q_lag = pubdelay(qvar_cutoff); 
q_update = abs(q_lag)-mismatch:mismatch:m_end;
m_name = mstart:m_end;

% Total amount of months 
m_num = length(mstart:m_end); 

% Create U-MIDAS availability
pubseq_midas = zeros(1,K_midas);
pubseq_midas = [];

% Number of unique groups 
num_groups = length(unique(pubseq));
unique_groups = unique(pubseq);


% Find Starting availability
pubcal_start = zeros(1,K_midas);

for i = 1:K
    ind_start = i*mlags - mlags +1 - pubdelay(i+1) - mstart + 1;
    ind_end = i*mlags;

    pubcal_start(ind_start:ind_end) = 1;
end


for j = 1:m_num

      for  i = 1:num_groups

          if j == m_num && i == qvar_cutoff % If clause for the cutoff point
               pubseq_midas = [pubseq_midas;pubcal_start];
          elseif j == m_num && i>qvar_cutoff
              [];
          elseif j ~=m_num
                
              if i == pubseq(1) %mvar_cutoff % If clause: when GDP need not need updating

                  if sum(m_name(j) == q_update)>0 % If clause: when GDP need updating (then just copy the line from before
                   pubseq_midas = [pubseq_midas;pubcal_start];
                  end

              else


              
            m_ind_update = find(pubseq == unique_groups(i))-1; %find(pubseq == pubseq(i))-1; %find(pubseq == i)-1;

           for m = 1:length(m_ind_update)
               ind_beg = m_ind_update(m)*mlags-mlags+1;
               ind_end = m_ind_update(m)*mlags;

               % Find which columns to update 
               colupdate = min(6 + mstart + pubdelay(m_ind_update(m)+1) + j -1,mlags);

               % Update the columns 
               pubcal_start(1,ind_end-colupdate+1:ind_end) = 1;

           end

        if mlags + mstart + pubdelay(m_ind_update(m)+1) + j -1 > mlags
        else

         pubseq_midas = [pubseq_midas;pubcal_start];

        end

              end

        
          end
      end


end


%{
    for i = 1:K
        ind_beg = i*mlags-mlags+1;
        ind_end = i*mlags;

        if sum(pubseq_midas(j,ind_beg:ind_end)<mlags)

            pubseq_midas(ind_beg)

        end

    end

%}





for i = 1:K
    groupall(:,(i-1)*mlags+1:mlags*i) = i;
end

end








