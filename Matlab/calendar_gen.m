function [pubseq_midas groupall] = calendar_gen(input)
%% Unpack information

K = input.K; % numger of higher frequency indicators
mismatch = input.mismatch; % mismatch in sampling frequency
mlags = input.mlags; % number of monthly observations used for estimation
mstart = input.mstart; % month index when nowcasting starts
m_end = input.mend; % month index when nowcasting ends
pubdelay = input.pubdelay; % publication delays for all K indicators
pubseq = input.pubseq; % publication number within the month GDP for the reference comes out.


qvar_cutoff = find(unique(pubseq(find((pubdelay) == pubdelay(1)))) == pubseq(1)); % publications in final month, given the monthly and quarterly publication delays


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

% Find Starting availability
pubcal_start = zeros(1,K_midas);

for i = 1:K
    ind_start = i*mlags - mlags +1 - pubdelay(i+1) - mstart + 1;
    ind_end = i*mlags;

    pubcal_start(ind_start:ind_end) = 1;
end


for j = 1:m_num

      for  i = unique(pubseq)'

          if j == m_num && i == pubseq(qvar_cutoff) % If clause for the cutoff point %%% Needs to change with the new loop definition
               pubseq_midas = [pubseq_midas;pubcal_start];
          elseif j == m_num && i>pubseq(qvar_cutoff)
              [];
          elseif j ~=m_num
                
              if i == pubseq(1) %mvar_cutoff % If clause: when GDP need not need updating

                  if sum(m_name(j) == q_update)>0 % If clause: when GDP need updating (then just copy the line from before
                   pubseq_midas = [pubseq_midas;pubcal_start];
                  end

              else


              
            m_ind_update = find(pubseq == i)-1; %find(pubseq == unique_groups(i))-1; %find(pubseq == pubseq(i))-1; %find(pubseq == i)-1;

           for m = 1:length(m_ind_update)
               ind_beg = m_ind_update(m)*mlags-mlags+1;
               ind_end = m_ind_update(m)*mlags;

               % Find which columns to update 
               colupdate = min(mlags + mstart + pubdelay(m_ind_update(m)+1) + j -1,mlags);

               % Update the columns 
               pubcal_start(1,ind_end-colupdate+1:ind_end) = 1;

           end

        if mlags + mstart + min(pubdelay(m_ind_update+1)) + j -1 > mlags % We assume here that those variables published together have the same publication lag. Ok, so this is really not to update anymore after the columns are full
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








