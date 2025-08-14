function [C,Xall] = midas_dat_r2_final(X,grouping,polyt)
grall = unique(grouping); % unique groups due to lags of individual monthly variables.
r = 2; % number of restrictions, for this function, 2.
poly = polyt - r; % polynomial degree adjusted for restrictions
lagl = histc(grouping, grall); % lag length of each covariate;
Call = [];%zeros(size(X,1),max(grall)*poly);
Xall = [];

% This is the correctly translated weights creator
for m = 1:size(grall,1)
        idx = find(grouping == grall(m));
        Xt =  X(:,idx);
        
        if lagl(m)>r
        C = zeros(size(Xt,1),poly);
    for j = 1:poly
        Ct= [];
    for i = 0:lagl(m)-1
       
        w = i^(j+1) + j*(lagl(m)-1)^(j+1) - (j+1)*(lagl(m)-1)^j*(i);
        Xtt = Xt(:,i+1)*w;
         Ct = [Ct Xtt];          
    end
          C(:,j) = sum(Ct,2);
    end
    
  
    

    else
        C = Xt;
end
    Call = [Call C];
    Xall = [Xall Xt];
end
C = Call;
