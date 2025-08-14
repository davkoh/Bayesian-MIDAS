function[Xm] = transf_umidas(y_q,y_m,input)

Xm = NaN(size(y_q,1),input.mlags*size(y_m,2));
idx = 1:input.mlags:size(Xm,2);
ml = input.mlags/(input.mismatch)-1;
for j = 1:size(y_m,2)
    xtemp0 = [];
    for i = input.mismatch:-1:1
        xtemp0 = [xtemp0 y_m(i:input.mismatch:end,j)];
    end


    if ml>0
        xtemp1 = [];
        xtemp1 = [xtemp1 xtemp0(1+ml:end,:)];
        for i = 1:ml
            xtemp1 = [xtemp1 xtemp0(1+ml-i:end-i,:);];
        end
    else
        xtemp1 = xtemp0;
    end
    Xm(:,idx(j):idx(j)+input.mlags-1)= xtemp1;
end
end