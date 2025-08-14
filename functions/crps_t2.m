function [crps_t2] = crps_t2(y, df)
    
    crps_t2 = NaN(size(y));

        df(df <= 1) = NaN;
        bfrac = beta(0.5, df - 0.5) / (beta(0.5, 0.5 * df))^2;
        crps_t2 = y .* (2 * tcdf(y, df) - 1) + 2 ./ (df - 1) .* (tpdf(y, df) .* (df + y.^2) - sqrt(df) .* bfrac);
    
end