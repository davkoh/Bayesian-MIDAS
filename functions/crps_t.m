function [crps_t] = crps_t(y, df, location, scale)
    if nargin < 3
        location = 0;
    end
    if nargin < 4
        scale = 1;
    end

    y = y - location;
    crps_t = NaN(size(y));

    if scale == 1
        df(df <= 1) = NaN;
        bfrac = beta(0.5, df - 0.5) / (beta(0.5, 0.5 * df))^2;
        crps_t = y .* (2 * tcdf(y, df) - 1) + 2 ./ (df - 1) .* (tpdf(y, df) .* (df + y.^2) - sqrt(df) .* bfrac);
    else
        if (scale < 0)
            scale = NaN;
        end
        if (all(scale > 0, 'all'))
            crps_t = scale .* crps_t2(y ./ scale, df);
        else
            crps_t = scale .* crps_t2(y ./ scale, df);
            ind1 = df == Inf;
            ind2 = scale == 0;
            crps_t(ind1) = fill(scale * crps_norm(y ./ scale), 1, length(crps_t)).*ind1;
            crps_t(ind2) = fill(abs(y), 1, length(crps_t)).*ind2;
        end
    end
end