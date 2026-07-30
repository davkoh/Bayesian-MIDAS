%% ===================================================================
%%  Helper: significance stars
%% ===================================================================
function s = sig_stars(ratio, pv)
    if ratio >= 1 || isnan(pv)
        s = '';
    elseif pv < 0.01
        s = '***';
    elseif pv < 0.05
        s = '**';
    elseif pv < 0.10
        s = '*';
    else
        s = '';
    end
end
