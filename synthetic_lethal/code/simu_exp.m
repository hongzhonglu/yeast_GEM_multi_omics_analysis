function [tp, tn, fp, fn] = simu_exp(simu, exp, tp, tn, fp, fn)
    % tp
    if simu < 0.000001 & (exp == 'SL' | ((-1000 <= exp) & (exp <= -0.35)))
        tp = tp + 1;
    end
    % tn
    if simu > 0.000001 & (exp == 'SS' | exp > -0.35) 
        tn = tn + 1;
    end
    % fp
    if simu < 0.000001 & (exp == 'SS' | exp > -0.35)
        fp = fp + 1;
    end
    % fn
    if simu > 0.000001 & (exp == 'SL' | ((-1000 <= exp) & (exp <= -0.35)))
        fn = fn + 1;
    end
    %fprintf('tp=%u\ntn=%u\nfp=%u\nfn=%u\n', tp, tn, fp, fn)
end