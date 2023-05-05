function [tp, tn, fp, fn] = simu_exp(simu, exp, tp, tn, fp, fn)
    % tp
    if simu < 0.00001 & exp == 'SL'
        tp = tp + 1;
    end
    % tn
    if simu > 0.00001 & exp == 'SS'
        tn = tn + 1;
    end
    % fp
    if simu < 0.00001 & exp == 'SS'
        fp = fp + 1;
    end
    % fn
    if simu > 0.00001 & exp == 'SL'
        fn = fn + 1;
    end
    fprintf('tp=%u\ntn=%u\nfp=%u\nfn=%u\n', tp, tn, fp, fn)
end