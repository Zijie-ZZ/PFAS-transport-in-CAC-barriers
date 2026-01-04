function Exp = load_PFOS_data_from_mat(matPath)

S = load(matPath);

% --- helper to read time+data pairs from a cell table ---
    function [t, c] = read_pairs(CELL, rowStart, colPairs)
        t = []; c = [];
        for k = 1:size(colPairs,1)
            tk = CELL{rowStart:end, colPairs(k,1)};
            ck = CELL{rowStart:end, colPairs(k,2)};

            tk = tk(~isnan(tk));
            ck = ck(~isnan(ck));

            n = min(numel(tk), numel(ck));
            t = [t; tk(1:n)];
            c = [c; ck(1:n)];
        end
    end

pairs3 = [4 6; 11 13; 18 20];
pairs4 = [4 6; 11 13; 18 20; 25 27];

% Predefine struct array with consistent fields
Exp = struct('Q_mL_h', {}, 't_days', {}, 'C_C0', {});

[t12, c12] = read_pairs(S.PFOS_Q12_Plot, 4, pairs3);
[t24, c24] = read_pairs(S.PFOS_Q24_Plot, 3, pairs4);
[t36, c36] = read_pairs(S.PFOS_Q36_Plot, 4, pairs3);

Exp(1) = struct('Q_mL_h',12,'t_days',t12,'C_C0',c12);
Exp(2) = struct('Q_mL_h',24,'t_days',t24,'C_C0',c24);
Exp(3) = struct('Q_mL_h',36,'t_days',t36,'C_C0',c36);

end
