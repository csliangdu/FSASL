function t = ms2tex(mean_val, std_val, ismax, sigs, sigs2, prefix)
if ~exist('ismax', 'var')
    ismax = 1;
end

sigs(isnan(sigs)) = 0; % failed to reject a=b, nan means a = b, of course sigs = 0
sigs2(isnan(sigs2)) = 1.0;

n = length(mean_val);
t = prefix;
if ismax
    [~, idx] = max(mean_val);
else
    [~, idx] = min(mean_val);
end
for i1 = 1:n
    if isempty(sigs2)
        if i1 == idx
            t = [t, '& \tabincell{c}{ \textbf{', num2str(mean_val(i1) * 100, '%4.2f'), '} \\ \textbf{$\pm$ ', num2str(std_val(i1) * 100, '%4.2f'), '}} '];
        else
            t = [t, '& \tabincell{c}{ ', num2str(mean_val(i1) * 100, '%4.2f'), ' \\ $\pm$ ', num2str(std_val(i1) * 100, '%4.2f'), '} '];
        end
    else
        if sigs(i1) == 0
            t = [t, '& \tabincell{c}{ \textbf{', num2str(mean_val(i1) * 100, '%4.2f'), '} \\ \textbf{$\pm$ ', num2str(std_val(i1) * 100, '%4.2f'), ' } \\ \textbf{', num2str(sigs2(i1), '%4.2f'), ' }} '];
        else
            t = [t, '& \tabincell{c}{ ', num2str(mean_val(i1) * 100, '%4.2f'), ' \\ $\pm$ ', num2str(std_val(i1) * 100, '%4.2f'), ' \\', num2str(sigs2(i1), '%4.2f'), ' } '];
        end
    end
    
end
t = [t, '\\ \hline'];