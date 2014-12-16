function plot_result(dataset, candiAlgs, plot_flag)
%==========================setup=======================================
% dataset = 'jaffe_213n_676d_10c';
if ~exist('plot_flag', 'var')
    plot_flag = 1;
end

if ~exist('candiAlgs', 'var') || isempty(candiAlgs)
    candiAlgs = { 'LapScore', 'MCFS',  'LLCFS', 'UDFS', 'NDFS',  'SPFS', 'RUFS',  'JELSR_lpp', 'GLSPFS', 'FSSL_11_11_5'};
end
% candiAlgs = {'AllFea', 'LapScore'};
candiLineStyles = {'-', '-.',    '-', '-.',    '-', '-.',    '-', '-.',     '-', '--',   '-', '--',    '-', '--',    '-', '--',    '-', '--',    '-', '--',    '-', '--'};
candiMarkers = {'o', '+', 's', 'd',               'o', '+', 's', 'd',        'o', '+', 's', 'd',        'o', '+', 's', 'd',  };
candiColors = [0 0 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1;           0 0 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1;           0 0 0; 1 0 1; 0 1 1; 1 0 0; 0 1 0; 0 0 1];
candiMarkerSpacing = [5,5;5,5;5,5;5,5;5,5;5,5;         5,5;5,5;5,5;5,5;5,5;5,5;  5,5;5,5;5,5;5,5;5,5;5,5; ];

%=====================================================================
res_algs = [];

algs = {};
lineStyles = {};
markers = {};
colors = [];
markerSpacing = [];

ii = 0;
for idx = 1:length(candiAlgs )
    res_file = [dataset, '_best_result_', candiAlgs{idx}, '.mat'];

    if exist(res_file, 'file')
        ii = ii + 1;
        if exist(res_file, 'file'); load(res_file); end
        if exist('res_gs', 'var')
            if ii == 1;
                res_algs= res_gs;
            else
                fn = fieldnames(res_gs);
                for i2 = 1:length(fn)
                    if ~isfield(res_algs, (fn{i2}))
                        res_algs.(fn{i2}) = []; %place holder, should be removed, some algos did not record time2
                    end
                    res_algs.(fn{i2}) = [res_algs.(fn{i2}); res_gs.(fn{i2})];
                end
            end
           
            algs{end+1} = candiAlgs{idx};
            lineStyles{end+1} = candiLineStyles{ii};
            markers{end+1} = candiMarkers{ii};
            colors = [colors; candiColors(ii, :)];
            markerSpacing = [markerSpacing; candiMarkerSpacing(ii, :)];
        end
        clear res_gs;
    end
end

if ~isempty(res_algs)
    res_gs_tt = [];
    
    % res_gs_tt = compute_ttest(dataset, algs, length(FeaNumCandi));
if isvector(FeaNumCandi ) && length(FeaNumCandi) > 10
    tmp = find(FeaNumCandi(1:end-1) - FeaNumCandi(2:end) > 0);
    tmp = [1; tmp(:); length(FeaNumCandi)];
    ids = cell(length(tmp) - 1, 1);
    for i1 = 1:length(ids)
        ids{i1} = tmp(i1):tmp(i1+1);
    end
    message = compute_message(algs, res_algs, dataset, ids);
else 
    message = [];
end
    save(['res_algs_', dataset, '.mat'], 'algs', 'res_algs', 'res_gs_tt', 'message');
    fns = {'mean_acc', 'mean_nmi_sqrt', 'mean_nmi_max', 'mean_purity', 'mean_prec', 'mean_recall', 'mean_f1', ...
        'best_obj_acc', 'best_obj_nmi_sqrt', 'best_obj_nmi_max', 'best_obj_purity', 'best_obj_prec', 'best_obj_recall', 'best_obj_f1',...
        'jac', 'red','loocv'};
    if plot_flag
        xData = (1:length(FeaNumCandi));
        %     figure;
        for i1 = 1:length(fns)
            figure;
            my_prettyPlot(xData, res_algs.(fns{i1}), colors, lineStyles, markers, markerSpacing, dataset, '# of features', fns{i1}, algs, 'SouthWest');
            if strcmp(fns{i1},'loocv')
                my_prettyPlot(xData, 1-res_algs.(fns{i1}), colors, lineStyles, markers, markerSpacing, dataset, '# of features', fns{i1}, algs, 'SouthWest');
            end
        end
    end
%     my_prettyPlot(xData, res_algs.mean_nmi_max, colors, lineStyles, markers, markerSpacing, dataset, '# of features', 'Normalized Mutual Information', algs, 'SouthWest');
%     figure;
%     my_prettyPlot(xData, res_algs.red, colors, lineStyles, markers, markerSpacing, dataset, '# of features', 'Redundancy', algs, 'SouthWest');
%     figure;
%     my_prettyPlot(xData, res_algs.f1, colors, lineStyles, markers, markerSpacing, dataset, '# of features', 'JAC', algs, 'SouthWest');
end
end

function my_prettyPlot(xData, yData, colors, lineStyles, markers, markerSpacing, title, xlabel, ylabel, legends, legendLoc)


options.colors = colors;
options.lineStyles = lineStyles;
options.markers = markers;
% options.markerSpacing = markerSpacing;
options.title = title;
options.xlabel = xlabel;
options.ylabel = ylabel;
options.legendStr = legends;
options.legend = legends;
options.legendLoc = legendLoc;
options.xlimits = [1, length(xData)];
% options.ylimits = [min(yData(:)), max(yData(:)) ];
prettyPlot(xData,yData,options);
hold off;
end


function res_gs_tt = compute_ttest(dataset, candiAlgs, nFeaNumCandi)
res_gs_tt = cell(1, nFeaNumCandi);
fns = {'aio_acc', 'aio_nmi_max', 'aio_nmi_sqrt', 'aio_purity', 'aio_prec', 'aio_recall', 'aio_f1'};
fns2 = {'mean_acc', 'mean_nmi_max', 'mean_nmi_sqrt', 'mean_purity', 'mean_prec', 'mean_recall', 'mean_f1'};
for i1 = 1:nFeaNumCandi
    for i2 = 1:length(fns);
        res_gs_tt{1, i1}.([fns{i2}, '_tt']) = ones(length(candiAlgs)) * -1;
        res_gs_tt{1, i1}.([fns{i2}, '_tt_p']) = ones(length(candiAlgs)) * -1;
    end
end

for i1 = 1:length(candiAlgs);
    res_file = [dataset, 'best_result_', dataset, '_', candiAlgs{i1}, '.mat'];
    if exist(res_file, 'file') 
        if exist(res_file, 'file'); load(res_file); end
        if exist('res_aio', 'var')
            res1 = res_aio;
            res1_ps = res_gs_ps;
            clear res_aio res_gs_ps;
            
            for i2 = i1+1:length(candiAlgs);
                res_file = [dataset, 'best_result_', dataset, '_', candiAlgs{i1}, '.mat'];
                if exist(res_file, 'file')
                    if exist(res_file, 'file'); load(res_file); end
                   
                    if exist('res_aio', 'var')
                        res2 = res_aio;
                        res2_ps = res_gs_ps;
                        clear res_aio res_gs_ps;
                        
                        for i3 = 1:nFeaNumCandi
                            for i4 = 1:length(fns)
                                tmp1 = res1_ps.(fns2{i4});
                                b1_idx = tmp1(i3);
                                tmp2 = res2_ps.(fns2{i4});
                                b2_idx = tmp2(i3);
                                r1 = res1{b1_idx, i3}.(fns{i4});
                                r2 = res2{b2_idx, i3}.(fns{i4});
                                [t1, t2] = ttest(r1, r2);
                                tmp1 = res_gs_tt{1, i3}.([fns{i4},'_tt']);
                                tmp1(i1, i2) = t1;
                                res_gs_tt{1, i3}.([fns{i4},'_tt']) = tmp1;
                                tmp1 = res_gs_tt{1, i3}.([fns{i4},'_tt_p']);
                                tmp1(i1, i2) = t2;
                                res_gs_tt{1, i3}.([fns{i4},'_tt_p']) = tmp1;
                            end
                        end
                    end
                end
            end
        end
    end
end
end

function message = compute_message(algs, res_algs, dataset, ids)
message = [];
fns = {'best_obj_acc', 'best_obj_nmi_sqrt', 'mean_acc', 'mean_nmi_sqrt', 'loocv', 'jac', 'red'};
ismax = [1, 1, 1, 1, 0, 1, 0];

tex_header = '\begin{table*}';
tex_header = [tex_header, char(13),'\caption{', dataset, 'all results}'];
tex_header = [tex_header, char(13),'\tiny \centering \label{table:res_aio}'];

tex_align = '| c ';

tex_title = 'Data Sets';
for i1 = 1:length(algs)
    tex_title = [tex_title,  ' & ', algs{i1}];
    tex_align = [tex_align,  ' | ', 'c'];
end
tex_title = [tex_title '\\ \hline'];
tex_header = [tex_header, char(13),'\begin{tabular}{', tex_align, '| }'];
tex_header = [tex_header, char(13),'\toprule'];

for i1 = 1:length(ids)
    tmp = ids{i1};
    for i2 = 1:length(fns)
        sigs = zeros(size(res_algs.(fns{i2}), 1), 1);
        sigs2 = sigs;
        if ismax(i2)
            [~, best_id] = max(mean(res_algs.(fns{i2})(:, ids{i1}), 2));
        else
            [~, best_id] = min( mean(1 - res_algs.(fns{i2})(:, ids{i1}), 2));
        end
        for i3 = 1:length(sigs)
            [sigs(i3), sigs2(i3)] = ttest(res_algs.(fns{i2})(i3, ids{i1}), res_algs.(fns{i2})(best_id, ids{i1}));
        end
        if ismax(i2)
            message = [message, char(13), ms2tex(mean(res_algs.(fns{i2})(:, ids{i1}), 2), std(res_algs.(fns{i2})(:, ids{i1}), 0, 2), ismax(i2), sigs, sigs2, [dataset(1:5), '_', fns{i2}, '_', num2str(tmp(1)), '_', num2str(tmp(end)) ])];
        else
            message = [message, char(13), ms2tex(mean(1 - res_algs.(fns{i2})(:, ids{i1}), 2), std(1 - res_algs.(fns{i2})(:, ids{i1}), 0, 2), ismax(i2), sigs, sigs2, [dataset(1:5), '_', fns{i2}, '_', num2str(tmp(1)), '_', num2str(tmp(end)) ])];
        end
    end
end
message = [tex_title, char(13), message ];


tex_end = [];
tex_end = [tex_end, char(13), '\bottomrule' ];
tex_end = [tex_end, char(13), '\end{tabular}' ];
tex_end = [tex_end, char(13), '\end{table*}' ];

message = [tex_header, char(13), message, char(13), tex_end];
message = strrep(message, '_', '-');
fid=fopen([dataset, '.tex'], 'w+');
fprintf(fid,  '%s', message);
fclose(fid);
end