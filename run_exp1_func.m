function [FeaNumCandi_aio, res_gs_aio, res_aio_aio, res_gs_ps_aio] = run_exp1_func(datasets, candiAlgs, username, password)

[flag_writeable, flag_uploadable, prefix] = mdcs_check(username, password);

if ~exist('datasets', 'var') || isempty(datasets)
    % datasets = {'test'};
    datasets = {'medical_706n_1449d_17c', 'PIE_Pose27_1428n_1024d_68c', 'USPS49_1673n_256d_2c', 'mfeat_pix_2000n_240d_10c'};
end
if ischar(datasets); datasets = {datasets}; end

if ~exist('candiAlgs', 'var') || isempty(candiAlgs)
    candiAlgs = {'AllFea', 'MaxVar', 'LapScore', 'TraceRatio', 'SPEC', 'LLCFS', 'SPFS', 'MCFS', 'UDFS', 'NDFS', 'RUFS',  'JELSR', 'GLSPFS', 'FSSL'};
end
if ischar(candiAlgs); candiAlgs = {candiAlgs}; end

if ~exist('exp_settings', 'var');  exp_settings = []; end
if ~isfield(exp_settings, 'FeaNumCandi')
    exp_settings.FeaNumCandi = [[5:5:50],[10:10:150],[50:50:300]];
end
if ~isfield(exp_settings, 'nKmeans')
    exp_settings.nKmeans = 20;
end
if ~isfield(exp_settings, 'prefix_mdcs')
    exp_settings.prefix_mdcs = prefix;
end

FeaNumCandi = exp_settings.FeaNumCandi;
FeaNumCandi_aio = cell(length(datasets), length(candiAlgs));
res_gs_aio = cell(length(datasets), length(candiAlgs));
res_aio_aio = cell(length(datasets), length(candiAlgs));
res_gs_ps_aio = cell(length(datasets), length(candiAlgs));

root_dir = pwd;
addpath(root_dir);
for id = 1:length(datasets)
    dataset = datasets{id};
    X = extractXY(dataset);
    exp_settings.FeaNumCandi = FeaNumCandi(FeaNumCandi < size(X, 2));
    clear X;
    
    disp(['data = ', dataset, ' ...']);
    try
        if ~exist([prefix, filesep, dataset], 'dir')
            mkdir([prefix, filesep, dataset]);
        end
        exp_settings.prefix_mdcs = [prefix, filesep, dataset];
    catch
        disp(['create dir: ', [prefix, filesep, dataset], 'failed, check the authorization']);
    end
    
    for iAlg = 1:length(candiAlgs)
        algo = candiAlgs{iAlg};
        disp(['algo = ', algo, ' ...']);
        switch lower(algo)
            case lower('AllFea')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_allfea_single_func(dataset, exp_settings);
            case lower('MaxVar')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_maxvar_single_func(dataset, exp_settings);
            case lower('LapScore')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_lapscore_single_func(dataset, exp_settings);
            case lower('SPEC')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_spec_single_func(dataset, exp_settings);
            case lower('TraceRatio')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_traceratio_single_func(dataset, exp_settings);
            case lower('LLCFS')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_llcfs_single_func(dataset, exp_settings);
            case lower('UDFS')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_udfs_single_func(dataset, exp_settings);
            case lower('SPFS')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_spfs_single_func(dataset, exp_settings);
            case lower('MCFS')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_mcfs_single_func(dataset, exp_settings);
            case lower('NDFS')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_ndfs_single_func(dataset, exp_settings);
            case lower('RUFS')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_rufs_single_func(dataset, exp_settings);
            case lower('JELSR_lpp')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_jelsr_lpp_single_func(dataset, exp_settings);
            case lower('JELSR_lle')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_jelsr_lle_single_func(dataset, exp_settings);
            case lower('JELSR_liang_lpp')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_jelsr_liang_lpp_single_func(dataset, exp_settings);
            case lower('JELSR_liang_lle')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_jelsr_liang_lle_single_func(dataset, exp_settings);
            case lower('CGSSL')
                disp('not supported yet');
            case lower('GLSPFS')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_glspfs_single_func(dataset, exp_settings);
            case lower('FSSL_11_11_1')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_fsasl_11_11_1_single_func(dataset, exp_settings);
            case lower('FSSL_11_11_5')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_fsasl_11_11_5_single_func(dataset, exp_settings);
             case lower('FSSL_11_5_5')
                [FeaNumCandi_aio{id, iAlg}, res_gs_aio{id, iAlg}, res_aio_aio{id, iAlg}, res_gs_ps_aio{id, iAlg}] = fs_unsup_fsasl_11_5_5_single_func(dataset, exp_settings);
            otherwise
                disp('not supported yet');
        end
        disp(['algo = ', algo, ' done']);
        email_notify(username, password, [username, '@ios.ac.cn'], [algo, ' on ', dataset, ' done']);
    end
    cd (exp_settings.prefix_mdcs);
    plot_result(dataset, candiAlgs, 0);
    email_notify(username, password, [username, '@ios.ac.cn'], ['all algo on ', dataset, ' done'], [dataset, '.tex']);
    cd(root_dir);
    disp(['data = ', dataset, ' done']);
end
rmpath(root_dir);