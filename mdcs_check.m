function [flag_writeable, flag_uploadable, prefix] = mdcs_check(username, password)
% This function check the availability of work and filedependencies
% directory
% pls use your username ....
%
flag_writeable = 1;
flag_uploadable = 1;
prefix = [];
if ispc
    prefix = pwd;
elseif isunix
    message = ['**************************************', ...
        'This is an automatically generated message. ',...
        'You received this email because you used the Matlab Distributed Computing Server (MDCS) ',...
        'in the Laboratory for Computer Science (LCS) recently, and some of  your directories and files ',...
        'are not automatically deleted by the stupid Matlab Job Scheduler (mjs). As a result, ', ...
        'othes failed to submit jobs on these workers, please help to delete them manually! ', ...
        'Thansk for your cooporation! **************************************'];
    message3 = 'Thanks again! Liang Du, a heavier user, from DMGroup@LCS.';
    
    [mdcs_ips, mdcs_dirs] = get_mdcs_ip_dir(100);
    for i1= 1:length(mdcs_ips)
        disp([' worker ip ', mdcs_ips{i1}]);
    end
    n_dir = zeros(length(mdcs_dirs), 1);
    n_dir2 = zeros(length(mdcs_dirs), 1);
    
    for i1= 1:length(mdcs_dirs)
        disp([' worker pwd ', mdcs_dirs{i1}]);
        unix('ls -l ../ |grep work');
        unix('ls -l ../ |grep filedependencies');
        [~, t1] = unix(['ls -l ../ |grep work | grep ', username, ' |wc -l']); % check the owner of work
        [~, t2] = unix(['ls -l ../ |grep filedependencies | grep ', username, ' |wc -l']); % check the owner of filedependencies
        n_dir(i1) = str2double(t1);
        n_dir2(i1) = str2double(t2);
        
        if str2double(t2) < 1
            message2 = [' Please login to the ip = ', ip, ' and manually delete the directory = ' dir];
            [~, dir2_owner] = unix('ls -l ../ |grep filedependencies |awk -F '' '' ''{print $3}''');
            email_notify(username, password, [dir2_owner, '@ios.ac.cn'], [message, message2, message3]);
        end
    end
    
    if sum(n_dir) < length(n_dir)
        warning('You are not the owner of some work directory ...');
        warning('     Write on this dir will failed ....');
        flag_writeable = 0;
    end
    
    if sum(n_dir2) < length(n_dir2)
        warning('You are not the owner of some filedependencies directory ...');
        warning('     upload dependent file on this dir will failed ....');
        flag_uploadable = 0;
    end
    
    root_dir = ['/home/', username];
    if exist(root_dir, 'dir')
        prefix = root_dir;
    end
    
end