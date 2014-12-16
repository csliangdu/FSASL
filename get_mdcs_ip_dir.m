function [mdcs_ips, mdcs_dirs] = get_mdcs_ip_dir(n)
mdcs_ips = cell(n,1);
mdcs_dirs = cell(n,1);
parfor i=1:n*1000
    localhost = java.net.InetAddress.getLocalHost();
    ip = localhost.getHostAddress();
    mdcs_ips{i} = char(ip);
    mdcs_dirs{i} = pwd;
end
mdcs_ips = unique(mdcs_ips);
mdcs_dirs = unique(mdcs_dirs);
end