function email_notify(user_name, password, email_to, email_message, attachments)
% send email notification of the progress
%
% usage:
%       send to myself
%       email_notify('my_username', 'my_pass', [], [], {'file1', 'file2'});
%       send to others
%       email_notify('my_username', 'my_pass', 'your_email', [], 'file1');
%
%
email_from = [user_name, '@ios.ac.cn'];
if ~exist('email_to', 'var') || isempty(email_to)
    email_to = email_from;
end
if ~exist('email_message', 'var') || isempty(email_message)
    email_message = ['hi, Matlab 任务已完成.'];
end
if ~exist('attachments', 'var') || isempty(attachments)
    attachments = '';
end

try
    setpref('Internet', 'E_mail', email_from);
    setpref('Internet', 'SMTP_Server', 'mail.ios.ac.cn');
    setpref('Internet', 'SMTP_Username', user_name);
    setpref('Internet', 'SMTP_Password', password);
    java_props = java.lang.System.getProperties;
    java_props.setProperty('mail.smtp.auth', 'true');
    java_props.setProperty('mail.smtp.socketFactory.port.class', 'javax.net.ssl.SSLSocketFactory');
    java_props.setProperty('mail.smtp.socketFactory.port', '465');
    email_subject = 'Matlab 运行报告';
    
    if ischar(attachments)
        attachments = {attachments};
    end
    idx = find(cellfun(@exist, attachments));
    atts = cell(length(idx),1);
    for i = 1:length(idx)
        atts{i} = attachments{idx(i)};
    end
    disp('send email ... ');
    sendmail(email_to, email_subject, email_message, atts);
    disp('send email done ');
catch
    disp('email failed ... ');
    disp(['message = ', email_message]);
end