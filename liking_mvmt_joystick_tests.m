
load('RESULTS_SUMMARY.mat');

subs={'I03','I04','I05','I06','I07','I08','I09','I10','I11','I12','I13','I14','I15','I16','I17','I18','I19','I20','I21','I22','I23','I24','I25','I26'};
runs=1:2;

TR=0.654;

for s=1:24
    for r=1:2
        
        bold = load(['./VS_ts/',subs{s},'_run',num2str(runs(r)),'_VS_ts.txt']); % time series of the average BOLD signal in the VS ROI
        
        mvmt = load(['./MvmtParams/',subs{s},'_run',num2str(runs(r)),'_mcf_abs.rms']); % time series of absolute RMS motion, from MCFLIRT
        
        times = load(['./EVs/',subs{s},'_run',num2str(runs(r)),'_norm_Liking.txt']); % onset times of each stimulus
        
        ts = 0:TR:(length(bold)-1)*TR;
        
        sname = str2num(subs{s}(2:3));
        
        ratings = RESULTS_SUMMARY.OnlineRatings(RESULTS_SUMMARY.Subject==sname & RESULTS_SUMMARY.Run==r); % on-line ratings for each stimulus, sampled at 10 Hz
        
        % Normalize ratings by run first (as the preprocessing did to the BOLD signals):
        for stim=1:size(ratings,1); r_no(stim) = length(ratings{stim}); end
        z_ratings = zscore(cell2mat(ratings));
        idx = 1; for stim=1:size(ratings,1); z_rats_stim{stim} = z_ratings(idx : (idx + r_no(stim) - 1)); idx = idx + r_no(stim); end
        
        if size(times,1) ~= size(ratings,1); disp('stim no. mismath'); end
        
        for i=1:size(times,1)
            [~,start_idx] = min(abs(ts - times(i,1)));
            [~,end_idx] = min(abs(ts - (times(i,1) + times(i,2))));
            
            vs_mean(s,r,i) = nanmean(bold(start_idx:end_idx)); % temporal mean of the VS BOLD signal during the stimulus
                        
            mvmt_mean(s,r,i) = nanmean(mvmt(start_idx:end_idx)); % temporal mean of the absolute RMS movement during the stimulus
            mvmt_var(s,r,i) = nanvar(mvmt(start_idx:end_idx)); % temporal variance of the absolute RMS movement during the stimulus
            
            rat_var(s,r,i) = nanvar(z_rats_stim{i}); % temporal variance of the liking ratings/joystick movements during the stimulus
            
            ids{s,r,i} = [num2str(s),'_',num2str(r),'_',num2str(i)]; % just to check the order of data points
            
            bold_mvmt_corr(s,r,i) = corr(mvmt(start_idx:end_idx), bold(start_idx:end_idx)); % correlation between the ROI-averaged BOLD time series and the absolute RMS motion time series for each stimulus
            
        end
    end
end

% tmp=ids(:); % To see the order

subs = repmat(1:24,1,50)'; runs = repmat([ones(25,1); 2*ones(25,1)],24,1); idx=1; for i=1:25; trials(idx:idx+47,1) = i; idx=idx+48; end

trials_tbl = table(subs, runs, trials, vs_mean(:), mvmt_mean(:), mvmt_var(:), rat_var(:), bold_mvmt_corr(:), 'VariableNames',{'sub','run','trial','vs_mean','mvmt_mean','mvmt_var','rat_var','bold_mvmt_corr'});

clearvars -except trials_tbl;


%% VS BOLD mean ~ joystick movement/liking variance

% Compare the likelihood of different random-effects structures with REML (per Diggle et al., 2002):
lme1a = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 + rat_var | sub) + (1 | sub:run)','FitMethod','REML');
lme1b = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 + rat_var | sub:run)','FitMethod','REML');
lme1c = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 | sub) + (1 | sub:run)','FitMethod','REML');
lme1d = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 | sub:run)','FitMethod','REML');
lme1e = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 | sub)','FitMethod','REML');
lme1f = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 + rat_var | sub)','FitMethod','REML');
lme1g = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 | run)','FitMethod','REML');
lme1h = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 + rat_var | run)','FitMethod','REML');
lme1i = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 + rat_var | run) + (1 + rat_var | sub)','FitMethod','REML');
lme1j = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 | run) + (1 + rat_var | sub)','FitMethod','REML');
lme1k = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 + rat_var | run) + (1 | sub)','FitMethod','REML');
lme1l = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 | run) + (1 | sub)','FitMethod','REML');
compare(lme1b, lme1a)
compare(lme1c, lme1b)
compare(lme1d, lme1c)
compare(lme1e, lme1d)
compare(lme1e, lme1f)
compare(lme1e, lme1g)
compare(lme1e, lme1h)
compare(lme1e, lme1i)
compare(lme1e, lme1j)
compare(lme1e, lme1k)
compare(lme1e, lme1l)
% d, e, and g are tied for best but d is most comprehensive

% Having settled on the optimal random-effects structure, now run the model with REML:
clear lme*
lme = fitlme(trials_tbl,'vs_mean ~ rat_var + (1 | sub:run)','fitmethod','reml');
% figure; plot(trials_tbl.rat_var, lme.fitted, '.');
disp(lme);

% Non-significant: p = 0.13907


%% Mean head movement ~ joystick/liking variance

% Compare the likelihood of different random-effects structures with REML (per Diggle et al., 2002):
lme1a = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 + rat_var | sub) + (1 | sub:run)','FitMethod','REML');
lme1b = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 + rat_var | sub:run)','FitMethod','REML');
lme1c = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 | sub) + (1 | sub:run)','FitMethod','REML');
lme1d = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 | sub:run)','FitMethod','REML');
lme1e = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 | sub)','FitMethod','REML');
lme1f = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 + rat_var | sub)','FitMethod','REML');
lme1g = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 | run)','FitMethod','REML');
lme1h = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 + rat_var | run)','FitMethod','REML');
lme1i = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 + rat_var | run) + (1 + rat_var | sub)','FitMethod','REML');
lme1j = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 | run) + (1 + rat_var | sub)','FitMethod','REML');
lme1k = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 + rat_var | run) + (1 | sub)','FitMethod','REML');
lme1l = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 | run) + (1 | sub)','FitMethod','REML');
compare(lme1b, lme1a)
compare(lme1c, lme1a)
compare(lme1d, lme1a)
compare(lme1e, lme1a)
compare(lme1f, lme1a)
compare(lme1g, lme1a)
compare(lme1h, lme1a)
compare(lme1a, lme1i)
compare(lme1j, lme1a)
compare(lme1k, lme1a)
compare(lme1l, lme1a)
% a wins

% Having settled on the optimal random-effects structure, now run the model with REML:
clear lme*
lme = fitlme(trials_tbl,'mvmt_mean ~ rat_var + (1 + rat_var | sub) + (1 | sub:run)','fitmethod','reml');
% figure; plot(trials_tbl.rat_var, lme.fitted, '.');
disp(lme);

% Significant positive effect: beta = 0.08, p = 0.021


%% VS BOLD mean ~ mean head movement

% Compare the likelihood of different random-effects structures with REML (per Diggle et al., 2002):
lme1a = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 + mvmt_mean | sub) + (1 | sub:run)','FitMethod','REML');
lme1b = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 + mvmt_mean | sub:run)','FitMethod','REML');
lme1c = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 | sub) + (1 | sub:run)','FitMethod','REML');
lme1d = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 | sub:run)','FitMethod','REML');
lme1e = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 | sub)','FitMethod','REML');
lme1f = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 + mvmt_mean | sub)','FitMethod','REML');
lme1g = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 | run)','FitMethod','REML');
lme1h = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 + mvmt_mean | run)','FitMethod','REML');
lme1i = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 + mvmt_mean | run) + (1 + mvmt_mean | sub)','FitMethod','REML');
lme1j = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 | run) + (1 + mvmt_mean | sub)','FitMethod','REML');
lme1k = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 + mvmt_mean | run) + (1 | sub)','FitMethod','REML');
lme1l = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 | run) + (1 | sub)','FitMethod','REML');
compare(lme1b, lme1a)
compare(lme1c, lme1b)
compare(lme1d, lme1c)
compare(lme1e, lme1d)
compare(lme1d, lme1f)
compare(lme1d, lme1g)
compare(lme1d, lme1h)
compare(lme1d, lme1i)
compare(lme1d, lme1j)
compare(lme1d, lme1k)
compare(lme1d, lme1l)
% d, e, and g are tied for best but d is most comprehensive

% Having settled on the optimal random-effects structure, now run the model with REML:
clear lme*
lme = fitlme(trials_tbl,'vs_mean ~ mvmt_mean + (1 | sub:run)','fitmethod','reml');
% figure; plot(trials_tbl.mvmt_mean, lme.fitted, '.');
disp(lme);

% Non-significant: p = 0.83124

%% VS BOLD mean ~ head movement variance

% Compare the likelihood of different random-effects structures with REML (per Diggle et al., 2002):
lme1a = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 + mvmt_var | sub) + (1 | sub:run)','FitMethod','REML');
lme1b = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 + mvmt_var | sub:run)','FitMethod','REML');
lme1c = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 | sub) + (1 | sub:run)','FitMethod','REML');
lme1d = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 | sub:run)','FitMethod','REML');
lme1e = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 | sub)','FitMethod','REML');
lme1f = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 + mvmt_var | sub)','FitMethod','REML');
lme1g = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 | run)','FitMethod','REML');
lme1h = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 + mvmt_var | run)','FitMethod','REML');
lme1i = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 + mvmt_var | run) + (1 + mvmt_var | sub)','FitMethod','REML');
lme1j = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 | run) + (1 + mvmt_var | sub)','FitMethod','REML');
lme1k = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 + mvmt_var | run) + (1 | sub)','FitMethod','REML');
lme1l = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 | run) + (1 | sub)','FitMethod','REML');
compare(lme1b, lme1a)
compare(lme1c, lme1b)
compare(lme1d, lme1c)
compare(lme1e, lme1d)
compare(lme1d, lme1f)
compare(lme1g, lme1d)
compare(lme1d, lme1h)
compare(lme1d, lme1i)
compare(lme1d, lme1j)
compare(lme1d, lme1k)
compare(lme1d, lme1l)
% d, e, and g are tied for best but d is most comprehensive

% Having settled on the optimal random-effects structure, now run the model with REML:
clear lme*
lme = fitlme(trials_tbl,'vs_mean ~ mvmt_var + (1 | sub:run)','fitmethod','reml');
% figure; plot(trials_tbl.mvmt_var, lme.fitted, '.');
disp(lme);

% Non-significant: p = 0.093954


%% BOLD-mvmt corr ~ joystick/liking variance

% Compare the likelihood of different random-effects structures with REML (per Diggle et al., 2002):
lme1a = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 + rat_var | sub) + (1 | sub:run)','FitMethod','REML');
lme1b = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 + rat_var | sub:run)','FitMethod','REML');
lme1c = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 | sub) + (1 | sub:run)','FitMethod','REML');
lme1d = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 | sub:run)','FitMethod','REML');
lme1e = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 | sub)','FitMethod','REML');
lme1f = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 + rat_var | sub)','FitMethod','REML');
lme1g = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 | run)','FitMethod','REML');
lme1h = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 + rat_var | run)','FitMethod','REML');
lme1i = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 + rat_var | run) + (1 + rat_var | sub)','FitMethod','REML');
lme1j = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 | run) + (1 + rat_var | sub)','FitMethod','REML');
lme1k = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 + rat_var | run) + (1 | sub)','FitMethod','REML');
lme1l = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 | run) + (1 | sub)','FitMethod','REML');
compare(lme1b, lme1a)
compare(lme1c, lme1b)
compare(lme1d, lme1c)
compare(lme1e, lme1d)
compare(lme1d, lme1f)
compare(lme1d, lme1g)
compare(lme1d, lme1h)
compare(lme1d, lme1i)
compare(lme1d, lme1j)
compare(lme1d, lme1k)
compare(lme1d, lme1l)
% d, e, and g are tied for best but d is most comprehensive

% Having settled on the optimal random-effects structure, now run the model with REML:
clear lme*
lme = fitlme(trials_tbl,'bold_mvmt_corr ~ rat_var + (1 | sub:run)','fitmethod','reml');
% figure; plot(trials_tbl.rat_var, lme.fitted, '.');
disp(lme);

% Non-significant: p = 0.880
