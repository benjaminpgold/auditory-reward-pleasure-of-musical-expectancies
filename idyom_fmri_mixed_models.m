
load('RESULTS_SUMMARY.mat');

%% Mixed-effects models (starting with mDW-IC)

% Start with the full, "beyond optimal" model (per Diggle et al., 2002), using REML to compare nested random-effects structures:

% Try different random-effects structures:
dwIClme1a=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC^2) + (1 | Subject)','FitMethod','REML');
dwIClme1b=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC^2) + (1 + zmDW_IC^2 | Subject)','FitMethod','REML');
dwIClme1c=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC^2) + (1 + zmDW_IC:zmDW_IC | Subject)','FitMethod','REML');
dwIClme1d=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC^2) + (1 + zmDW_IC | Subject)','FitMethod','REML');
compare(dwIClme1a,dwIClme1b)
compare(dwIClme1c,dwIClme1b)
compare(dwIClme1d,dwIClme1b)
% dwIClme1d wins

% Having settled on the optimal random-effects structure, now optimize the fixed effects with ML estimation and LR tests:
clearvars dwIClme*
dwIClme2a=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC + zmDW_IC:zmDW_IC) + (1 + zmDW_IC | Subject)');
dwIClme2b=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC:zmDW_IC) + (1 + zmDW_IC | Subject)');
dwIClme2c=fitlme(RESULTS_SUMMARY,'AvgRating~ zmDW_IC + (1 + zmDW_IC | Subject)');
compare(dwIClme2b,dwIClme2a) % lqm vs. qm
compare(dwIClme2c,dwIClme2a) % lqm vs. lm
% dwIClme2c wins

lmedwIC2=fitlme(RESULTS_SUMMARY,'AvgRating~ zmDW_IC + (1 + zmDW_IC | Subject)','fitmethod','reml'); % 2c with reml
clearvars dwIClme*

lmedwIC2.coefTest
lmedwIC2.Rsquared
% figure; plot(RESULTS_SUMMARY.zmDW_IC,lmedwIC2.fitted,'x');

% Compare nested models only with ML estimation (not REML):
lmedwIC=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC) + (1 + zmDW_IC | Subject)');
lmedwIC2only=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC:zmDW_IC) + (1 + zmDW_IC | Subject)');
lmedwIC2=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_IC + zmDW_IC:zmDW_IC) + (1 + zmDW_IC | Subject)');
compare(lmedwIC,lmedwIC2)
compare(lmedwIC,lmedwIC2only)


load('sub_info.mat');
[b,bnames]=randomEffects(lmedwIC2);

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([1:2:15,19:2:end])); end % corrs between BMRQ or Gold-MSI and sub intercepts
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([1:2:15,19:2:end]),'type','spearman'); end % corrs between BMRQ or Gold-MSI and sub intercepts
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

figure; histogram(b([1:2:15,19:2:end]));
figure; histogram(sub_info{[1:8,10:end],10});

% Only ~sig: Gold-MSI Musical Training vs. intercepts, Spearman, rho = 0.47, p = 0.025, FDR-corrected p = 0.304
% All other ps >= 0.071

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([2:2:16,20:2:end])); end % corrs between BMRQ or Gold-MSI and sub slopes
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([2:2:16,20:2:end]),'type','spearman'); end % corrs between BMRQ or Gold-MSI and sub slopes
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

% All ps >= 0.232


%% The same, but now with normalized liking:

for s=3:26; idx=RESULTS_SUMMARY.Subject==s; normratings(idx)= zscore(RESULTS_SUMMARY.AvgRating(idx)); end; RESULTS_SUMMARY.normLiking = normratings(:);

% Start with the full, "beyond optimal" model (per Diggle et al., 2002), using REML to compare nested random-effects structures:

% Try different random-effects structures:
dwIClme1a=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC^2)','FitMethod','REML');
dwIClme1b=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC^2) + (1 + zmDW_IC^2 | Subject)','FitMethod','REML');
dwIClme1c=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC^2) + (1 + zmDW_IC:zmDW_IC | Subject)','FitMethod','REML');
dwIClme1d=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC^2) + (1 + zmDW_IC | Subject)','FitMethod','REML');
compare(dwIClme1a,dwIClme1b)
compare(dwIClme1a,dwIClme1c)
compare(dwIClme1c,dwIClme1d)
% dwIClme1d wins (by likelihood)

% Having settled on the optimal random-effects structure, now optimize the fixed effects with ML estimation and LR tests:
clearvars dwIClme*
dwIClme2a=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC + zmDW_IC:zmDW_IC) + (1 + zmDW_IC | Subject)');
dwIClme2b=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC:zmDW_IC) + (1 + zmDW_IC | Subject)');
dwIClme2c=fitlme(RESULTS_SUMMARY,'normLiking~ zmDW_IC + (1 + zmDW_IC | Subject)');
compare(dwIClme2b,dwIClme2a) % lqm vs. qm
compare(dwIClme2c,dwIClme2a) % lqm vs. lm
% dwIClme2c wins

lmedwIC2=fitlme(RESULTS_SUMMARY,'normLiking~ zmDW_IC + (1 + zmDW_IC | Subject)','fitmethod','reml'); % 2c with reml
clearvars dwIClme*

lmedwIC2.coefTest
lmedwIC2.Rsquared
% figure; plot(RESULTS_SUMMARY.zmDW_IC,lmedwIC2.fitted,'x');

% Compare nested models only with ML estimation (not REML):
lmedwIC=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC) + (1 + zmDW_IC | Subject)');
lmedwIC2only=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC:zmDW_IC) + (1 + zmDW_IC | Subject)');
lmedwIC2=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_IC + zmDW_IC:zmDW_IC) + (1 + zmDW_IC | Subject)');
compare(lmedwIC,lmedwIC2)
compare(lmedwIC,lmedwIC2only)


%% mDW_Entropy

% Start with the full, "beyond optimal" model (per Diggle et al., 2002), using REML to compare nested random-effects structures:

% Try different random-effects structures:
dwEntlme1a=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 | Subject)','FitMethod','REML');
dwEntlme1b=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 + zmDW_Entropy^2 | Subject)','FitMethod','REML');
dwEntlme1c=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 + zmDW_Entropy:zmDW_Entropy | Subject)','FitMethod','REML');
dwEntlme1d=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 + zmDW_Entropy | Subject)','FitMethod','REML');
compare(dwEntlme1a,dwEntlme1b)
compare(dwEntlme1a,dwEntlme1c)
compare(dwEntlme1a,dwEntlme1d)
% dwEntlme1a wins

% Having settled on the optimal random-effects structure, now optimize the fixed effects with ML estimation and LR tests:
clearvars dwEntlme*
dwEntlme2a=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 | Subject)');
dwEntlme2b=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy:zmDW_Entropy) + (1 | Subject)');
dwEntlme2c=fitlme(RESULTS_SUMMARY,'AvgRating~ zmDW_Entropy + (1 | Subject)');
compare(dwEntlme2b,dwEntlme2a) % lqm vs. qm
compare(dwEntlme2c,dwEntlme2a) % lqm vs. lm
% dwIClme2a wins

lmedwEnt2=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 | Subject)','fitmethod','reml'); % 2a with REML

clearvars dwEntlme*

lmedwEnt2.coefTest
lmedwEnt2.Rsquared
% figure; plot(RESULTS_SUMMARY.zmDW_Entropy,lmedwEnt2.fitted,'x');

% Compare nested models only with ML estimation (not REML):
lmedwEnt=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy) + (1 | Subject)');
lmedwEnt2only=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy:zmDW_Entropy) + (1 | Subject)');
lmedwEnt2=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 | Subject)');
compare(lmedwEnt,lmedwEnt2)
compare(lmedwEnt2only,lmedwEnt2)


load('sub_info.mat');
[b,bnames]=randomEffects(lmedwEnt2);

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([1:8,10:end])); end % corrs between BMRQ or Gold-MSI and sub intercepts
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([1:8,10:end]),'type','spearman'); end % corrs between BMRQ or Gold-MSI and sub intercepts
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

figure; histogram(b([1:8,10:end]));
figure; histogram(sub_info{[1:8,10:end],10});

% Only ~sig: Gold-MSI Musical Training vs. intercepts, Spearman, rho = 0.47, p = 0.025, FDR-corrected p = 0.304
% All other ps >= 0.071


%% The same, but now with normalized liking:

for s=3:26; idx=RESULTS_SUMMARY.Subject==s; normratings(idx)= zscore(RESULTS_SUMMARY.AvgRating(idx)); end; RESULTS_SUMMARY.normLiking = normratings(:);

% Start with the full, "beyond optimal" model (per Diggle et al., 2002), using REML to compare nested random-effects structures:

% Try different random-effects structures:
dwEntlme1a=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_Entropy^2)','FitMethod','REML');
dwEntlme1b=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_Entropy^2) + (-1 + zmDW_Entropy^2 | Subject)','FitMethod','REML');
dwEntlme1c=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_Entropy^2) + (-1 + zmDW_Entropy:zmDW_Entropy | Subject)','FitMethod','REML');
dwEntlme1d=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_Entropy^2) + (-1 + zmDW_Entropy | Subject)','FitMethod','REML');
compare(dwEntlme1a,dwEntlme1b)
compare(dwEntlme1a,dwEntlme1c)
compare(dwEntlme1a,dwEntlme1d)
% dwEntlme1a wins

% Having settled on the optimal random-effects structure, now optimize the fixed effects with ML estimation and LR tests:
clearvars dwEntlme*
dwEntlme2a=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_Entropy^2)');
dwEntlme2b=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_Entropy:zmDW_Entropy)');
dwEntlme2c=fitlme(RESULTS_SUMMARY,'normLiking~ zmDW_Entropy');
compare(dwEntlme2b,dwEntlme2a) % lqm vs. qm
compare(dwEntlme2c,dwEntlme2a) % lqm vs. lm
% dwIClme2a wins

lmedwEnt2=fitlme(RESULTS_SUMMARY,'normLiking~ (zmDW_Entropy^2)','fitmethod','reml'); % 2a with REML

clearvars dwEntlme*

lmedwEnt2.coefTest
lmedwEnt2.Rsquared
% figure; plot(RESULTS_SUMMARY.zmDW_Entropy,lmedwEnt2.fitted,'x');

% Compare nested models only with ML estimation (not REML):
lmedwEnt=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy) + (1 | Subject)');
lmedwEnt2only=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy:zmDW_Entropy) + (1 | Subject)');
lmedwEnt2=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 | Subject)');
compare(lmedwEnt,lmedwEnt2)
compare(lmedwEnt2only,lmedwEnt2)


%% mDW_IC x mDW_Entropy

% Start with the full, "beyond optimal" model (per Diggle et al., 2002), using REML to compare nested random-effects structures:

% Try different random-effects structures:
dwALLlme1a=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 | Subject)','FitMethod','REML');
dwALLlme1b=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy^2)*(mDW_IC^2) | Subject)','FitMethod','REML');
dwALLlme1c=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy^2)*(mDW_IC:mDW_IC) | Subject)','FitMethod','REML');
dwALLlme1d=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy^2)*mDW_IC | Subject)','FitMethod','REML');
dwALLlme1e=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy^2) | Subject)','FitMethod','REML');
dwALLlme1f=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy:mDW_Entropy)*(mDW_IC^2) | Subject)','FitMethod','REML');
dwALLlme1g=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy:mDW_Entropy)*(mDW_IC:mDW_IC) | Subject)','FitMethod','REML');
dwALLlme1h=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy:mDW_Entropy)*mDW_IC | Subject)','FitMethod','REML');
dwALLlme1i=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy:mDW_Entropy) | Subject)','FitMethod','REML');
dwALLlme1j=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_Entropy*(mDW_IC^2) | Subject)','FitMethod','REML');
dwALLlme1k=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_Entropy*(mDW_IC:mDW_IC) | Subject)','FitMethod','REML');
dwALLlme1l=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_Entropy*mDW_IC | Subject)','FitMethod','REML');
dwALLlme1m=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme1n=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC^2)*(mDW_Entropy^2) | Subject)','FitMethod','REML');
dwALLlme1o=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC^2)*(mDW_Entropy:mDW_Entropy) | Subject)','FitMethod','REML');
dwALLlme1p=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC^2)*mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme1q=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC^2) | Subject)','FitMethod','REML');
dwALLlme1r=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC:mDW_IC)*(mDW_Entropy^2) | Subject)','FitMethod','REML');
dwALLlme1s=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC:mDW_IC)*(mDW_Entropy:mDW_Entropy) | Subject)','FitMethod','REML');
dwALLlme1t=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC:mDW_IC)*mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme1u=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC:mDW_IC) | Subject)','FitMethod','REML');
dwALLlme1v=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC*(mDW_Entropy^2) | Subject)','FitMethod','REML');
dwALLlme1w=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC*(mDW_Entropy:mDW_Entropy) | Subject)','FitMethod','REML');
dwALLlme1x=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC*mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme1y=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC | Subject)','FitMethod','REML');
vars = who;
mdl_idx=find(~cellfun(@isempty,regexp(vars,'dwALLlme1')));
for i=1:numel(mdl_idx)
    eval(cell2mat(strcat('mdlAICs(i)=',vars(mdl_idx(i)),'.ModelCriterion.AIC;')));
end
[~,idx]=min(mdlAICs); disp(['winning model is ',cell2mat(vars(mdl_idx(idx)))]); clear vars mdl_idx idx mdlAICs
% dwALLlme1y wins by AIC

dwALLlme2a=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC + mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2b=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC + mDW_Entropy:mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2c=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC + mDW_Entropy^2 | Subject)','FitMethod','REML');
dwALLlme2d=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC:mDW_IC + mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2e=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC:mDW_IC + mDW_Entropy:mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2f=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC:mDW_IC + mDW_Entropy^2 | Subject)','FitMethod','REML');
dwALLlme2g=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC^2 + mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2h=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC^2 + mDW_Entropy:mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2i=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC^2 + mDW_Entropy^2 | Subject)','FitMethod','REML');
vars = who;
mdl_idx=find(~cellfun(@isempty,regexp(vars,'dwALLlme2')));
for i=1:numel(mdl_idx)
    eval(cell2mat(strcat('mdlAICs(i)=',vars(mdl_idx(i)),'.ModelCriterion.AIC;')));
end
[~,idx]=min(mdlAICs); disp(['winning model is ',cell2mat(vars(mdl_idx(idx)))]); clear vars mdl_idx idx mdlAICs
compare(dwALLlme1y,dwALLlme2e)
% dwALLlme1y wins

% Having settled on the optimal random-effects structure, now optimize the fixed effects with ML estimation and LR tests:
clearvars dwALLlme*
dwALLlme3a=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_IC | Subject)');
dwALLlme3b=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC:mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3c=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3d=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy:mDW_Entropy)*(mDW_IC^2) + (1 + mDW_IC | Subject)');
dwALLlme3e=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy:mDW_Entropy)*(mDW_IC:mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3f=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy:mDW_Entropy)*(mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3g=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy)*(mDW_IC^2) + (1 + mDW_IC | Subject)');
dwALLlme3h=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy)*(mDW_IC:mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3i=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy)*(mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3j=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2) + (mDW_IC^2) + (1 + mDW_IC | Subject)');
dwALLlme3k=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2) + (mDW_IC:mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3l=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2) + (mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3m=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy:mDW_Entropy) + (mDW_IC^2) + (1 + mDW_IC | Subject)');
dwALLlme3n=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy:mDW_Entropy) + (mDW_IC:mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3o=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy:mDW_Entropy) + (mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3p=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy) + (mDW_IC^2) + (1 + mDW_IC | Subject)');
dwALLlme3q=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy) + (mDW_IC:mDW_IC) + (1 + mDW_IC | Subject)');
dwALLlme3r=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy) + (mDW_IC) + (1 + mDW_IC | Subject)');
compare(dwALLlme3b,dwALLlme3a)
compare(dwALLlme3c,dwALLlme3a)
compare(dwALLlme3d,dwALLlme3a)
compare(dwALLlme3e,dwALLlme3a)
compare(dwALLlme3f,dwALLlme3a)
compare(dwALLlme3g,dwALLlme3a)
compare(dwALLlme3h,dwALLlme3a)
compare(dwALLlme3i,dwALLlme3a)
compare(dwALLlme3j,dwALLlme3a)
compare(dwALLlme3k,dwALLlme3a)
compare(dwALLlme3l,dwALLlme3a)
compare(dwALLlme3m,dwALLlme3a)
compare(dwALLlme3n,dwALLlme3a)
compare(dwALLlme3o,dwALLlme3a)
compare(dwALLlme3p,dwALLlme3a)
compare(dwALLlme3q,dwALLlme3a)
compare(dwALLlme3r,dwALLlme3a)
% dwALLlme3a wins (by LR/AIC)
lmeALL=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_IC | Subject)','fitmethod','reml');
clearvars dwALLlme*

lmeALL.coefTest
lmeALL.Rsquared


load('sub_info.mat');
[b,bnames]=randomEffects(lmeALL);

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([1:2:15,19:2:end])); end % corrs between BMRQ or Gold-MSI and sub intercepts
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([1:2:15,19:2:end]),'type','spearman'); end % corrs between BMRQ or Gold-MSI and sub intercepts
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

% All ps >= 0.057

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([2:2:16,20:2:end])); end % corrs between BMRQ or Gold-MSI and sub slopes
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

clear r p
for i=1:12; [r(i),p(i)]=corr(sub_info{[1:8,10:end],i}, b([2:2:16,20:2:end]),'type','spearman'); end % corrs between BMRQ or Gold-MSI and sub slopes
disp([sub_info.Properties.VariableNames(1:12); num2cell(p)]);

% All ps >= 0.218

%% The same, but now with normalized liking:

for s=3:26; idx=RESULTS_SUMMARY.Subject==s; normratings(idx)= zscore(RESULTS_SUMMARY.AvgRating(idx)); end; RESULTS_SUMMARY.normLiking = normratings(:);

% Start with the full, "beyond optimal" model (per Diggle et al., 2002), using REML to compare nested random-effects structures:

% Try different random-effects structures:
dwALLlme1a=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 | Subject)','FitMethod','REML');
dwALLlme1b=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy^2)*(mDW_IC^2) | Subject)','FitMethod','REML');
dwALLlme1c=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy^2)*(mDW_IC:mDW_IC) | Subject)','FitMethod','REML');
dwALLlme1d=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy^2)*mDW_IC | Subject)','FitMethod','REML');
dwALLlme1e=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy^2) | Subject)','FitMethod','REML');
dwALLlme1f=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy:mDW_Entropy)*(mDW_IC^2) | Subject)','FitMethod','REML');
dwALLlme1g=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy:mDW_Entropy)*(mDW_IC:mDW_IC) | Subject)','FitMethod','REML');
dwALLlme1h=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy:mDW_Entropy)*mDW_IC | Subject)','FitMethod','REML');
dwALLlme1i=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_Entropy:mDW_Entropy) | Subject)','FitMethod','REML');
dwALLlme1j=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_Entropy*(mDW_IC^2) | Subject)','FitMethod','REML');
dwALLlme1k=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_Entropy*(mDW_IC:mDW_IC) | Subject)','FitMethod','REML');
dwALLlme1l=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_Entropy*mDW_IC | Subject)','FitMethod','REML');
dwALLlme1m=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme1n=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC^2)*(mDW_Entropy^2) | Subject)','FitMethod','REML');
dwALLlme1o=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC^2)*(mDW_Entropy:mDW_Entropy) | Subject)','FitMethod','REML');
dwALLlme1p=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC^2)*mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme1q=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC^2) | Subject)','FitMethod','REML');
dwALLlme1r=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC:mDW_IC)*(mDW_Entropy^2) | Subject)','FitMethod','REML');
dwALLlme1s=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC:mDW_IC)*(mDW_Entropy:mDW_Entropy) | Subject)','FitMethod','REML');
dwALLlme1t=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC:mDW_IC)*mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme1u=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + (mDW_IC:mDW_IC) | Subject)','FitMethod','REML');
dwALLlme1v=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC*(mDW_Entropy^2) | Subject)','FitMethod','REML');
dwALLlme1w=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC*(mDW_Entropy:mDW_Entropy) | Subject)','FitMethod','REML');
dwALLlme1x=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC*mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme1y=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC | Subject)','FitMethod','REML');
vars = who;
mdl_idx=find(~cellfun(@isempty,regexp(vars,'dwALLlme1')));
for i=1:numel(mdl_idx)
    eval(cell2mat(strcat('mdlAICs(i)=',vars(mdl_idx(i)),'.ModelCriterion.AIC;')));
end
[~,idx]=min(mdlAICs); disp(['winning model is ',cell2mat(vars(mdl_idx(idx)))]); clear vars mdl_idx idx mdlAICs
% dwALLlme1u wins by AIC

dwALLlme2a=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC + mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2b=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC + mDW_Entropy:mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2c=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC + mDW_Entropy^2 | Subject)','FitMethod','REML');
dwALLlme2d=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC:mDW_IC + mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2e=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC:mDW_IC + mDW_Entropy:mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2f=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC:mDW_IC + mDW_Entropy^2 | Subject)','FitMethod','REML');
dwALLlme2g=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC^2 + mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2h=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC^2 + mDW_Entropy:mDW_Entropy | Subject)','FitMethod','REML');
dwALLlme2i=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_IC^2)*(mDW_Entropy^2) + (1 + mDW_IC^2 + mDW_Entropy^2 | Subject)','FitMethod','REML');
vars = who;
mdl_idx=find(~cellfun(@isempty,regexp(vars,'dwALLlme2')));
for i=1:numel(mdl_idx)
    eval(cell2mat(strcat('mdlAICs(i)=',vars(mdl_idx(i)),'.ModelCriterion.AIC;')));
end
[~,idx]=min(mdlAICs); disp(['winning model is ',cell2mat(vars(mdl_idx(idx)))]); clear vars mdl_idx idx mdlAICs
compare(dwALLlme1u,dwALLlme2e)
% dwALLlme1u wins

% Having settled on the optimal random-effects structure, now optimize the fixed effects with ML estimation and LR tests:
clearvars dwALLlme*
dwALLlme3a=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3b=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC:mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3c=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3d=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy:mDW_Entropy)*(mDW_IC^2) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3e=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy:mDW_Entropy)*(mDW_IC:mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3f=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy:mDW_Entropy)*(mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3g=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy)*(mDW_IC^2) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3h=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy)*(mDW_IC:mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3i=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy)*(mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3j=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2) + (mDW_IC^2) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3k=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2) + (mDW_IC:mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3l=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2) + (mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3m=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy:mDW_Entropy) + (mDW_IC^2) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3n=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy:mDW_Entropy) + (mDW_IC:mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3o=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy:mDW_Entropy) + (mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3p=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy) + (mDW_IC^2) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3q=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy) + (mDW_IC:mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
dwALLlme3r=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy) + (mDW_IC) + (1 + (mDW_IC:mDW_IC) | Subject)');
compare(dwALLlme3b,dwALLlme3a)
compare(dwALLlme3c,dwALLlme3a)
compare(dwALLlme3d,dwALLlme3a)
compare(dwALLlme3e,dwALLlme3a)
compare(dwALLlme3f,dwALLlme3a)
compare(dwALLlme3g,dwALLlme3a)
compare(dwALLlme3h,dwALLlme3a)
compare(dwALLlme3i,dwALLlme3a)
compare(dwALLlme3j,dwALLlme3a)
compare(dwALLlme3k,dwALLlme3a)
compare(dwALLlme3l,dwALLlme3a)
compare(dwALLlme3m,dwALLlme3a)
compare(dwALLlme3n,dwALLlme3a)
compare(dwALLlme3o,dwALLlme3a)
compare(dwALLlme3p,dwALLlme3a)
compare(dwALLlme3q,dwALLlme3a)
compare(dwALLlme3r,dwALLlme3a)
% dwALLlme3a wins (by LR/AIC)
lmeALL=fitlme(RESULTS_SUMMARY,'normLiking~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + (mDW_IC:mDW_IC) | Subject)','fitmethod','reml');
clearvars dwALLlme*

lmeALL.coefTest
lmeALL.Rsquared


%% Figures

% % % % Main effects:

% Figure 1A
lmedwIC2=fitlme(RESULTS_SUMMARY,'AvgRating~ zmDW_IC + (1 + zmDW_IC | Subject)','fitmethod','reml'); % winning model from above

figure; set(gcf,'color','w'); plot(RESULTS_SUMMARY.zmDW_IC,fitted(lmedwIC2),'.','MarkerSize',9); xlabel('Standardized Mean Duration-Weighted IC'); ylabel('Fitted Mean Liking Rating');
xs=get(gca,'xlim'); xvals=xs(1):.1:xs(2); box off;
yvals=lmedwIC2.Coefficients{1,2}+lmedwIC2.Coefficients{2,2}.*xvals; hold on; plot(xvals,yvals,'linewidth',1.2);
legend('Adjusted Responses','Fitted Model');

% Figure 1B
lmedwEnt2=fitlme(RESULTS_SUMMARY,'AvgRating~ (zmDW_Entropy^2) + (1 | Subject)','fitmethod','reml'); % winning model from above

figure; set(gcf,'color','w'); plot(RESULTS_SUMMARY.zmDW_Entropy,fitted(lmedwEnt2),'.','MarkerSize',9); xlabel('Standardized Mean Duration-Weighted Entropy'); ylabel('Fitted Mean Liking Rating');
xs=get(gca,'xlim'); xvals=xs(1):.1:xs(2); box off;
yvals=lmedwEnt2.Coefficients{1,2}+lmedwEnt2.Coefficients{2,2}.*xvals+lmedwEnt2.Coefficients{3,2}.*xvals.*xvals; hold on; plot(xvals,yvals,'linewidth',1.2);

% Figure 1C
lmeALL=fitlme(RESULTS_SUMMARY,'AvgRating~ (mDW_Entropy^2)*(mDW_IC^2) + (1 + mDW_IC | Subject)','fitmethod','reml'); % winning model from above

[X,Y]=meshgrid(min(RESULTS_SUMMARY.mDW_IC)-.2:.01:max(RESULTS_SUMMARY.mDW_IC)+.2,...
    linspace(min(RESULTS_SUMMARY.mDW_Entropy)-.2,max(RESULTS_SUMMARY.mDW_Entropy)+.2,numel(min(RESULTS_SUMMARY.mDW_IC)-.2:.01:max(RESULTS_SUMMARY.mDW_IC)-.2)));

% Create the fitted response:
Z_ALL=lmeALL.Coefficients{1,2}+lmeALL.Coefficients{2,2}.*X+lmeALL.Coefficients{3,2}.*Y+lmeALL.Coefficients{4,2}.*X.*Y...
    +lmeALL.Coefficients{5,2}.*(X.^2)+lmeALL.Coefficients{6,2}.*(Y.^2)+lmeALL.Coefficients{7,2}.*(X.^2).*Y+lmeALL.Coefficients{8,2}.*X.*(Y.^2)...
    +lmeALL.Coefficients{9,2}.*(X.^2).*(Y.^2);

cmap=[[linspace(0,1,128)', linspace(0,1,128)', ones(128,1)]; [ones(128,1), linspace(1,0,128)', linspace(1,0,128)']];

figure; set(gcf,'color','w'); s=surf(X,Y,Z_ALL,'EdgeColor','interp'); colormap(cmap); h1=colorbar('location','northoutside'); set(h1,'position',[.13,.84,.775,.0508]);
view(90,-90); hold on; contour(X,Y,Z_ALL,7,'k','linewidth',.2); set(gca,'position',[.13,.11,.775,.7])
title(h1,'Fitted Mean Liking','fontsize',11);
xlabel('Raw mDW-IC'); ylabel('Raw mDW-Entropy'); 

% Figure 1D

norm_liking = 0; % 0 to use non-normalized liking ratings; 1 to use normalized liking ratings
if norm_liking == 0; liking_data = RESULTS_SUMMARY.AvgRating; elseif norm_liking == 1; liking_data = RESULTS_SUMMARY.normLiking; end

% Choose only the test stimuli (excluding the practice stimuli):
exp_mDW_IC=idyom_data.zdwmIC(1:50);
exp_mDW_Ent=idyom_data.zdwmEntropy(1:50);

% Establish starting points for k-means clustering (approximated from visually inspecting the 2D distribution) to achieve more balanced clusters:
% figure; plot(exp_mDW_IC, exp_mDW_Ent,'.');
start=[-1.2,-1.5; -1.2,1.5; 0.4,-1; 0.4,1.5; 2,0; 2,1.5];

% Run the k-means clustering:
[master_idx,ctrs]=kmeans([exp_mDW_IC,exp_mDW_Ent],6,'start',start);

% Plot the stimuli and the starting points
figure; set(gcf,'color','w'); plot(exp_mDW_IC,exp_mDW_Ent,'o','MarkerFaceColor','b');
box off;
xlabel('Standardized mDW-IC'); ylabel('Standardized mDW-Entropy');
hold on; plot(start(:,1),start(:,2),'o','markerfacecolor','r');
legend('Stimuli','Clustering starting points','location','southeast');

% Find the stimuli in each cluster, plot them, and save their mDW-IC values in the order of clusters:
figure; set(gcf,'color','w');
h1=plot(start(:,1),start(:,2),'kd','MarkerSize',8,'MarkerFaceColor','k'); hold on;
colors={'b','r','y','g','c','m'}; clear idx
for i=1:6
    plot(exp_mDW_IC(master_idx==i),exp_mDW_Ent(master_idx==i),'o','MarkerFaceColor',colors{i},'MarkerEdgeColor','k'); hold on;
    idx{i}=round(exp_mDW_IC(master_idx==i),3);
end
box off; legend('Cluster start points','Low mDW-IC, Low mDW-Ent','Low mDW-IC, High mDW-Ent','Medium mDW-IC, Low mDW-Ent','Medium mDW-IC, High mDW-Ent','High mDW-IC, Low mDW-Ent','High mDW-IC, High mDW-Ent','location','southeast');
xlabel('Standardized mDW-IC'); ylabel('Standardized mDW-Entropy');

% Collect and save each subject's average liking rating for each stimulus category:
cat_data=nan(numel(unique(RESULTS_SUMMARY.Subject)),6);
for s=1:numel(unique(RESULTS_SUMMARY.Subject))
    sub_temp=nan(52,6); % NaN matrix to place stimulus ratings
    a=1; b=1; c=1; d=1; e=1; f=1; % counters for the 6 stimulus categories
    for i=1:size(RESULTS_SUMMARY,1)
        if RESULTS_SUMMARY.Subject(i)==s+2
            if any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{1})
                sub_temp(a,1)=liking_data(i); a=a+1;
            elseif any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{2})
                sub_temp(b,2)=liking_data(i); b=b+1;
            elseif any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{3})
                sub_temp(c,3)=liking_data(i); c=c+1;
            elseif any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{4})
                sub_temp(d,4)=liking_data(i); d=d+1;
            elseif any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{5})
                sub_temp(e,5)=liking_data(i); e=e+1;
            elseif any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{6})
                sub_temp(f,6)=liking_data(i); f=f+1;
            end
        end
    end
    cat_data(s,:) = nanmean(sub_temp,1); clear sub_temp
end
clear a b c d e f g master_idx idx s i

figure; set(gcf,'color','w'); errorbar(1:3,[nanmean(cat_data(:,1)),nanmean(cat_data(:,3)),nanmean(cat_data(:,5))],...
    [nanstd(cat_data(:,1))/sqrt(numel(find(~isnan(cat_data(:,1))))),nanstd(cat_data(:,3))/numel(find(~isnan(cat_data(:,3)))),nanstd(cat_data(:,5))/sqrt(numel(find(~isnan(cat_data(:,5)))))],...
    'Color',[0.12,0.722,0.114],'LineWidth',1.2); hold on;
errorbar(1:3,[nanmean(cat_data(:,2)),nanmean(cat_data(:,4)),nanmean(cat_data(:,6))],...
    [nanstd(cat_data(:,2))/sqrt(numel(find(~isnan(cat_data(:,2))))),nanstd(cat_data(:,4))/sqrt(numel(find(~isnan(cat_data(:,4))))),nanstd(cat_data(:,6))/sqrt(numel(find(~isnan(cat_data(:,6)))))],...
    'Color',[0.622,0.333,0.753],'LineWidth',1.2);
box off; set(gca,'xtick',[1,2,3]); set(gca,'xticklabel',{'Low','Medium','High'}); xlabel('mDW-IC');
set(gca,'xlim',[.7,3.3]); ylabel('Mean Liking Rating \pm S.E.M.')
legend('Low mDW-Entropy','High mDW-Entropy','location','northeast')

