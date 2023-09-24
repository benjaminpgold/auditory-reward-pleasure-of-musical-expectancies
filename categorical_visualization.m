%% Main ROI effects of liking

load('likingROIs.mat');

% Test the main effect of liking in each ROI:
[~,p]=ttest(likingROIs{:,:},0);
disp([likingROIs.Properties.VariableNames;num2cell(p)]);

% Figure 2B
figure; set(gcf,'color','w'); set(gcf,'position',[725,600,560,250]);
h=notBoxPlot(likingROIs{:,:}); hold on; set(gca,'units','pixels');
set(gca,'xticklabel',{'VS','R STG','aPFC'}); ylabel('Participant {\it\beta}s');
plot([0,4],[0,0],'k--'); legend([h(1).data, h(1).mu, h(1).sdPtch],{'Raw data','Mean','Mean \pm S.D.'},'location','northeastoutside');
set(gca,'position',[60,40,370,166.075]); t=title('Liking effects in {\ita priori} ROIs'); set(t,'position',[2,50,0]);
text(2,40,'*','fontsize',20,'horizontalalign','center');

% Look for relationships with the BMRQ or the Goldsmith Musical Sophistication Index:
load('sub_info.mat');

for i=1:size(likingROIs,2)
    for j=1:size(sub_info,2)-1
        tmp_idx = ~isnan(sub_info{:,j}); % find subjects with available data
        [r_sub_likingROI(i,j), p_sub_likingROI(i,j)] = corr(likingROIs{tmp_idx,i}, sub_info{tmp_idx,j}); % correlate their BOLD liking responses vs. questionnaire scores
        [r_sub_likingROI_spear(i,j), p_sub_likingROI_spear(i,j)] = corr(likingROIs{tmp_idx,i}, sub_info{tmp_idx,j},'type','spearman'); % correlate their BOLD liking responses vs. questionnaire scores
    end
end
% No significant effects (all ps ≥ 0.082)

% Look for relationships between the BMRQ or Goldsmiths MSI and average liking:
tmp_idx1 = ~isnan(sub_info{:,:}); % find subjects with available data
tmp_idx2 = find(mean(tmp_idx,2)==1);

for i=1:12; [r_sub_avgliking(i),p_sub_avgliking(i)] = corr(sub_info{tmp_idx2,i}, sub_info{tmp_idx2,13},'type','spearman'); end
% Only ~sig: Gold-MSI Musical Training vs. avg. liking, Spearman, rho = 0.47, p = 0.025, FDR-corrected p = 0.304
% All other ps ≥ 0.071


% Figure 2C
figure(1); set(gcf,'color','w');
figure(2); tmp_lm=fitlm(sub_info{:,13}, likingROIs{:,1}); tmp_grph=plot(tmp_lm); bestline = [tmp_grph(2).XData', tmp_grph(2).YData']; lowlim = [tmp_grph(3).XData', tmp_grph(3).YData']; highlim = [tmp_grph(4).XData', tmp_grph(4).YData']; close figure 2
subplot(131); hold on;
d=plot(sub_info{:,13}, likingROIs{:,1}, 'o', 'markerfacecolor','b', 'markeredgecolor','b'); xlabel('Average liking rating'); ylabel('Liking {\it\beta} (a.u.)'); box off; title({'VS',['F(1,22) = ',num2str(anova(tmp_lm).F(1),'%.2f'),', {\it\beta} = ',num2str(tmp_lm.Coefficients.Estimate(2),'%.2f'),', {\itp} = ',num2str(tmp_lm.coefTest,'%.3f')],''}); set(gca,'ytick',[]);
mainLine=plot(bestline(:,1),bestline(:,2),'color','b','linestyle','--'); patchSaturation=0.15; col=get(mainLine,'color'); edgeColor=col+(1-col)*0.55; patchColor=col+(1-col)*(1-patchSaturation);set(gcf,'renderer','painters'); yP=[lowlim(:,2)',fliplr(highlim(:,2)')]; xP=[lowlim(:,1)',fliplr(lowlim(:,1)')];
p=patch(xP,yP,1,'facecolor',patchColor,'edgecolor','none','facealpha',1); edge(1)=plot(lowlim(:,1),lowlim(:,2),'-','color',edgeColor); edge(2)=plot(highlim(:,1),highlim(:,2),'-','color',edgeColor);
uistack(p,'bottom'); uistack(d,'top');

figure(2); tmp_lm=fitlm(sub_info{:,13}, likingROIs{:,2}); tmp_grph=plot(tmp_lm); bestline = [tmp_grph(2).XData', tmp_grph(2).YData']; lowlim = [tmp_grph(3).XData', tmp_grph(3).YData']; highlim = [tmp_grph(4).XData', tmp_grph(4).YData']; close figure 2
subplot(132); hold on;
d=plot(sub_info{:,13}, likingROIs{:,2}, 'o', 'markerfacecolor','r', 'markeredgecolor','r'); xlabel('Average liking rating'); ylabel('Liking {\it\beta} (a.u.)'); box off; title({'R STG',['F(1,22) = ',num2str(anova(tmp_lm).F(1),'%.2f'),', {\it\beta} = ',num2str(tmp_lm.Coefficients.Estimate(2),'%.2f'),', {\itp} = ',num2str(tmp_lm.coefTest,'%.3f')],''}); set(gca,'ytick',[]);
mainLine=plot(bestline(:,1),bestline(:,2),'color','r','linestyle','--'); patchSaturation=0.15; col=get(mainLine,'color'); edgeColor=col+(1-col)*0.55; patchColor=col+(1-col)*(1-patchSaturation);set(gcf,'renderer','painters'); yP=[lowlim(:,2)',fliplr(highlim(:,2)')]; xP=[lowlim(:,1)',fliplr(lowlim(:,1)')];
p=patch(xP,yP,1,'facecolor',patchColor,'edgecolor','none','facealpha',1); edge(1)=plot(lowlim(:,1),lowlim(:,2),'-','color',edgeColor); edge(2)=plot(highlim(:,1),highlim(:,2),'-','color',edgeColor);
uistack(p,'bottom'); uistack(d,'top');

figure(2); tmp_lm=fitlm(sub_info{:,13}, likingROIs{:,3}); tmp_grph=plot(tmp_lm); bestline = [tmp_grph(2).XData', tmp_grph(2).YData']; lowlim = [tmp_grph(3).XData', tmp_grph(3).YData']; highlim = [tmp_grph(4).XData', tmp_grph(4).YData']; close figure 2
subplot(133); hold on;
d=plot(sub_info{:,13}, likingROIs{:,3}, 'o', 'markerfacecolor',[.2,.6,.1], 'markeredgecolor',[.2,.6,.1]); xlabel('Average liking rating'); ylabel('Liking {\it\beta} (a.u.)'); box off; title({'aPFC',['F(1,22) = ',num2str(anova(tmp_lm).F(1),'%.2f'),', {\it\beta} = ',num2str(tmp_lm.Coefficients.Estimate(2),'%.2f'),', {\itp} = ',num2str(tmp_lm.coefTest,'%.3f')],''}); set(gca,'ytick',[]);
mainLine=plot(bestline(:,1),bestline(:,2),'color',[.2,.6,.1],'linestyle','--'); patchSaturation=0.15; col=get(mainLine,'color'); edgeColor=col+(1-col)*0.55; patchColor=col+(1-col)*(1-patchSaturation);set(gcf,'renderer','painters'); yP=[lowlim(:,2)',fliplr(highlim(:,2)')]; xP=[lowlim(:,1)',fliplr(lowlim(:,1)')];
p=patch(xP,yP,1,'facecolor',patchColor,'edgecolor','none','facealpha',1); edge(1)=plot(lowlim(:,1),lowlim(:,2),'-','color',edgeColor); edge(2)=plot(highlim(:,1),highlim(:,2),'-','color',edgeColor);
uistack(p,'bottom'); uistack(d,'top');
set(gcf,'position',[680,558,835,388]);


% Outlier exploration for Figure 2C:
figure(1); set(gcf,'color','w');
figure(2); tmp_lm=fitlm(sub_info{:,13}, likingROIs{:,1}); tmp_grph=plot(tmp_lm); bestline = [tmp_grph(2).XData', tmp_grph(2).YData']; lowlim = [tmp_grph(3).XData', tmp_grph(3).YData']; highlim = [tmp_grph(4).XData', tmp_grph(4).YData']; close figure 2
subplot(141); hold on;
d=plot(sub_info{:,13}, likingROIs{:,1}, 'o', 'markerfacecolor','b', 'markeredgecolor','b'); xlabel('Average liking rating'); ylabel('Liking {\it\beta} (a.u.)'); box off;
title({'VS (all points)',['F(1,22) = ',num2str(anova(tmp_lm).F(1),'%.2f'),', {\it\beta} = ',num2str(tmp_lm.Coefficients.Estimate(2),'%.2f'),', {\itp} = ',num2str(tmp_lm.coefTest,'%.3f')],''}); set(gca,'ytick',[]);
mainLine=plot(bestline(:,1),bestline(:,2),'color','b','linestyle','--'); patchSaturation=0.15; col=get(mainLine,'color'); edgeColor=col+(1-col)*0.55; patchColor=col+(1-col)*(1-patchSaturation);set(gcf,'renderer','painters'); yP=[lowlim(:,2)',fliplr(highlim(:,2)')]; xP=[lowlim(:,1)',fliplr(lowlim(:,1)')];
p=patch(xP,yP,1,'facecolor',patchColor,'edgecolor','none','facealpha',1); edge(1)=plot(lowlim(:,1),lowlim(:,2),'-','color',edgeColor); edge(2)=plot(highlim(:,1),highlim(:,2),'-','color',edgeColor);
uistack(p,'bottom'); uistack(d,'top');
plot(sub_info{1, 13}, likingROIs{1,1}, 'o', 'markerfacecolor','r', 'markeredgecolor','r');
plot(sub_info{3, 13}, likingROIs{3,1}, 'o', 'markerfacecolor','g', 'markeredgecolor','g');
xs=get(gca,'xlim'); ys=get(gca,'ylim');

tmp_idx = logical([0;ones(23,1)]); % set an index for all but the first data outlier
figure(2); tmp_lm=fitlm(sub_info{tmp_idx, 13}, likingROIs{tmp_idx,1}); tmp_grph=plot(tmp_lm); bestline = [tmp_grph(2).XData', tmp_grph(2).YData']; lowlim = [tmp_grph(3).XData', tmp_grph(3).YData']; highlim = [tmp_grph(4).XData', tmp_grph(4).YData']; close figure 2
subplot(142); hold on;
d=plot(sub_info{tmp_idx, 13}, likingROIs{tmp_idx,1}, 'o', 'markerfacecolor','b', 'markeredgecolor','b'); xlabel('Average liking rating'); ylabel('Liking {\it\beta} (a.u.)'); box off;
title({'VS (no red point)',['F(1,21) = ',num2str(anova(tmp_lm).F(1),'%.2f'),', {\it\beta} = ',num2str(tmp_lm.Coefficients.Estimate(2),'%.2f'),', {\itp} = ',num2str(tmp_lm.coefTest,'%.3f')],''}); set(gca,'ytick',[]);
mainLine=plot(bestline(:,1),bestline(:,2),'color','b','linestyle','--'); patchSaturation=0.15; col=get(mainLine,'color'); edgeColor=col+(1-col)*0.55; patchColor=col+(1-col)*(1-patchSaturation);set(gcf,'renderer','painters'); yP=[lowlim(:,2)',fliplr(highlim(:,2)')]; xP=[lowlim(:,1)',fliplr(lowlim(:,1)')];
p=patch(xP,yP,1,'facecolor',patchColor,'edgecolor','none','facealpha',1); edge(1)=plot(lowlim(:,1),lowlim(:,2),'-','color',edgeColor); edge(2)=plot(highlim(:,1),highlim(:,2),'-','color',edgeColor);
uistack(p,'bottom'); uistack(d,'top');
set(gca,'xlim',xs); set(gca,'ylim',ys);

tmp_idx(1) = 1; tmp_idx(3)=0; % now set an index for all but the second data outlier
figure(2); tmp_lm=fitlm(sub_info{tmp_idx, 13}, likingROIs{tmp_idx,1}); tmp_grph=plot(tmp_lm); bestline = [tmp_grph(2).XData', tmp_grph(2).YData']; lowlim = [tmp_grph(3).XData', tmp_grph(3).YData']; highlim = [tmp_grph(4).XData', tmp_grph(4).YData']; close figure 2
subplot(143); hold on;
d=plot(sub_info{tmp_idx, 13}, likingROIs{tmp_idx,1}, 'o', 'markerfacecolor','b', 'markeredgecolor','b'); xlabel('Average liking rating'); ylabel('Liking {\it\beta} (a.u.)'); box off;
title({'VS (no green point)',['F(1,21) = ',num2str(anova(tmp_lm).F(1),'%.2f'),', {\it\beta} = ',num2str(tmp_lm.Coefficients.Estimate(2),'%.2f'),', {\itp} = ',num2str(tmp_lm.coefTest,'%.3f')],''}); set(gca,'ytick',[]);
mainLine=plot(bestline(:,1),bestline(:,2),'color','b','linestyle','--'); patchSaturation=0.15; col=get(mainLine,'color'); edgeColor=col+(1-col)*0.55; patchColor=col+(1-col)*(1-patchSaturation);set(gcf,'renderer','painters'); yP=[lowlim(:,2)',fliplr(highlim(:,2)')]; xP=[lowlim(:,1)',fliplr(lowlim(:,1)')];
p=patch(xP,yP,1,'facecolor',patchColor,'edgecolor','none','facealpha',1); edge(1)=plot(lowlim(:,1),lowlim(:,2),'-','color',edgeColor); edge(2)=plot(highlim(:,1),highlim(:,2),'-','color',edgeColor);
uistack(p,'bottom'); uistack(d,'top');
set(gca,'xlim',xs); set(gca,'ylim',ys);

tmp_idx(1) = 0; % and finally, set an index for all but both the first and second outliers
figure(2); tmp_lm=fitlm(sub_info{tmp_idx, 13}, likingROIs{tmp_idx,1}); tmp_grph=plot(tmp_lm); bestline = [tmp_grph(2).XData', tmp_grph(2).YData']; lowlim = [tmp_grph(3).XData', tmp_grph(3).YData']; highlim = [tmp_grph(4).XData', tmp_grph(4).YData']; close figure 2
subplot(144); hold on;
d=plot(sub_info{tmp_idx, 13}, likingROIs{tmp_idx,1}, 'o', 'markerfacecolor','b', 'markeredgecolor','b'); xlabel('Average liking rating'); ylabel('Liking {\it\beta} (a.u.)'); box off;
title({'VS (no red or green points)',['F(1,20) = ',num2str(anova(tmp_lm).F(1),'%.2f'),', {\it\beta} = ',num2str(tmp_lm.Coefficients.Estimate(2),'%.2f'),', {\itp} = ',num2str(tmp_lm.coefTest,'%.3f')],''}); set(gca,'ytick',[]);
mainLine=plot(bestline(:,1),bestline(:,2),'color','b','linestyle','--'); patchSaturation=0.15; col=get(mainLine,'color'); edgeColor=col+(1-col)*0.55; patchColor=col+(1-col)*(1-patchSaturation);set(gcf,'renderer','painters'); yP=[lowlim(:,2)',fliplr(highlim(:,2)')]; xP=[lowlim(:,1)',fliplr(lowlim(:,1)')];
p=patch(xP,yP,1,'facecolor',patchColor,'edgecolor','none','facealpha',1); edge(1)=plot(lowlim(:,1),lowlim(:,2),'-','color',edgeColor); edge(2)=plot(highlim(:,1),highlim(:,2),'-','color',edgeColor);
uistack(p,'bottom'); uistack(d,'top');
set(gca,'xlim',xs); set(gca,'ylim',ys);
set(gcf,'position',[1161,737,1124,420]);



%% mDW-IC x Liking ROI effects:

load('RESULTS_SUMMARY.mat');

% Choose only the test stimuli (excluding the practice stimuli):
exp_mDW_IC=idyom_data.zdwmIC(1:50);
exp_mDW_Ent=idyom_data.zdwmEntropy(1:50);

% Divide IC and Ent by kmeans:

% Establish starting points for k-means clustering (approximated from visually inspecting the 2D distribution) to achieve more balanced clusters:
% figure; plot(exp_mDW_IC, exp_mDW_Ent,'.');
start=[-1.2,-1.5; -1.2,1.5; 0.4,-1; 0.4,1.5];

% Run the k-means clustering:
master_idx=kmeans([exp_mDW_IC,exp_mDW_Ent],4,'start',start);

% Find the stimuli in each cluster and save their mDW-IC values in the order of clusters:
for i=1:4
    idx{i}=round(exp_mDW_IC(master_idx==i),3);
end

% Create a NaN matrix to hold each subject's average liking rating for each stimulus category:
liking_data=nan(numel(unique(RESULTS_SUMMARY.Subject)),4);

% Collect and save each subject's average liking rating for each stimulus category:
for s=1:numel(unique(RESULTS_SUMMARY.Subject))
    sub_temp=nan(52,4); % NaN matrix to place stimulus ratings
    a=1; b=1; c=1; d=1; % counters for the 4 stimulus categories
    for i=1:size(RESULTS_SUMMARY,1)
        if RESULTS_SUMMARY.Subject(i)==s+2
            if any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{1})
                sub_temp(a,1)=RESULTS_SUMMARY.AvgRating(i); a=a+1;
            elseif any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{2})
                sub_temp(b,2)=RESULTS_SUMMARY.AvgRating(i); b=b+1;
            elseif any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{3})
                sub_temp(c,3)=RESULTS_SUMMARY.AvgRating(i); c=c+1;
            elseif any(round(RESULTS_SUMMARY.zmDW_IC(i),3)==idx{4})
                sub_temp(d,4)=RESULTS_SUMMARY.AvgRating(i); d=d+1;
            end
        end
    end
    liking_data(s,:) = nanmean(sub_temp,1); clear sub_temp % average the ratings by stimulus category and save them in liking_data
end

load('4cat_ROI_stimcats.mat'); % load each subject's average BOLD response to each stimulus category in each ROI

% For the IC x Liking interaction, convert 0s to NaNs (0s indicate no data for that category, e.g., no stimuli high in mDW-IC that the subject rated positively):
ICxLiking_prior_VS(ICxLiking_prior_VS==0)=NaN; ICxLiking_prior_vmPFC(ICxLiking_prior_vmPFC==0)=NaN; ICxLiking_prior_RSTG(ICxLiking_prior_RSTG==0)=NaN; ICxLiking_ph_VS(ICxLiking_ph_VS==0)=NaN; ICxLiking_ph_aPFC(ICxLiking_ph_aPFC==0)=NaN; ICxLiking_ph_RSTG(ICxLiking_ph_RSTG==0)=NaN;

% Enter BOLD data into a table:
ICxLiking_ROI_stimcats = table(ICxLiking_prior_VS, ICxLiking_prior_vmPFC, ICxLiking_prior_RSTG, ICxLiking_ph_VS, ICxLiking_ph_aPFC, ICxLiking_ph_RSTG);

% For the mDW-IC x mDW-Entropy interaction, each subject encountered all four categories. So simply enter BOLD data into a table:
ICxEnt_ROI_stimcats = table(ICxEnt_prior_VS, ICxEnt_prior_RSTG, ICxEnt_prior_vmPFC, ICxEnt_ph_R_STR, ICxEnt_ph_L_STR);

% Plot the mean +/- SEM BOLD response to each stimulus category in each ROI:
for j=1:size(ICxLiking_ROI_stimcats,2)

    [a,b]=find(ICxLiking_ROI_stimcats{:,j}==0);
    for i=1:numel(a); ICxLiking_ROI_stimcats{:,j}(a(i),b(i))=NaN; end

    figure; set(gcf,'color','w'); hold on;
    errorbar(1:2, nanmean(ICxLiking_ROI_stimcats{:,j}(:,[1,3])), ...
        nanstd(ICxLiking_ROI_stimcats{:,j}(:,[1,3]))./sqrt([numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,1)))), numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,3))))]),'r');
    errorbar(1:2, nanmean(ICxLiking_ROI_stimcats{:,j}(:,[2,4])), ...
        nanstd(ICxLiking_ROI_stimcats{:,j}(:,[2,4]))./sqrt([numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,2)))), numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,4))))]),'b');
    set(gca,'xlim',[0,3]); box off; legend('Low Liking','High Liking'); set(gca,'xtick',1:2); set(gca,'xticklabel',{'Low IC','High IC'}); ylabel('Mean \beta \pm S.E.M.');
    title(strrep(ICxLiking_ROI_stimcats.Properties.VariableNames{j},'ICxLiking_',''),'interpreter','none');
end

% Plot the axiomatic RPE case:
figure; set(gcf,'color','w'); hold on;
plot(1:2, [-.01,-1],'Color',[0.722,0.114,0.12],'LineWidth',1.2); plot(1:2, [.01,1],'Color',[0.333,0.622,0.753],'LineWidth',1.2);
set(gca,'xlim',[0,3]); box off; legend('Low Liking','High Liking'); set(gca,'xtick',1:2); set(gca,'xticklabel',{'Low IC','High IC'}); ylabel('Mean \beta \pm S.E.M.'); set(gca,'ytick',[]);
title('Idealized RPEs');


% Figure 3 B to F
load('ICxLikingROIs.mat');

[~,p]=ttest(ICxLikingROIs{:,:},0); disp([ICxLikingROIs.Properties.VariableNames;num2cell(p)]);
% A priori R STG p = 0.008, other ps ≥ 0.392

figure; set(gcf,'color','w'); set(gcf,'position',[820,328,740,910]);
subplot(311); h=notBoxPlot(ICxLikingROIs{:,:});hold on;
set(gca,'xticklabel',{'VS','R STG','aPFC'}); ylabel('Participant {\it\beta}s');
plot([0,4],[0,0],'k--'); legend([h(1).data, h(1).mu, h(1).sdPtch],{'Raw data','Mean','Mean \pm S.D.'},'location','northeastoutside');
set(gca,'position',[.25,.75,.5,.1825]); t=title('Surprise X Liking effects in {\ita priori} ROIs'); set(t,'position',[2,60,0]);
text(2,50,'**','fontsize',20,'horizontalalign','center');

subplot(323); hold on;
plot(1:2, [.01,1],'Color',[0.333,0.622,0.753],'LineWidth',1.2); plot(1:2, [-.01,-1],'Color',[0.722,0.114,0.12],'LineWidth',1.2);
set(gca,'xlim',[0,3]); box off; set(gca,'xtick',1:2); set(gca,'xticklabel',{'Low surprise','High surprise'}); ylabel('BOLD response (a.u.)'); set(gca,'ytick',[]); set(gca,'ylim',[-1.2,1.2]);
title('Axiomatic RPEs'); set(gca,'position',[.08,.4,.35,.25]);

subplot(324); hold on; j=4;
[a,b]=find(ICxLiking_ROI_stimcats{:,j}==0);
for i=1:numel(a); ICxLiking_ROI_stimcats{:,j}(a(i),b(i))=NaN; end
errorbar(1:2, nanmean(ICxLiking_ROI_stimcats{:,j}(:,[2,4])), ...
    nanstd(ICxLiking_ROI_stimcats{:,j}(:,[2,4]))./sqrt([numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,2)))), numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,4))))]),'Color',[0.333,0.622,0.753],'linewidth',1.2);
errorbar(1:2, nanmean(ICxLiking_ROI_stimcats{:,j}(:,[1,3])), ...
    nanstd(ICxLiking_ROI_stimcats{:,j}(:,[1,3]))./sqrt([numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,1)))), numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,3))))]),'Color',[0.722,0.114,0.12],'LineWidth',1.2);
set(gca,'xlim',[0,3]); box off; l=legend('Liked stimuli','Disliked stimuli'); set(gca,'xtick',1:2); set(gca,'xticklabel',{'Low surprise','High surprise'}); ylabel({'Mean BOLD response (\beta)','\pm S.E.M. (a.u.)'}); set(gca,'ytick',[]); set(gca,'ylim',[-13,0])
title('VS activation cluster'); set(gca,'position',[.5,.4,.35,.25]); set(l,'position',[.78,.6,.1716,.0663]);

subplot(325); hold on; j=3;
[a,b]=find(ICxLiking_ROI_stimcats{:,j}==0);
for i=1:numel(a); ICxLiking_ROI_stimcats{:,j}(a(i),b(i))=NaN; end
errorbar(1:2, nanmean(ICxLiking_ROI_stimcats{:,j}(:,[2,4])), ...
    nanstd(ICxLiking_ROI_stimcats{:,j}(:,[2,4]))./sqrt([numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,2)))), numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,4))))]),'Color',[0.333,0.622,0.753],'linewidth',1.2);
errorbar(1:2, nanmean(ICxLiking_ROI_stimcats{:,j}(:,[1,3])), ...
    nanstd(ICxLiking_ROI_stimcats{:,j}(:,[1,3]))./sqrt([numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,1)))), numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,3))))]),'Color',[0.722,0.114,0.12],'LineWidth',1.2);
set(gca,'xlim',[0,3]); box off; set(gca,'xtick',1:2); set(gca,'xticklabel',{'Low surprise','High surprise'}); ylabel({'Mean BOLD response (\beta)','\pm S.E.M. (a.u.)'}); set(gca,'ytick',[]);
title('R STG {\ita priori} ROI'); set(gca,'position',[.08,.03,.35,.25]);

subplot(326); hold on; j=6;
[a,b]=find(ICxLiking_ROI_stimcats{:,j}==0);
for i=1:numel(a); ICxLiking_ROI_stimcats{:,j}(a(i),b(i))=NaN; end
errorbar(1:2, nanmean(ICxLiking_ROI_stimcats{:,j}(:,[2,4])), ...
    nanstd(ICxLiking_ROI_stimcats{:,j}(:,[2,4]))./sqrt([numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,2)))), numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,4))))]),'Color',[0.333,0.622,0.753],'linewidth',1.2);
errorbar(1:2, nanmean(ICxLiking_ROI_stimcats{:,j}(:,[1,3])), ...
    nanstd(ICxLiking_ROI_stimcats{:,j}(:,[1,3]))./sqrt([numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,1)))), numel(find(~isnan(ICxLiking_ROI_stimcats{:,j}(:,3))))]),'Color',[0.722,0.114,0.12],'LineWidth',1.2);
set(gca,'xlim',[0,3]); box off; set(gca,'xtick',1:2); set(gca,'xticklabel',{'Low surprise','High surprise'}); ylabel({'Mean BOLD response (\beta)','\pm S.E.M. (a.u.)'}); set(gca,'ytick',[]); set(gca,'ylim',[25,43]);
title('R STG activation cluster'); set(gca,'position',[.5,.03,.35,.25]);

% Look for relationships with the BMRQ or the Goldsmith Musical Sophistication Index:
for i=1:size(ICxLikingROIs,2)
    for j=1:size(sub_info,2)-1
        tmp_idx = ~isnan(sub_info{:,j}); % find subjects with available data
        [r_sub_ICxLikingROI(i,j), p_sub_ICxLikingROI(i,j)] = corr(ICxLikingROIs{tmp_idx,i}, sub_info{tmp_idx,j}); % correlate their BOLD ICxLiking responses vs. questionnaire scores
        [r_sub_ICxLikingROI_spear(i,j), p_sub_ICxLikingROI_spear(i,j)] = corr(ICxLikingROIs{tmp_idx,i}, sub_info{tmp_idx,j}); % correlate their BOLD ICxLiking responses vs. questionnaire scores
    end
end
% No significant effects (all ps ≥ 0.067)


%% mDW-IC x mDW-Entropy ROI effects:

load('ICxEntROIs.mat')

% Figure 4 B to D

figure; set(gcf,'color','w'); set(gcf,'units','pixels'); set(gcf,'position',[820,631,740,607]);
subplot(311); h=notBoxPlot(ICxEntROIs{:,1:3});hold on; set(gca,'units','pixels');
set(gca,'xticklabel',{'VS','R STG','aPFC'}); ylabel('Participant {\it\beta}s');
plot([0,4],[0,0],'k--'); legend([h(1).data, h(1).mu, h(1).sdPtch],{'Raw data','Mean','Mean \pm S.D.'},'location','northeastoutside');
set(gca,'fontsize',9); set(gca,'position',[186,380,370,166.075]); t=title('Surprise X Uncertainty effects in {\ita priori} ROIs'); set(t,'position',[2,115,0]);

[~,p]=ttest(ICxEntROIs{:,:},0); disp([ICxEntROIs.Properties.VariableNames;num2cell(p)]);
% All ps ≥ 0.297

subplot(223); set(gca,'units','pixels'); set(gca,'position',[60.2,61.5,259,227.5]); hold on;
errorbar(1:2, nanmean(liking_data(:,[1,3])), ...
    nanstd(liking_data(:,[1,3]))./sqrt([numel(find(~isnan(liking_data(:,1)))), numel(find(~isnan(liking_data(:,3))))]),'Color',[0.12,0.722,0.114],'linewidth',1.2);
errorbar(1:2, nanmean(liking_data(:,[2,4])), ...
    nanstd(liking_data(:,[2,4]))./sqrt([numel(find(~isnan(liking_data(:,2)))), numel(find(~isnan(liking_data(:,4))))]),'Color',[0.622,0.333,0.753],'LineWidth',1.2);
set(gca,'xlim',[0,3]); box off; set(gca,'xtick',1:2); set(gca,'xticklabel',{'Low surprise','High surprise'}); ylabel('Mean liking \pm S.E.M.');
title('Liking ratings');

subplot(224); set(gca,'units','pixels'); set(gca,'position',[371,61.5,259,227.5]); hold on;
errorbar(1:2, nanmean(ICxEnt_ROI_stimcats{:,4}(:,[1,3])), ...
    nanstd(ICxEnt_ROI_stimcats{:,4}(:,[1,3]))./sqrt([numel(find(~isnan(ICxEnt_ROI_stimcats{:,4}(:,1)))), numel(find(~isnan(ICxEnt_ROI_stimcats{:,4}(:,3))))]),'Color',[0.12,0.722,0.114],'linewidth',1.2);
errorbar(1:2, nanmean(ICxEnt_ROI_stimcats{:,4}(:,[2,4])), ...
    nanstd(ICxEnt_ROI_stimcats{:,4}(:,[2,4]))./sqrt([numel(find(~isnan(ICxEnt_ROI_stimcats{:,4}(:,2)))), numel(find(~isnan(ICxEnt_ROI_stimcats{:,4}(:,4))))]),'Color',[0.622,0.333,0.753],'LineWidth',1.2);
set(gca,'xlim',[0,3]); box off; l=legend('Low uncertainty','High uncertainty'); set(gca,'xtick',1:2); set(gca,'xticklabel',{'Low surprise','High surprise'}); ylabel({'Mean BOLD response ({\beta})','\pm S.E.M. (a.u.)'});
title('VS activation cluster'); set(l,'units','pixels'); set(l,'position',[576.192,245.055,133.667,60.333]);

% Look for relationships with the BMRQ or the Goldsmith Musical Sophistication Index:
for i=1:size(ICxEntROIs,2)
    for j=1:size(sub_info,2)-1
        tmp_idx = ~isnan(sub_info{:,j}); % find subjects with available data
        [r_sub_ICxEntROI(i,j), p_sub_ICxEntROI(i,j)] = corr(ICxEntROIs{tmp_idx,i}, sub_info{tmp_idx,j}); % correlate their BOLD ICxEnt responses vs. questionnaire scores
        [r_sub_ICxEntROI_spear(i,j), p_sub_ICxEntROI_spear(i,j)] = corr(ICxEntROIs{tmp_idx,i}, sub_info{tmp_idx,j},'type','spearman'); % correlate their BOLD ICxEnt responses vs. questionnaire scores
    end
end
% Gold-MSI Musical Training vs. a priori aPFC: Spearman's rho = -0.44, p = 0.037, FDR-corrected p = 0.793
