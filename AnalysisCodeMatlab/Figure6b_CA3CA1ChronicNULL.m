load('~/Documents/MATLAB/allSessions.mat')
%%
sessnum = 80:87;
figureName = 'Figure6b_CA3CA1Chronic_null.pdf';
statName = 'Figure6b_CA3CA1Chronic_null.xls';
propName = 'Figure6b_CA3CA1Chronic_proportions_null.xls';
fStatName = 'Figure6b_CA3CA1Chronic_ftest_null.xls';

%% delta FR
sessionID = [];
cellShank = [];
qdaClass = [];
pvalue = [];
nullvalue = [];
laserResponse = [];
nullResponse = [];
deltaFR = [];
nullFR = [];
%%
for sess = sessnum
    if isfile([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics2.cellinfo.mat'])
        
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics2.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.stimResp2.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.nullResp.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.spikes.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.laser.events.mat'])
    else
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.stimResp2.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.nullResp.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.spikes.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.laser.events.mat'])
    end
    
    pvalue = [pvalue stimResp.pvalue];
    nullvalue = [nullvalue nullResp.pvalue];    
    laserResponse = [laserResponse stimResp.laserResponse];
    nullResponse = [nullResponse nullResp.laserResponse];
    
    trigTimes = laser.timestamps(:,1); %laser on times
    nullTimes = laser.timestamps(:,2); %laser off times
    dur = round(mean(laser.timestamps(:,2)-laser.timestamps(:,1))*1000); %stim time in ms
    
    for cellID = 1:length(spikes.UID)
            spikeTimes = spikes.times{cellID};
            
            [psth,trialspx] = mpsth(spikeTimes,nullTimes,'pre',-dur,'post',3*dur,'chart',0);
            
            index = cellfun(@isempty, trialspx) == 0;
            newTrialSpx = trialspx(index);
            for n = 1:length(newTrialSpx)
                preSpikes{n} = newTrialSpx{n}(newTrialSpx{n}<2*dur-10);
                stimSpikes{n} = newTrialSpx{n}(newTrialSpx{n}>=2*dur+10);
            end
            
            totalDur = (dur/1000)*length(newTrialSpx);
            
            PreFR = sum(cellfun(@length,preSpikes))/totalDur;
            StimFR = sum(cellfun(@length,stimSpikes))/totalDur;
            
            nullFR = [nullFR (StimFR-PreFR)];
            
            clear preSpikes stimSpikes
    end

    for cellID = 1:length(spikes.UID)
        spikeTimes = spikes.times{cellID};
        
        [psth,trialspx] = mpsth(spikeTimes,trigTimes,'pre',dur,'post',dur,'chart',0);
        
        index = cellfun(@isempty, trialspx) == 0;
        newTrialSpx = trialspx(index);
        for n = 1:length(newTrialSpx)
            preSpikes{n} = newTrialSpx{n}(newTrialSpx{n}<-10);
            stimSpikes{n} = newTrialSpx{n}(newTrialSpx{n}>=10);
        end
        
        totalDur = (dur/1000)*length(newTrialSpx);
        
        PreFR = sum(cellfun(@length,preSpikes))/totalDur;
        StimFR = sum(cellfun(@length,stimSpikes))/totalDur;
        
        deltaFR = [deltaFR (StimFR-PreFR)];
        
        clear preSpikes stimSpikes
    end    
    C = cell(1,length(cell_metrics.UID));
    C(:) = {spikes.sessionName};
        
    sessionID = [sessionID C];
    qdaNum = zeros(length(cell_metrics.UID),1)';
    
    for n = 1:length(cell_metrics.UID)
        if cell_metrics.PostClassification(n) == "Classified Interneuron"
            qdaNum(n) = 1;
        elseif cell_metrics.PostClassification(n) == "Classified Pyramidal"
            qdaNum(n) = 0;
        else
            qdaNum(n) = 0.5;
        end
    end
        
    qdaClass = [qdaClass qdaNum];
        
    clear qdaNum
    clear stimResp spikes cell_metrics 
end
%%
pyramidal = [];
interneuron = [];
ratID = {'RatDD' 'RatDF' 'RatDG' 'RatDH'};

intIDX = find(qdaClass>=0.99);
pyrIDX = find(qdaClass<=0.01);

interneuron.deltaFR = deltaFR(intIDX);
pyramidal.deltaFR = deltaFR(pyrIDX);

interneuron.nullFR = nullFR(intIDX);
pyramidal.nullFR = nullFR(pyrIDX);

interneuron.pValue = pvalue(intIDX);
pyramidal.pValue = pvalue(pyrIDX);

interneuron.nullValue = nullvalue(intIDX);
pyramidal.nullValue = nullvalue(pyrIDX);

interneuron.response = laserResponse(intIDX);
pyramidal.response = laserResponse(pyrIDX);

interneuron.nullResponse = nullResponse(intIDX);
pyramidal.nullResponse = nullResponse(pyrIDX);

interneuron.SessID = sessionID(intIDX);
pyramidal.SessID = sessionID(pyrIDX);



% for n = 1:size(ratID,2)
%     ii = find(contains(interneuron.SessID,ratID{n}));
%     pp = find(contains(pyramidal.SessID,ratID{n}));
%     
%     interneuron.ratID(ii) = {ratID{n}};
%     pyramidal.ratID(pp) = {ratID{n}};
% end

%%
clr = [0.3 0.3 0.3; 0.7 0.1 0.1; 0.9 0.3 0.2; 1 0.5 0; 1 0.8 0.2];

PPpval = [];
PPmean = [];
PPn = [];
PPstdev = [];
PPsem = [];

PIpval = [];
PImean = [];
PIn = [];
PIstdev = [];
PIsem = [];

%%
% P = cell(1,length(pyramidal.deltaFR));
% P(:) = {'All'};
% PratID = [pyramidal.ratID P];
% PdeltaFR = [pyramidal.deltaFR pyramidal.deltaFR];
%%
subplot(3,2,1)
h1 = cdfplot([pyramidal.pValue interneuron.pValue]);
h1.LineWidth = 1;
h1.Color = [0 0 0];
hold on
h2 = cdfplot([pyramidal.nullValue interneuron.nullValue]);
h2.LineWidth = 1;
h2.Color = [.5 .5 .5];
h2.LineStyle = '-.';
text(.7,.25,'-- Stim','color','k','fontsize',6,'fontweight','bold','horizontalalignment','left')
text(.7,.1,'-.- Null','color',[.5 .5 .5],'fontsize',6,'fontweight','bold','horizontalalignment','left')
title('p-value CDF','fontsize',8,'FontWeight','normal')
ylabel('CDF')
hold off
box off
ax = gca;
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
%%
subplot(3,2,3)
yline(0,'--')
ylim([-2 7])
hold on
vpout = violinplot([pyramidal.nullFR pyramidal.deltaFR],[ones(length(pyramidal.nullFR),1)' 2*ones(length(pyramidal.deltaFR),1)'],'showdata',false,'showmean',false);
for n = 1:length(vpout)
    vpout(n).EdgeColor = 'none';
    vpout(n).BoxColor = 'none';
    vpout(n).ViolinColor = [0.3 0.3 0.3];
    vpout(n).MedianPlot.MarkerEdgeColor = 'none';
    vpout(n).MedianPlot.MarkerFaceColor = 'none';
    vpout(n).MeanPlot.Color = [1 1 1];
    vpout(n).MeanPlot.XData = [n-.35 n+.35];
    vpout(n).MeanPlot.LineWidth = 2;    
    if rem(n,2) == 0
        scatter(vpout(n).ScatterPlot.XData,vpout(n).ScatterPlot.YData,4,[0,0,0],'filled')   
    else
        scatter(vpout(n).ScatterPlot.XData,vpout(n).ScatterPlot.YData,4,[.5,.5,.5],'filled')
    end
end

vpout = violinplot([pyramidal.nullFR pyramidal.deltaFR],[ones(length(pyramidal.nullFR),1)' 2*ones(length(pyramidal.deltaFR),1)'],'showdata',false,'showmean',false);
for n = 1:length(vpout)
    vpout(n).EdgeColor = 'none';
    vpout(n).BoxColor = 'none';
    vpout(n).ViolinColor = [0.3 0.3 0.3];
    vpout(n).ViolinAlpha = .3;
    vpout(n).MedianPlot.MarkerEdgeColor = 'none';
    vpout(n).MedianPlot.MarkerFaceColor = 'none';
    vpout(n).MeanPlot.Color = [0 0 0];
    vpout(n).MeanPlot.XData = [n-.35 n+.35];
    vpout(n).MeanPlot.LineWidth = 2;    
end

%text(5,max(ylim),num2str(wpyrFR),'fontsize',6,'horizontalalignment','right','verticalalignment','top')
xlim([0 3])
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
ax.XTick = [1 2];
ax.XTickLabel = {'null','stim'};
%xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Pyramidal Cells','fontsize',8,'FontWeight','normal')
%%
subplot(3,2,4)
yline(0,'--')
ylim([-2 7])
hold on
vpout = violinplot([interneuron.nullFR interneuron.deltaFR],[ones(length(interneuron.nullFR),1)' 2*ones(length(interneuron.deltaFR),1)'],'showdata',false,'showmean',false);
for n = 1:length(vpout)
    vpout(n).EdgeColor = 'none';
    vpout(n).BoxColor = 'none';
    vpout(n).ViolinColor = [0.3 0.3 0.3];
    vpout(n).MedianPlot.MarkerEdgeColor = 'none';
    vpout(n).MedianPlot.MarkerFaceColor = 'none';
    vpout(n).MeanPlot.Color = [1 1 1];
    vpout(n).MeanPlot.XData = [n-.35 n+.35];
    vpout(n).MeanPlot.LineWidth = 2;    
    if rem(n,2) == 0
        scatter(vpout(n).ScatterPlot.XData,vpout(n).ScatterPlot.YData,4,[0,0,0],'filled')   
    else
        scatter(vpout(n).ScatterPlot.XData,vpout(n).ScatterPlot.YData,4,[.5,.5,.5],'filled')
    end
end

vpout = violinplot([interneuron.nullFR interneuron.deltaFR],[ones(length(interneuron.nullFR),1)' 2*ones(length(interneuron.deltaFR),1)'],'showdata',false,'showmean',false);
for n = 1:length(vpout)
    vpout(n).EdgeColor = 'none';
    vpout(n).BoxColor = 'none';
    vpout(n).ViolinColor = [0.3 0.3 0.3];
    vpout(n).ViolinAlpha = .3;
    vpout(n).MedianPlot.MarkerEdgeColor = 'none';
    vpout(n).MedianPlot.MarkerFaceColor = 'none';
    vpout(n).MeanPlot.Color = [0 0 0];
    vpout(n).MeanPlot.XData = [n-.35 n+.35];
    vpout(n).MeanPlot.LineWidth = 2;    
end

%text(5,max(ylim),num2str(wpyrFR),'fontsize',6,'horizontalalignment','right','verticalalignment','top')
xlim([0 3])
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
ax.XTick = [1 2];
ax.XTickLabel = {'null','stim'};
%xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Interneurons','fontsize',8,'FontWeight','normal')

%%
ipV = interneuron.nullValue;
iFR = interneuron.nullFR;
    
iInc = length(find(ipV <= 0.01 & iFR >= 0));
iDec = length(find(ipV <= 0.01 & iFR < 0));
iNR = length(find(ipV > 0.01));
    
ppV = pyramidal.nullValue;
pFR = pyramidal.nullFR;
    
pInc = length(find(ppV <= 0.01 & pFR >= 0));
pDec = length(find(ppV <= 0.01 & pFR < 0));
pNR = length(find(ppV > 0.01));
    
intTotal = [iInc iDec iNR];
pyrTotal = [pInc pDec pNR];
    
clear iInc iDec iNR pInc pDec pNR ipV iFR ppV pFR

%%
clr = [0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8];
%% proportion across sessions
subplot(3,2,5)
pTotal = nansum(pyrTotal,1);
b = bar(1:3,[pTotal(1:3)/sum(pTotal(1:3))], ...
    'facecolor', 'flat','edgecolor','none');
b.CData = clr;
box off
ylim([0 1])
set(gca,'ytick',[0 .5 1])
set(gca,'xtick',[])
%set(gca,'xticklabels',{'0','200','400','600','800','1000','1200','1400'})
text(4,.95,'Increase','color',clr(1,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(4,.85,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(4,.75,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

%xlabel('distance from fiber (\mum)')
ax = gca; 
ax.YAxis.FontSize = 7;
ylabel('Proportion','fontsize',7)
title('Pyramidal Cells','fontsize',8,'FontWeight','normal')

subplot(3,2,6)
iTotal = nansum(intTotal,1);
b = bar(1:3,[iTotal(1:3)/sum(iTotal(1:3))], ...
    'facecolor', 'flat','edgecolor','none');
b.CData = clr;
box off
ylim([0 1])
set(gca,'ytick',[0 .5 1])
set(gca,'xtick',[])
%set(gca,'xticklabels',{'0','200','400','600','800','1000','1200','1400'})
text(4,.95,'Increase','color',clr(1,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(4,.85,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(4,.75,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.YAxis.FontSize = 7;
ylabel('Proportion','fontsize',7)
title('Interneurons','fontsize',8,'FontWeight','normal')
%%

ResponseType = ["inc" "dec" "nr"]';
PyramidalProp = pTotal';
InterneuronProp = iTotal';

CA3prop = table(ResponseType,PyramidalProp,InterneuronProp);

writetable(CA3prop,['~/Documents/MATLAB/GenerateFigs/PDFstats/',propName])
%%
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3, 4], 'PaperUnits', 'Inches', 'PaperSize', [4, 5])
print(gcf, '-dpdf',['~/Documents/MATLAB/GenerateFigs/PDFfigs/',figureName])
%%
[h1,pyr,CI,statsp] = vartest2(pyramidal.deltaFR,pyramidal.nullFR);
[h2,int,CI,statsi] = vartest2(interneuron.deltaFR,interneuron.nullFR);

f_test = table(statsp.df1,statsp.fstat,pyr,statsi.df1,statsi.fstat,int);

writetable(f_test,['~/Documents/MATLAB/GenerateFigs/PDFstats/',fStatName]);
%%
close 
clear
