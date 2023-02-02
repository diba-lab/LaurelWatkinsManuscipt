load('~/Documents/MATLAB/allSessions.mat')
%%
sessnum = 80:87;
figureName = 'Figure6b_CA3CA1Chronic.pdf';
statName = 'Figure6b_CA3CA1Chronic.xls';
propName = 'Figure6b_CA3CA1Chronic_proportions.xls';

%% delta FR
sessionID = [];
cellShank = [];
qdaClass = [];
pvalue = [];
laserResponse = [];
deltaFR = [];
%%
for sess = sessnum
    if isfile([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics2.cellinfo.mat'])
        
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics2.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.stimResp2.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.spikes.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.laser.events.mat'])
    else
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.stimResp2.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.spikes.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.laser.events.mat'])
    end
    
    pvalue = [pvalue stimResp.pvalue];
    laserResponse = [laserResponse stimResp.laserResponse];
    
    trigTimes = laser.timestamps(:,1); %laser on times
    dur = round(mean(laser.timestamps(:,2)-laser.timestamps(:,1))*1000); %stim time in ms
    
    for cellID = 1:length(spikes.UID)
        spikeTimes = spikes.times{cellID};
        
        [psth,trialspx] = mpsth(spikeTimes,trigTimes,'pre',dur,'post',dur,'chart',0);
        
        index = cellfun(@isempty, trialspx) == 0;
        newTrialSpx = trialspx(index);
        
        if ~isempty(newTrialSpx)
            
            for n = 1:length(newTrialSpx)
                preSpikes{n} = newTrialSpx{n}(newTrialSpx{n}<-10);
                stimSpikes{n} = newTrialSpx{n}(newTrialSpx{n}>=10);
            end
            
            totalDur = (dur/1000)*length(newTrialSpx);
            
            PreFR = sum(cellfun(@length,preSpikes))/totalDur;
            StimFR = sum(cellfun(@length,stimSpikes))/totalDur;
            
            deltaFR = [deltaFR (StimFR-PreFR)];
        else
            deltaFR = [deltaFR nan];
        end
        
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
    clear stimResp spikes cell_metrics shank
end
%%
pyramidal = [];
interneuron = [];
ratID = {'NintendoDS' 'RatDX' 'RatDZ' 'RatEA'};

intIDX = find(qdaClass>=0.99);
pyrIDX = find(qdaClass<=0.01);

interneuron.deltaFR = deltaFR(intIDX);
pyramidal.deltaFR = deltaFR(pyrIDX);

interneuron.pValue = pvalue(intIDX);
pyramidal.pValue = pvalue(pyrIDX);

interneuron.response = laserResponse(intIDX);
pyramidal.response = laserResponse(pyrIDX);

interneuron.SessID = sessionID(intIDX);
pyramidal.SessID = sessionID(pyrIDX);

for n = 1:size(ratID,2)
    ii = find(contains(interneuron.SessID,ratID{n}));
    pp = find(contains(pyramidal.SessID,ratID{n}));
    
    interneuron.ratID(ii) = {ratID{n}};
    pyramidal.ratID(pp) = {ratID{n}};
end

%%
clr = [0.3 0.3 0.3; 0.3 0.9 0.9; 0.3 0.7 0.7; .3 0.5 0.9; 0.2 .3 .7];

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
P = cell(1,length(pyramidal.deltaFR));
P(:) = {'All'};
PratID = [pyramidal.ratID P];
PdeltaFR = [pyramidal.deltaFR pyramidal.deltaFR];
Ppvalue = [pyramidal.pValue pyramidal.pValue];
%%
vPlotratID = {'All' 'NintendoDS' 'RatDX' 'RatDZ' 'RatEA'};

subplot(2,2,1)
yline(0,'--')
ylim([-2 7])
hold on
vpout = violinplot(PdeltaFR,PratID,'showdata',false,'showmean',true);
for n = 1:length(vpout)
    vpout(n).EdgeColor = 'none';
    vpout(n).BoxColor = 'none';
    vpout(n).ViolinColor = clr(n,:);
    vpout(n).ViolinAlpha = 0.7;
    vpout(n).MedianPlot.MarkerEdgeColor = 'none';
    vpout(n).MedianPlot.MarkerFaceColor = 'none';
    vpout(n).MeanPlot.Color = [0 0 0];
    vpout(n).MeanPlot.XData = [n-.35 n+.35];
    vpout(n).MeanPlot.LineWidth = 2;
    pvals = Ppvalue(contains(PratID,vPlotratID(n)));
    pSess = PratID(contains(PratID,vPlotratID(n)));
    for nn = 1:length(vPlotratID)
        sessPvals = find(contains(pSess,vPlotratID(nn)));
        sigP = sessPvals(pvals(sessPvals) <=0.01);
        nsP  = sessPvals(pvals(sessPvals)>0.01);  
        if ~isempty(sessPvals)
        scatter(vpout(n).ScatterPlot.XData(sigP),vpout(n).ScatterPlot.YData(sigP),6,clr(nn,:),'filled')
        scatter(vpout(n).ScatterPlot.XData(nsP),vpout(n).ScatterPlot.YData(nsP),6,clr(nn,:))
        end
    end
    Prngns = diff(ylim)/20;
    Prngss = diff(ylim)/20+diff(ylim)/40;
    text(n,min(ylim)+Prngns/4,num2str(length(vpout(n).ScatterPlot.XData)),'fontsize',5,'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom')
    [h,p] = ttest(vpout(n).ScatterPlot.YData);
    PPpval(n) = p;
    PPn(n) = length(vpout(n).ScatterPlot.XData);
    PPmean(n) = vpout(n).MeanPlot.YData(1);
    PPstdev(n) = std(vpout(n).ScatterPlot.YData);
    PPsem(n) = PPstdev(n)/sqrt(length(vpout(n).ScatterPlot.YData));

    Prngns = diff(ylim)/10;
    Prngss = (diff(ylim)/10)+.1;
    if p > 0.05
        text(n,max(ylim)-Prngns,'ns','fontsize',6,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
    elseif p <= 0.05 && p > 0.01
        text(n,max(ylim)-Prngss,'*','fontsize',6,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
    elseif p <= 0.01 && p > 0.001
        text(n,max(ylim)-Prngss,'**','fontsize',6,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
    elseif p <= 0.001
        text(n,max(ylim)-Prngss,'***','fontsize',6,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
    end
end

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
ax.XTickLabel = {'All','RatDS','RatDX','RatDZ','RatEA'};
%xticklabels({'All','RatDS','RatDX','RatDZ','RatEA'},'fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Pyramidal Cells')
%%
I = cell(1,length(interneuron.deltaFR));
I(:) = {'All'};
IratID = [interneuron.ratID I];
IdeltaFR = [interneuron.deltaFR interneuron.deltaFR];
Ipvalue = [interneuron.pValue interneuron.pValue];
%%
subplot(2,2,2)
yline(0,'--')
ylim([-2 7])
hold on
viout = violinplot(IdeltaFR,IratID,'showdata',false,'showmean',true);
for n = 1:length(viout)
    viout(n).EdgeColor = 'none';
    viout(n).BoxColor = 'none';
    viout(n).ViolinColor = clr(n,:);
    viout(n).ViolinAlpha = 0.7;
    viout(n).MedianPlot.MarkerEdgeColor = 'none';
    viout(n).MedianPlot.MarkerFaceColor = 'none';
    viout(n).MeanPlot.Color = [0 0 0];
    viout(n).MeanPlot.XData = [n-.35 n+.35];
    viout(n).MeanPlot.LineWidth = 2;
    pvals = Ipvalue(contains(IratID,vPlotratID(n)));
    pSess = IratID(contains(IratID,vPlotratID(n)));
    for nn = 1:length(vPlotratID)
        sessPvals = find(contains(pSess,vPlotratID(nn)));
        sigP = sessPvals(pvals(sessPvals) <=0.01);
        nsP  = sessPvals(pvals(sessPvals)>0.01);  
        if ~isempty(sessPvals)
        scatter(viout(n).ScatterPlot.XData(nsP),viout(n).ScatterPlot.YData(nsP),6,[.5 .5 .5])
        scatter(viout(n).ScatterPlot.XData(sigP),viout(n).ScatterPlot.YData(sigP),6,clr(nn,:),'filled')
        end
    end
    Irngns = diff(ylim)/20;
    Irngss = diff(ylim)/20+diff(ylim)/40;
    text(n,min(ylim)+Irngns/4,num2str(length(viout(n).ScatterPlot.XData)),'fontsize',5,'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom')
    [h,p] = ttest(viout(n).ScatterPlot.YData);
    PIpval(n) = p;
    PIn(n) = length(viout(n).ScatterPlot.XData);
    PImean(n) = viout(n).MeanPlot.YData(1);
    PIstdev(n) = std(viout(n).ScatterPlot.YData);
    PIsem(n) = PIstdev(n)/sqrt(length(viout(n).ScatterPlot.YData));
    Irngns = diff(ylim)/10;
    Irngss = (diff(ylim)/10)+.1;
    if p > 0.05
        text(n,max(ylim)-Irngns,'ns','fontsize',6,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
    elseif p <= 0.05 && p > 0.01
        text(n,max(ylim)-Irngss,'*','fontsize',6,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
    elseif p <= 0.01 && p > 0.001
        text(n,max(ylim)-Irngss,'**','fontsize',6,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
    elseif p <= 0.001
        text(n,max(ylim)-Irngss,'***','fontsize',6,'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom')
    end    
end
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
ax.XTickLabel = {'All','RatDS','RatDX','RatDZ','RatEA'};
ylabel('\Delta FR','fontsize',7)
title('Interneurons')

%% make pvalue table
rID = unique(PratID)';
PyramidalMean = PPmean';
PyramidalSD = PPstdev';
PyramidalSEM = PPsem';
PyramidalN = PPn';
PyramidalPValue = PPpval';

InterneuronMean = PImean';
InterneuronSD = PIstdev';
InterneuronSEM = PIsem';
InterneuronN = PIn';
InterneuronPValue = PIpval';


CA3FRpvals = table(rID,PyramidalMean,PyramidalSD,PyramidalSEM,PyramidalN,PyramidalPValue,...
    InterneuronMean,InterneuronSD,InterneuronSEM,InterneuronN,InterneuronPValue);

writetable(CA3FRpvals,['~/Documents/MATLAB/GenerateFigs/PDFstats/',statName])
%%
ipV = interneuron.pValue;
iFR = interneuron.deltaFR;
    
iInc = length(find(ipV <= 0.01 & iFR >= 0));
iDec = length(find(ipV <= 0.01 & iFR < 0));
iNR = length(find(ipV > 0.01));
    
ppV = pyramidal.pValue;
pFR = pyramidal.deltaFR;
    
pInc = length(find(ppV <= 0.01 & pFR >= 0));
pDec = length(find(ppV <= 0.01 & pFR < 0));
pNR = length(find(ppV > 0.01));
    
intTotal = [iInc iDec iNR];
pyrTotal = [pInc pDec pNR];
    
clear iInc iDec iNR pInc pDec pNR ipV iFR ppV pFR

%%
clr = [0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8];
%% proportion across sessions
subplot(2,2,3)
pTotal = nansum(pyrTotal,1);
b = bar(1:3,[pTotal(1:3)/sum(pTotal(1:3))], ...
    'facecolor', 'flat','edgecolor','none');
b.CData = clr;
box off
ylim([0 1])
set(gca,'ytick',[0 .5 1])
set(gca,'xtick',[])
%set(gca,'xticklabels',{'0','200','400','600','800','1000','1200','1400'})
text(4,.975,'Increase','color',clr(1,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(4,.9,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(4,.825,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

%xlabel('distance from fiber (\mum)')
ax = gca; 
ax.YAxis.FontSize = 7;
ylabel('Proportion','fontsize',7)
%title('Pyramidal Cells','fontsize',8)

subplot(2,2,4)
iTotal = nansum(intTotal,1);
b = bar(1:3,[iTotal(1:3)/sum(iTotal(1:3))], ...
    'facecolor', 'flat','edgecolor','none');
b.CData = clr;
box off
ylim([0 1])
set(gca,'ytick',[0 .5 1])
set(gca,'xtick',[])
%set(gca,'xticklabels',{'0','200','400','600','800','1000','1200','1400'})
text(4,.975,'Increase','color',clr(1,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(4,.9,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(4,.825,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.YAxis.FontSize = 7;
ylabel('Proportion','fontsize',7)
%title('Interneurons','fontsize',8)
%%

ResponseType = ["inc" "dec" "nr"]';
PyramidalProp = pTotal';
InterneuronProp = iTotal';

CA3prop = table(ResponseType,PyramidalProp,InterneuronProp);

writetable(CA3prop,['~/Documents/MATLAB/GenerateFigs/PDFstats/',propName])
%%
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3, 3], 'PaperUnits', 'Inches', 'PaperSize', [3, 3])
print(gcf, '-dpdf',['~/Documents/MATLAB/GenerateFigs/PDFfigs/',figureName]);
%%
close 
clear
