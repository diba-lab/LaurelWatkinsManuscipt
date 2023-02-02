%% Optrode plots
load('~/Documents/MATLAB/allSessions.mat')
%%
sessnum = 32:39;
figureName = 'Figure4a_CA1LEDstim.pdf';
statName = 'Figure4a_CA1LEDstim_deltaFR.xls';
propName = 'Figure4a_CA1LEDstim_proportions.xls';

%% delta FR
sessionID = [];
cellShank = [];
stimShank = [];
qdaClass = [];
pvalue = [];
laserResponse = []; 
deltaFR = [];
%%
for sess = sessnum
    if isfile([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics4.cellinfo.mat'])
        
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics4.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.stimResp2.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.spikes.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.laser.events.mat'])
        
        fields = fieldnames(laser);
        
        for uLED = 1:(length(fields)-1)
        
        lightLED = str2double(fields{uLED}(5:end));
        if lightLED == 1 || lightLED == 2 || lightLED == 3
            sS = 1;
        elseif lightLED == 4 || lightLED == 5 || lightLED == 6
            sS = 2;
        elseif lightLED == 7 || lightLED == 8 || lightLED == 9
            sS = 3;
        elseif lightLED == 10 || lightLED == 11 || lightLED == 12
            sS = 4;
        end
            
        sessionID = [sessionID cell_metrics.sessionName];
        cellShank = [cellShank spikes.shankID];
        qdaClass = [qdaClass cell_metrics.qdaClassification];
        stimShank = [stimShank sS*ones(length(spikes.shankID),1)'];
        pvalue = [pvalue stimResp.(fields{uLED}).pvalue];
        laserResponse = [laserResponse stimResp.(fields{uLED}).laserResponse];
        
        trigTimes = laser.(fields{uLED})(:,1); %laser on times
        if size(laser.timestamps,2) == 1
            dur = 250;
        else
            dur = round(mean(laser.(fields{uLED})(:,2)-laser.(fields{uLED})(:,1))*1000); %stim time in ms
        end
        
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
            clear preSpikes stimSpikes lsrPWR
        end
        end
        clear stimResp spikes cell_metrics shank
    end
end
%%
cellDist = abs(stimShank-cellShank)*200;

pyramidal = [];
interneuron = [];

intIDX = find(qdaClass>=0.99);
pyrIDX = find(qdaClass<=0.01);

interneuron.deltaFR = deltaFR(intIDX);
pyramidal.deltaFR = deltaFR(pyrIDX);

interneuron.distance =cellDist(intIDX);
pyramidal.distance =cellDist(pyrIDX);

interneuron.pValue = pvalue(intIDX);
pyramidal.pValue = pvalue(pyrIDX);

interneuron.response = laserResponse(intIDX);
pyramidal.response = laserResponse(pyrIDX);

interneuron.SessID = sessionID(intIDX);
pyramidal.SessID = sessionID(pyrIDX);

%%
a1 = 1; a2 = 8; a3 = 16; a4 = 24; 
wpyrFR = nanmean([a1*pyramidal.deltaFR(pyramidal.distance == 0) a2*pyramidal.deltaFR(pyramidal.distance == 200) ...
    a3*pyramidal.deltaFR(pyramidal.distance == 400) a4*pyramidal.deltaFR(pyramidal.distance == 600)])/49;
wintFR = nanmean([a1*interneuron.deltaFR(interneuron.distance == 0) a2*interneuron.deltaFR(interneuron.distance == 200) ...
    a3*interneuron.deltaFR(interneuron.distance == 400) a4*interneuron.deltaFR(interneuron.distance == 600)])/49;
%%
allDist = unique(cellDist);

allSess = unique(pyramidal.SessID);

clr = [0.7 0.1 0.1; 0.7 0.1 0.1; 0.7 0.1 0.1; 0.7 0.1 0.1; ...
    0.9 0.3 0.2; ...
    1 0.5 0;...
    0.7 0.3 0.4; 0.7 0.3 0.4];

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
subplot(1,4,1)
yline(0,'--')
ylim([-15 25])
hold on
vpout = violinplot(pyramidal.deltaFR,pyramidal.distance,'showdata',false,'showmean',false);
for n = 1:length(vpout)
    vpout(n).EdgeColor = 'none';
    vpout(n).BoxColor = 'none';
    vpout(n).ViolinColor = [0.3 0.3 0.3];
    vpout(n).MedianPlot.MarkerEdgeColor = 'none';
    vpout(n).MedianPlot.MarkerFaceColor = 'none';
    vpout(n).MeanPlot.Color = [1 1 1];
    vpout(n).MeanPlot.XData = [n-.35 n+.35];
    vpout(n).MeanPlot.LineWidth = 2;    
    pvals = pyramidal.pValue(pyramidal.distance == allDist(n) ...
        & ~isnan(pyramidal.deltaFR));
    pSess = pyramidal.SessID(pyramidal.distance == allDist(n) ...
        & ~isnan(pyramidal.deltaFR));
    for nn = 1:length(allSess)
        sessPvals = find(contains(pSess,allSess(nn)));
        sigP = sessPvals(pvals(sessPvals) <=0.01);
        nsP  = sessPvals(pvals(sessPvals)>0.01);  
        if ~isempty(sessPvals)
        scatter(vpout(n).ScatterPlot.XData(sigP),vpout(n).ScatterPlot.YData(sigP),2,clr(nn,:),'filled')
        scatter(vpout(n).ScatterPlot.XData(nsP),vpout(n).ScatterPlot.YData(nsP),2,clr(nn,:))
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

    if p > 0.05
        text(n,max(ylim)-Prngns,'ns','fontsize',5,'HorizontalAlignment','center', ...
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

vpout = violinplot(pyramidal.deltaFR,pyramidal.distance,'showdata',false,'showmean',true);
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
xlim([0 5])
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Pyramidal Cells','fontsize',8)
%%
subplot(1,4,2)
yline(0,'--')
hold on
ylim([-15 25])
viout = violinplot(interneuron.deltaFR,interneuron.distance,'showdata',false,'showmean',false);
for n = 1:length(viout)
    viout(n).EdgeColor = 'none';
    viout(n).BoxColor = 'none';
    viout(n).ViolinColor = [0.3 0.3 0.3];
    viout(n).MedianPlot.MarkerEdgeColor = 'none';
    viout(n).MedianPlot.MarkerFaceColor = 'none';
    viout(n).MeanPlot.Color = [1 1 1];
    viout(n).MeanPlot.XData = [n-.35 n+.35];
    viout(n).MeanPlot.LineWidth = 2;
    pvals = interneuron.pValue(interneuron.distance == allDist(n) ...
        & ~isnan(interneuron.deltaFR));
    pSess = interneuron.SessID(interneuron.distance == allDist(n) ...
        & ~isnan(interneuron.deltaFR));
    for nn = 1:length(allSess)
        sessPvals = find(contains(pSess,allSess(nn)));
        sigP = sessPvals(pvals(sessPvals) <=0.01);
        nsP  = sessPvals(pvals(sessPvals)>0.01);  
        if ~isempty(sessPvals)
        scatter(viout(n).ScatterPlot.XData(sigP),viout(n).ScatterPlot.YData(sigP),2,clr(nn,:),'filled')
        scatter(viout(n).ScatterPlot.XData(nsP),viout(n).ScatterPlot.YData(nsP),2,clr(nn,:))
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

    if p > 0.05
        text(n,max(ylim)-Irngns,'ns','fontsize',5,'HorizontalAlignment','center', ...
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

viout = violinplot(interneuron.deltaFR,interneuron.distance,'showdata',false,'showmean',true);
for n = 1:length(viout)
    viout(n).EdgeColor = 'none';
    viout(n).BoxColor = 'none';
    viout(n).ViolinColor = [0.3 0.3 0.3];
    viout(n).ViolinAlpha = 0.3;
    viout(n).MedianPlot.MarkerEdgeColor = 'none';
    viout(n).MedianPlot.MarkerFaceColor = 'none';
    viout(n).MeanPlot.Color = [0 0 0];
    viout(n).MeanPlot.XData = [n-.35 n+.35];
    viout(n).MeanPlot.LineWidth = 2;
end

%text(5,max(ylim),num2str(wintFR),'fontsize',6,'horizontalalignment','right','verticalalignment','top')
xlim([0 5])
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Interneurons','fontsize',8)

%% make pvalue table
distance_um = [0 200 400 600]';
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


CA3FRpvals = table(distance_um,PyramidalMean,PyramidalSD,PyramidalSEM,PyramidalN,PyramidalPValue,...
    InterneuronMean,InterneuronSD,InterneuronSEM,InterneuronN,InterneuronPValue);

writetable(CA3FRpvals,['~/Documents/MATLAB/GenerateFigs/PDFstats/',statName])
%%
DistanceList = 0:200:600;

pyrTotal = [];
intTotal = [];
%%
for nn = 1:length(DistanceList)
    ipV = interneuron.pValue(interneuron.distance == DistanceList(nn));
    iFR = interneuron.deltaFR(interneuron.distance == DistanceList(nn));
    
    iInc = length(find(ipV <= 0.01 & iFR >= 0));
    iDec = length(find(ipV <= 0.01 & iFR < 0));
    iNR = length(find(ipV > 0.01));
    
    ppV = pyramidal.pValue(pyramidal.distance == DistanceList(nn));
    pFR = pyramidal.deltaFR(pyramidal.distance == DistanceList(nn));
    
    pInc = length(find(ppV <= 0.01 & pFR >= 0));
    pDec = length(find(ppV <= 0.01 & pFR < 0));
    pNR = length(find(ppV > 0.01));
    
    intTotal = [intTotal iInc iDec iNR];
    pyrTotal = [pyrTotal pInc pDec pNR];
    
    clear iInc iDec iNR pInc pDec pNR ipV iFR ppV pFR

end
%%
clr = [0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8; ...
    0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8; ...
    0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8; ...
    0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8];
%% proportion across sessions
subplot(1,4,3)
pTotal = nansum(pyrTotal,1);
b = bar(1:12,[pTotal(1:3)/sum(pTotal(1:3)) pTotal(4:6)/sum(pTotal(4:6)) ...
    pTotal(7:9)/sum(pTotal(7:9)) pTotal(10:12)/sum(pTotal(10:12))], ...
    'facecolor', 'flat','edgecolor','none');
b.CData = clr;
box off
ylim([0 1])
set(gca,'ytick',[0 .25 .5 .75 1])
set(gca,'xtick',[2 5 8 11])
set(gca,'xticklabels',{'0','200','400','600'})
text(13,.95,'Increase','color',clr(1,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(13,.9,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(13,.85,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('Proportion','fontsize',7)
title('Pyramidal Cells','fontsize',8)

subplot(1,4,4)
iTotal = nansum(intTotal,1);
b = bar(1:12,[iTotal(1:3)/sum(iTotal(1:3)) iTotal(4:6)/sum(iTotal(4:6)) ...
    iTotal(7:9)/sum(iTotal(7:9)) iTotal(10:12)/sum(iTotal(10:12))], ...
    'facecolor', 'flat','edgecolor','none');
b.CData = clr;
box off
ylim([0 1])
set(gca,'ytick',[0 .25 .5 .75 1])
set(gca,'xtick',[2 5 8 11])
set(gca,'xticklabels',{'0','200','400','600'})
text(13,.95,'Increase','color',clr(1,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(13,.9,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(13,.85,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('Proportion','fontsize',7)
title('Interneurons','fontsize',8)
%%
distance_um = [0 0 0 200 200 200 400 400 400 600 600 600]';
ResponseType = ["inc" "dec" "nr" "inc" "dec" "nr" "inc" "dec" "nr" "inc" "dec" "nr"]';
PyramidalProp = pTotal';
InterneuronProp = iTotal';

CA3prop = table(distance_um,ResponseType,PyramidalProp,InterneuronProp);

writetable(CA3prop,['~/Documents/MATLAB/GenerateFigs/PDFstats/',propName])
%%
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7, 2], 'PaperUnits', 'Inches', 'PaperSize', [7, 3])
print(gcf, '-dpdf',['~/Documents/MATLAB/GenerateFigs/PDFfigs/',figureName]);
%%
close 
clear