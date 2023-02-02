%% Optrode plots
load('~/Documents/MATLAB/allSessions.mat')
%%
sessnum = 51:54;
figureName = 'Figure5b_CA3AllFiberGAD.pdf';
statName = 'Figure5b_CA3AllFiberGAD_deltaFR.xls';
propName = 'Figure5b_CA3AllFiberGAD_proportions.xls';

%% delta FR
sessionID = [];
cellShank = [];
qdaClass = [];
intensity = [];
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
        
        for pwr = 1:(length(fields)-1)
        
        lsrPWR = str2double(fields{pwr}(5:end));
        sessionID = [sessionID cell_metrics.sessionName];
        cellShank = [cellShank spikes.shankID];
        qdaClass = [qdaClass cell_metrics.qdaClassification];
        intensity = [intensity lsrPWR*ones(length(spikes.shankID),1)'];
        pvalue = [pvalue stimResp.(fields{pwr}).pvalue];
        laserResponse = [laserResponse stimResp.(fields{pwr}).laserResponse];
        
        trigTimes = laser.(fields{pwr})(:,1); %laser on times
        dur = round(mean(laser.(fields{pwr})(:,2)-laser.(fields{pwr})(:,1))*1000); %stim time in ms
        
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
pyramidal = [];
interneuron = [];

intIDX = find(qdaClass>=0.99);
pyrIDX = find(qdaClass<=0.01);

interneuron.deltaFR = deltaFR(intIDX);
pyramidal.deltaFR = deltaFR(pyrIDX);

interneuron.intensity =intensity(intIDX);
pyramidal.intensity =intensity(pyrIDX);

interneuron.pValue = pvalue(intIDX);
pyramidal.pValue = pvalue(pyrIDX);

interneuron.response = laserResponse(intIDX);
pyramidal.response = laserResponse(pyrIDX);

interneuron.SessID = sessionID(intIDX);
pyramidal.SessID = sessionID(pyrIDX);
%%
%clr = [0.3 0.3 0.3; 0.7 0.1 0.1; 0.9 0.3 0.2; 1 0.5 0; 1 0.8 0.2];
clr = [0.7 0.1 0.1; 0.9 0.3 0.2; 0.9 0.3 0.2; 1 0.5 0];


%%
allIntensity = sort(-1*(unique(intensity)));
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
allSess = unique(pyramidal.SessID);
%%
subplot(1,4,1)
yline(0,'--')
hold on
vpout = violinplot(pyramidal.deltaFR,-1*(pyramidal.intensity),'showdata',false,'showmean',false);
ylim([-6 6])
for n = 1:length(vpout)
    vpout(n).EdgeColor = 'none';
    vpout(n).BoxColor = 'none';
    vpout(n).ViolinColor = [0.3 0.3 0.3];
    vpout(n).MedianPlot.MarkerEdgeColor = 'none';
    vpout(n).MedianPlot.MarkerFaceColor = 'none';
    vpout(n).MeanPlot.Color = [1 1 1];
    vpout(n).MeanPlot.XData = [n-.35 n+.35];
    vpout(n).MeanPlot.LineWidth = 2;    
    pvals = pyramidal.pValue(pyramidal.intensity == -1*allIntensity(n) ...
        & ~isnan(pyramidal.deltaFR));
    pSess = pyramidal.SessID(pyramidal.intensity == -1*allIntensity(n) ...
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

vpout = violinplot(pyramidal.deltaFR,-1*(pyramidal.intensity),'showdata',false,'showmean',true);
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

xticklabels({'.001','.003','.01','.17','.75','3.3','6'})
xlim([0 8])
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('Laser Intensity (mW)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Pyramidal Cells','fontsize',8)
%%
subplot(1,4,2)
yline(0,'--')
hold on
viout = violinplot(interneuron.deltaFR,-1*(interneuron.intensity),'showdata',false,'showmean',false);
ylim([-10 10])
for n = 1:length(viout)
    viout(n).EdgeColor = 'none';
    viout(n).BoxColor = 'none';
    viout(n).ViolinColor = [0.3 0.3 0.3];
    viout(n).MedianPlot.MarkerEdgeColor = 'none';
    viout(n).MedianPlot.MarkerFaceColor = 'none';
    viout(n).MeanPlot.Color = [1 1 1];
    viout(n).MeanPlot.XData = [n-.35 n+.35];
    viout(n).MeanPlot.LineWidth = 2;
    pvals = interneuron.pValue(interneuron.intensity == -1*allIntensity(n) ...
        & ~isnan(interneuron.deltaFR));
    pSess = interneuron.SessID(interneuron.intensity == -1*allIntensity(n) ...
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

viout = violinplot(interneuron.deltaFR,-1*(interneuron.intensity),'showdata',false,'showmean',true);
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

xticklabels({'.001','.003','.01','.17','.75','3.3','6'})
xlim([0 8])
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('Laser Intensity (mW)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Interneurons','fontsize',8)

%% make pvalue table
laserIntensity = [.001 .003 .01 .17 .75 3.3 6]';
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


CA3FRpvals = table(laserIntensity,PyramidalMean,PyramidalSD,PyramidalSEM,PyramidalN,PyramidalPValue,...
    InterneuronMean,InterneuronSD,InterneuronSEM,InterneuronN,InterneuronPValue);

writetable(CA3FRpvals,['~/Documents/MATLAB/GenerateFigs/PDFstats/',statName])
%%
allIntensity = -sort(allIntensity);

pyrTotal = [];
intTotal = [];
%%
for nn = 1:length(allIntensity)
    ipV = interneuron.pValue(interneuron.intensity == allIntensity(nn));
    iFR = interneuron.deltaFR(interneuron.intensity == allIntensity(nn));
    
    iInc = length(find(ipV <= 0.01 & iFR >= 0));
    iDec = length(find(ipV <= 0.01 & iFR < 0));
    iNR = length(find(ipV > 0.01));
    
    ppV = pyramidal.pValue(pyramidal.intensity == allIntensity(nn));
    pFR = pyramidal.deltaFR(pyramidal.intensity == allIntensity(nn));
    
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
    0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8; ...
    0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8; ...
    0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8; ...
    0.2 0.4 0.7; 0.4 0.8 1; 0.8 0.8 0.8];
%% proportion across sessions
subplot(1,4,3)
pTotal = nansum(pyrTotal,1);
b = bar(1:21,[pTotal(1:3)/sum(pTotal(1:3)) pTotal(4:6)/sum(pTotal(4:6)) ...
    pTotal(7:9)/sum(pTotal(7:9)) pTotal(10:12)/sum(pTotal(10:12))  pTotal(13:15)/sum(pTotal(13:15)) ...
     pTotal(16:18)/sum(pTotal(16:18))  pTotal(19:21)/sum(pTotal(19:21))], ...
    'facecolor', 'flat','edgecolor','none');
b.CData = clr;
box off
ylim([0 1])
set(gca,'ytick',[0 .25 .5 .75 1])
set(gca,'xtick',[2 5 8 11 14 17 20])
set(gca,'xticklabels',{'.001','.003','.01','.17','.75','3.3','6'})
text(22,.95,'Increase','color',clr(1,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(22,.9,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(22,.85,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('Laser Intensity (mW)','fontsize',7)
ylabel('Proportion','fontsize',7)
title('Pyramidal Cells','fontsize',8)

subplot(1,4,4)
iTotal = nansum(intTotal,1);
b = bar(1:21,[iTotal(1:3)/sum(iTotal(1:3)) iTotal(4:6)/sum(iTotal(4:6)) ...
    iTotal(7:9)/sum(iTotal(7:9)) iTotal(10:12)/sum(iTotal(10:12)) iTotal(13:15)/sum(iTotal(13:15)) ...
    iTotal(16:18)/sum(iTotal(16:18)) iTotal(19:21)/sum(iTotal(19:21))], ...
    'facecolor', 'flat','edgecolor','none');
b.CData = clr;
box off
ylim([0 1])
set(gca,'ytick',[0 .25 .5 .75 1])
set(gca,'xtick',[2 5 8 11 14 17 20])
set(gca,'xticklabels',{'.001','.003','.01','.17','.75','3.3','6'})
text(22,.95,'Increase','color',clr(1,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(22,.9,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(22,.85,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('Laser Intensity (mW)','fontsize',7)
ylabel('Proportion','fontsize',7)
title('Interneurons','fontsize',8)
%%
laserIntensity = [.001 .001 .001 .003 .003 .003 .01 .01 .01 ...
    .17 .17 .17 .75 .75 .75 3.3 3.3 3.3 6 6 6]';
ResponseType = ["inc" "dec" "nr" "inc" "dec" "nr" "inc" "dec" "nr" "inc" "dec" "nr" ...
    "inc" "dec" "nr" "inc" "dec" "nr" "inc" "dec" "nr"]';
PyramidalProp = pTotal';
InterneuronProp = iTotal';

CA3prop = table(laserIntensity,ResponseType,PyramidalProp,InterneuronProp);

writetable(CA3prop,['~/Documents/MATLAB/GenerateFigs/PDFstats/',propName])
%%
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 7, 2], 'PaperUnits', 'Inches', 'PaperSize', [7, 3])
print(gcf, '-dpdf',['C:/Users/Watkinsdejong/OneDrive/Documents/MATLAB/GenerateFigs/PDFfigs/',figureName]);

%%
close all 
clear 
