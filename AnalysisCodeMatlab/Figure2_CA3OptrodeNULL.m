%% Optrode plots
load('~/Documents/MATLAB/allSessions.mat')
%%
sessnum = 7:18;
figureName = 'Figure2_CA3Optrode_null.pdf';
statName = 'Figure2_CA3Optrode_deltaFR_null.xls';
propName = 'Figure2_CA3Optrode_proportions_null.xls';
fStatName = 'Figure2_CA3Optrode_ftest_null.xls';

%% delta FR
sessionID = [];
stimShank =[];
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
    if isfile([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics4.cellinfo.mat'])
        
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.cell_metrics4.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.stimResp2.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.nullResp.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '-STIMshank.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.spikes.cellinfo.mat'])
        load([sessions{sess}.FileBase, sessions{sess}.FileName, '.laser.events.mat'])
        
        sessionID = [sessionID cell_metrics.sessionName];
        stimShank = [stimShank shank.fiber(1)*ones(length(spikes.UID),1)']; %need to fix this for >1 fiber shank
        cellShank = [cellShank spikes.shankID];
        qdaClass = [qdaClass cell_metrics.qdaClassification];
        pvalue = [pvalue stimResp.pvalue];
        nullvalue = [nullvalue nullResp.pvalue];
        laserResponse = [laserResponse stimResp.laserResponse];
        nullResponse = [nullResponse nullResp.laserResponse];
        
        nullTimes = laser.timestamps(:,2); %laser off times
        trigTimes = laser.timestamps(:,1); %laser off times

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

interneuron.nullFR = nullFR(intIDX);
pyramidal.nullFR = nullFR(pyrIDX);

interneuron.distance =cellDist(intIDX);
pyramidal.distance =cellDist(pyrIDX);

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


%%
for m = 1:length(interneuron.distance)
    if interneuron.distance(m) > 600
        interneuron.distance(m) = 600;
    end
end

for mm = 1:length(pyramidal.distance)
    if pyramidal.distance(mm) > 600
        pyramidal.distance(mm) = 600;
    end
end
%% %% Weighted FR calculation
a1 = 1; a2 = 8; a3 = 16; a4 = 24; 
wpyrFR = mean([a1*pyramidal.deltaFR(pyramidal.distance == 0) a2*pyramidal.deltaFR(pyramidal.distance == 200) ...
    a3*pyramidal.deltaFR(pyramidal.distance == 400) a4*pyramidal.deltaFR(pyramidal.distance == 600)])/49;
wintFR = mean([a1*interneuron.deltaFR(interneuron.distance == 0) a2*interneuron.deltaFR(interneuron.distance == 200) ...
    a3*interneuron.deltaFR(interneuron.distance == 400) a4*interneuron.deltaFR(interneuron.distance == 600)])/49;
%%
allDist = unique(cellDist);

allSess = unique(pyramidal.SessID);

clr = [0.7 0.1 0.1; ...
    0.9 0.3 0.2; 0.9 0.3 0.2;...
    1 0.5 0;...
    0.7 0.3 0.4; 0.7 0.3 0.4; ...
    1 0.8 0.2; 1 0.8 0.2; 1 0.8 0.2; 1 0.8 0.2; ...
    0.9 0.6 0.5;
    0.9 0.7 0.1];

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
ylim([-10 5])
hold on
vpout = violinplot([pyramidal.nullFR pyramidal.deltaFR],[pyramidal.distance-50 pyramidal.distance+50],'showdata',false,'showmean',false);
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

vpout = violinplot([pyramidal.nullFR pyramidal.deltaFR],[pyramidal.distance-50 pyramidal.distance+50],'showdata',false,'showmean',false);
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
xlim([0 9])
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
ax.XTick = [1.5 3.5 5.5 7.5];
ax.XTickLabel = {'0','200','400','600'};
xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Pyramidal Cells','fontsize',8,'FontWeight','normal')

%%
subplot(3,2,4)
yline(0,'--')
ylim([-35 15])
hold on
viout = violinplot([interneuron.nullFR interneuron.deltaFR],[interneuron.distance-50 interneuron.distance+50],'showdata',false,'showmean',false);
for n = 1:length(viout)
    viout(n).EdgeColor = 'none';
    viout(n).BoxColor = 'none';
    viout(n).ViolinColor = [0.3 0.3 0.3];
    viout(n).MedianPlot.MarkerEdgeColor = 'none';
    viout(n).MedianPlot.MarkerFaceColor = 'none';
    viout(n).MeanPlot.Color = [1 1 1];
    viout(n).MeanPlot.XData = [n-.35 n+.35];
    viout(n).MeanPlot.LineWidth = 2;    
    if rem(n,2) == 0
        scatter(viout(n).ScatterPlot.XData,viout(n).ScatterPlot.YData,4,[0,0,0],'filled')   
    else
        scatter(viout(n).ScatterPlot.XData,viout(n).ScatterPlot.YData,4,[.5,.5,.5],'filled')
    end
end

viout = violinplot([interneuron.nullFR interneuron.deltaFR],[interneuron.distance-50 interneuron.distance+50],'showdata',false,'showmean',false);
for n = 1:length(viout)
    viout(n).EdgeColor = 'none';
    viout(n).BoxColor = 'none';
    viout(n).ViolinColor = [0.3 0.3 0.3];
    viout(n).ViolinAlpha = .3;
    viout(n).MedianPlot.MarkerEdgeColor = 'none';
    viout(n).MedianPlot.MarkerFaceColor = 'none';
    viout(n).MeanPlot.Color = [0 0 0];
    viout(n).MeanPlot.XData = [n-.35 n+.35];
    viout(n).MeanPlot.LineWidth = 2;    
end

%text(5,max(ylim),num2str(wpyrFR),'fontsize',6,'horizontalalignment','right','verticalalignment','top')
xlim([0 9])
ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
ax.XTick = [1.5 3.5 5.5 7.5];
ax.XTickLabel = {'0','200','400','600'};
xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Interneurons','fontsize',8,'FontWeight','normal')

%%
DistanceList = 0:200:600;

pyrTotal = [];
intTotal = [];
%%
for nn = 1:length(DistanceList)
    ipV = interneuron.nullValue(interneuron.distance == DistanceList(nn));
    iFR = interneuron.nullFR(interneuron.distance == DistanceList(nn));
    
    iInc = length(find(ipV <= 0.01 & iFR >= 0));
    iDec = length(find(ipV <= 0.01 & iFR < 0));
    iNR = length(find(ipV > 0.01));
    
    ppV = pyramidal.nullValue(pyramidal.distance == DistanceList(nn));
    pFR = pyramidal.nullFR(pyramidal.distance == DistanceList(nn));
    
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
subplot(3,2,5)
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
text(13,.85,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(13,.75,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('Null Proportion','fontsize',7)
title('Pyramidal Cells','fontsize',8,'FontWeight','normal')

subplot(3,2,6)
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
text(13,.85,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(13,.75,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('distance from fiber (\mum)','fontsize',7)
ylabel('Null Proportion','fontsize',7)
title('Interneurons','fontsize',8,'FontWeight','normal')
%%
distance_um = [0 0 0 200 200 200 400 400 400 600 600 600]';
ResponseType = ["inc" "dec" "nr" "inc" "dec" "nr" "inc" "dec" "nr" "inc" "dec" "nr"]';
PyramidalProp = pTotal';
InterneuronProp = iTotal';

CA3prop = table(distance_um,ResponseType,PyramidalProp,InterneuronProp);

writetable(CA3prop,['~/Documents/MATLAB/GenerateFigs/PDFstats/',propName])
%%
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3, 4], 'PaperUnits', 'Inches', 'PaperSize', [4, 5])
print(gcf, '-dpdf',['~/Documents/MATLAB/GenerateFigs/PDFfigs/',figureName]);
%%
distance_um = [0 200 400 600]';

[h1,Pzero,CI,statsp1] = vartest2(pyramidal.deltaFR(pyramidal.distance == 0),pyramidal.nullFR(pyramidal.distance == 0));
[h2,Ptwo,CI,statsp2] = vartest2(pyramidal.deltaFR(pyramidal.distance == 200),pyramidal.nullFR(pyramidal.distance == 200));
[h3,Pfour,CI,statsp3] = vartest2(pyramidal.deltaFR(pyramidal.distance == 400),pyramidal.nullFR(pyramidal.distance == 400));
[h4,Psix,CI,statsp4] = vartest2(pyramidal.deltaFR(pyramidal.distance == 600),pyramidal.nullFR(pyramidal.distance == 600));

PyramidalFTest = [Pzero Ptwo Pfour Psix]';
PyramidalFStat = [statsp1.fstat statsp2.fstat statsp3.fstat statsp4.fstat]';
PyramidalDF = [statsp1.df1 statsp2.df1 statsp3.df1 statsp4.df1]';

[h1,Izero,CI,statsi1] = vartest2(interneuron.deltaFR(interneuron.distance == 0),interneuron.nullFR(interneuron.distance == 0));
[h2,Itwo,CI,statsi2] = vartest2(interneuron.deltaFR(interneuron.distance == 200),interneuron.nullFR(interneuron.distance == 200));
[h3,Ifour,CI,statsi3] = vartest2(interneuron.deltaFR(interneuron.distance == 400),interneuron.nullFR(interneuron.distance == 400));
[h4,Isix,CI,statsi4] = vartest2(interneuron.deltaFR(interneuron.distance == 600),interneuron.nullFR(interneuron.distance == 600));

InterneuronFTest = [Izero Itwo Ifour Isix]';
InterneuronFStat = [statsi1.fstat statsi2.fstat statsi3.fstat statsi4.fstat]';
InterneuronDF = [statsi1.df1 statsi2.df1 statsi3.df1 statsi4.df1]';

f_test = table(distance_um,PyramidalDF,PyramidalFStat,PyramidalFTest,...
    InterneuronDF,InterneuronFStat,InterneuronFTest);

writetable(f_test,['~/Documents/MATLAB/GenerateFigs/PDFstats/',fStatName]);

%%
clear
close
