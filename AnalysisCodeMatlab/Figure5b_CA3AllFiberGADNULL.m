%% Optrode plots
load('~/Documents/MATLAB/allSessions.mat')
%%
sessnum = 51:54;
figureName = 'Figure5b_CA3AllFiberGAD_null.pdf';
statName = 'Figure5b_CA3AllFiberGAD_deltaFR_null.xls';
propName = 'Figure5b_CA3AllFiberGAD_proportions_null.xls';
fStatName = 'Figure5b_CA3AllFiberGAD_ftest_null.xls';
%% delta FR
sessionID = [];
stimShank =[];
cellShank = [];
qdaClass = [];
intensity = [];
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

        for pwr = 1:(length(fields)-1)
        
        lsrPWR = str2double(fields{pwr}(5:end));
        %sessionID = [sessionID cell_metrics.sessionName];
        %cellShank = [cellShank spikes.shankID];
        %qdaClass = [qdaClass cell_metrics.qdaClassification];
        %intensity = [intensity lsrPWR*ones(length(spikes.shankID),1)'];
        nullvalue = [nullvalue nullResp.(fields{pwr}).pvalue];
        nullResponse = [nullResponse nullResp.(fields{pwr}).laserResponse];
        
        nullTimes = laser.(fields{pwr})(:,2); %laser off times
        dur = round(mean(laser.(fields{pwr})(:,2)-laser.(fields{pwr})(:,1))*1000); %stim time in ms
        
        for cellID = 1:length(spikes.UID)
            spikeTimes = spikes.times{cellID};
            
            [psth,trialspx] = mpsth(spikeTimes,nullTimes,'pre',-dur,'post',3*dur,'chart',0);
            
            index = cellfun(@isempty, trialspx) == 0;
            newTrialSpx = trialspx(index);
            if ~isempty(newTrialSpx)
                for n = 1:length(newTrialSpx)
                    preSpikes{n} = newTrialSpx{n}(newTrialSpx{n}<2*dur-10);
                    stimSpikes{n} = newTrialSpx{n}(newTrialSpx{n}>=2*dur+10);
                end
                
                totalDur = (dur/1000)*length(newTrialSpx);
                
                PreFR = sum(cellfun(@length,preSpikes))/totalDur;
                StimFR = sum(cellfun(@length,stimSpikes))/totalDur;
                
                nullFR = [nullFR (StimFR-PreFR)];
            else
                nullFR = [nullFR nan];
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

interneuron.nullFR = nullFR(intIDX);
pyramidal.nullFR = nullFR(pyrIDX);

interneuron.intensity =intensity(intIDX);
pyramidal.intensity =intensity(pyrIDX);

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
%clr = [0.3 0.3 0.3; 0.7 0.1 0.1; 0.9 0.3 0.2; 1 0.5 0; 1 0.8 0.2];
clr = [0.7 0.1 0.1; 0.7 0.1 0.1; 0.9 0.3 0.2; 0.9 0.3 0.2; 1 0.5 0];


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
ylim([-6 6])
hold on
vpout = violinplot([pyramidal.nullFR pyramidal.deltaFR],[pyramidal.intensity-5 pyramidal.intensity+5],'showdata',false,'showmean',false);
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

vpout = violinplot([pyramidal.nullFR pyramidal.deltaFR],[pyramidal.intensity-5 pyramidal.intensity+5],'showdata',false,'showmean',false);
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
xlim([0 15])
ax = gca; 
ax.XTickLabel ={'.001','.003','.01','.17','.75','3.3','6'};
ax.XTick = [1.5 3.5 5.5 7.5 9.5 11.5 13.5];
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('Laser Intensity (mW)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Pyramidal Cells','fontsize',8,'FontWeight','normal')
%%
subplot(3,2,4)
yline(0,'--')
ylim([-10 10])
hold on
vpout = violinplot([interneuron.nullFR interneuron.deltaFR],[interneuron.intensity-5 interneuron.intensity+5],'showdata',false,'showmean',false);
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

vpout = violinplot([interneuron.nullFR interneuron.deltaFR],[interneuron.intensity-5 interneuron.intensity+5],'showdata',false,'showmean',false);
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
xlim([0 15])
ax = gca; 
ax.XTickLabel ={'.001','.003','.01','.17','.75','3.3','6'};
ax.XTick = [1.5 3.5 5.5 7.5 9.5 11.5 13.5];
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('Laser Intensity (mW)','fontsize',7)
ylabel('\Delta FR','fontsize',7)
title('Interneurons','fontsize',8,'FontWeight','normal')

%% make pvalue table
%laserIntensity = [.001 .003 .01 .17 .75 3.3 6]';
%PyramidalMean = PPmean';
% PyramidalSD = PPstdev';
% PyramidalSEM = PPsem';
% PyramidalN = PPn';
% PyramidalPValue = PPpval';
% 
% InterneuronMean = PImean';
% InterneuronSD = PIstdev';
% InterneuronSEM = PIsem';
% InterneuronN = PIn';
% InterneuronPValue = PIpval';
% 
% 
% CA3FRpvals = table(laserIntensity,PyramidalMean,PyramidalSD,PyramidalSEM,PyramidalN,PyramidalPValue,...
%     InterneuronMean,InterneuronSD,InterneuronSEM,InterneuronN,InterneuronPValue);
% 
% writetable(CA3FRpvals,['~/Documents/MATLAB/GenerateFigs/PDFstats/',statName])
%%
allIntensity = -sort(allIntensity);

pyrTotal = [];
intTotal = [];
%%
for nn = 1:length(allIntensity)
    ipV = interneuron.nullValue(interneuron.intensity == allIntensity(nn));
    iFR = interneuron.nullFR(interneuron.intensity == allIntensity(nn));
    
    iInc = length(find(ipV <= 0.01 & iFR >= 0));
    iDec = length(find(ipV <= 0.01 & iFR < 0));
    iNR = length(find(ipV > 0.01));
    
    ppV = pyramidal.nullValue(pyramidal.intensity == allIntensity(nn));
    pFR = pyramidal.nullFR(pyramidal.intensity == allIntensity(nn));
    
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
subplot(3,2,5)
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
text(22,.85,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(22,.75,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('Laser Intensity (mW)','fontsize',7)
ylabel('Proportion','fontsize',7)
title('Pyramidal Cells','fontsize',8,'FontWeight','normal')

subplot(3,2,6)
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
text(22,.85,'Decrease','color',clr(2,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')
text(22,.75,'NS','color',clr(3,:),'fontsize',6,'fontweight','bold','horizontalalignment','right')

ax = gca; 
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
xlabel('Laser Intensity (mW)','fontsize',7)
ylabel('Proportion','fontsize',7)
title('Interneurons','fontsize',8,'FontWeight','normal')
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
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 3, 4], 'PaperUnits', 'Inches', 'PaperSize', [4, 5])
print(gcf, '-dpdf',['~/Documents/MATLAB/GenerateFigs/PDFfigs/',figureName]);

%%
intensity_uW = [.001 .003 .01 .17 .75 3 6]';

[h1,PooOne,CI,statsp1] = vartest2(pyramidal.deltaFR(pyramidal.intensity == 190),pyramidal.nullFR(pyramidal.intensity == 190));
[h2,PooThree,CI,statsp2] = vartest2(pyramidal.deltaFR(pyramidal.intensity == 170),pyramidal.nullFR(pyramidal.intensity == 170));
[h3,PoOne,CI,statsp3] = vartest2(pyramidal.deltaFR(pyramidal.intensity == 150),pyramidal.nullFR(pyramidal.intensity == 150));
[h4,PoneSeven,CI,statsp4] = vartest2(pyramidal.deltaFR(pyramidal.intensity == 130),pyramidal.nullFR(pyramidal.intensity == 130));
[h5,PsevenFive,CI,statsp5] = vartest2(pyramidal.deltaFR(pyramidal.intensity == 110),pyramidal.nullFR(pyramidal.intensity == 110));
[h6,PThree,CI,statsp6] = vartest2(pyramidal.deltaFR(pyramidal.intensity == 90),pyramidal.nullFR(pyramidal.intensity == 90));
[h7,Psix,CI,statsp7] = vartest2(pyramidal.deltaFR(pyramidal.intensity == 70),pyramidal.nullFR(pyramidal.intensity == 70));

PyramidalFTest = [PooOne PooThree PoOne PoneSeven PsevenFive PThree Psix]';
PyramidalFStat = [statsp1.fstat statsp2.fstat statsp3.fstat statsp4.fstat ...
    statsp5.fstat statsp6.fstat statsp7.fstat]';
PyramidalDF = [statsp1.df1 statsp2.df1 statsp3.df1 statsp4.df1 ...
    statsp5.df1 statsp6.df1 statsp7.df1]';

[h1,IooOne,CI,statsi1] = vartest2(interneuron.deltaFR(interneuron.intensity == 190),interneuron.nullFR(interneuron.intensity == 190));
[h2,IooThree,CI,statsi2] = vartest2(interneuron.deltaFR(interneuron.intensity == 170),interneuron.nullFR(interneuron.intensity == 170));
[h3,IoOne,CI,statsi3] = vartest2(interneuron.deltaFR(interneuron.intensity == 150),interneuron.nullFR(interneuron.intensity == 150));
[h4,IoneSeven,CI,statsi4] = vartest2(interneuron.deltaFR(interneuron.intensity == 130),interneuron.nullFR(interneuron.intensity == 130));
[h5,IsevenFive,CI,statsi5] = vartest2(interneuron.deltaFR(interneuron.intensity == 110),interneuron.nullFR(interneuron.intensity == 110));
[h6,IThree,CI,statsi6] = vartest2(interneuron.deltaFR(interneuron.intensity == 90),interneuron.nullFR(interneuron.intensity == 90));
[h7,Isix,CI,statsi7] = vartest2(interneuron.deltaFR(interneuron.intensity == 70),interneuron.nullFR(interneuron.intensity == 70));

InterneuronFTest = [IooOne IooThree IoOne IoneSeven IsevenFive IThree Isix]';
InterneuronFStat = [statsi1.fstat statsi2.fstat statsi3.fstat statsi4.fstat ...
    statsi5.fstat statsi6.fstat statsi7.fstat]';
InterneuronDF = [statsi1.df1 statsi2.df1 statsi3.df1 statsi4.df1 ...
    statsi5.df1 statsi6.df1 statsi7.df1]';

f_test = table(intensity_uW,PyramidalDF,PyramidalFStat,PyramidalFTest,...
    InterneuronDF,InterneuronFStat,InterneuronFTest);

writetable(f_test,['~/Documents/MATLAB/GenerateFigs/PDFstats/',fStatName]);
%%
close all 
clear 
