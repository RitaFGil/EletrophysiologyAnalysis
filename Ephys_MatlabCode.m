%% Electrophysiology Data Analysis
% Written by Rita Gil
% Last modified 15.12.2023

% PDS calculation per animal:
% Load Data
% Notch filter to remove power line noise (and harmonics)
% PSD calculation per animal per recording

% Combine animals and calculate final Power spectrum for LFP and MUA frequency ranges.
% Calculate fina Spectrogram for each group

% Steady-state boxplots for MUA and LFP

% Correlation of steady-state between MUA/LFP and fMRI 

% Electrophysiology Neurometric curve:
% Load PSD peak integral value per trial
% Categorize trials as a 'pulsating' or 'continuous' light report (done for early and late steady-state intervals)
% Fit psichometric curve
% Boostraping with replacement for final threshold calculation (meand and standard deviation)


%% Load data and remove power line noise and harmonics

clear all; clc; close all;
filename='path'; %path where raw data from Open Ephys is
D=load_open_ephys_binary(filename,'continuous',1,'mmap');
Data.raw=D.Data.Data(1).mapped;
clear D;

% Separation of different frequency conditions
Fs=30000; %sampling frequency Hz
numChannels=64; %number of channels in the probe
numShanks=8; %number of probe shanks
numCycles=10; %number of cycles per condition
stimTime=15; %Stimulation time in sec
restTime=45; %Resting period in sec

%Condition X
frequency=X; %Frequency in Hertz
timePulseX=1/frequency; %Period for a specific frequency
durX=Fs*timePulseX; %Period for a specific frequency in samples
numpulseX=Y;%number of light pulses during stimulation period

%Condition Y (following condition X)
frequency=X; %Frequency in Hertz
timePulseY=1/frequency; %Period for a specific frequency
durY=Fs*timePulseY; %Period for a specific frequency in samples
numpulseY=Y;%number of light pulses during stimulation period

%probe map check
orderVec=[48,35,44,46,43,38,40,39,42,54,33,34,37,52,50,36,49,63,56,55,51,61,58,53,57,41,59,...
    47,64,45,62,60,23,7,17,5,19,2,6,4,1,15,9,10,3,13,11,8,12,24,32,31,14,27,30,16,29,18,20,22,28,21,25,26];

%Creating stimulation vector and resizing it
stim_vec=double(Data.raw(67,:)); %67 is the trigger channel
stim_vec(stim_vec<2000)=0;
stim_vec(stim_vec~=0)=3000;
[PKS,LOCS] = findpeaks(stim_vec);

%Sample vector for each condition
vecX=LOCS(1)-(restTime*Fs):LOCS(numCycles)+((stimTime+restTime)*Fs);
vecY=LOCS(numCycles+1)-(restTime*Fs):LOCS(2*numCycles)+((stimTime+restTime)*Fs);

%stimulation vector for each frequency condition
stim_vecX=double(stim_vec(LOCS(1)-(restTime*Fs):LOCS(numCycles)+((stimTime+restTime)*Fs)));
stim_vecY=double(stim_vec(LOCS(numCycles+1)-(restTime*Fs):LOCS(2*numCycles)+((stimTime+restTime)*Fs)));

cd('path') %path to save data
save('workspace.mat','filename','Fs','numChannels','numShanks','numCycles','restTime','stimTime','stim_vec','timePulseX','timePulseY',...
    'vecX','vecY','durX','durY','numpulseX','numpulseY','stim_vecX','stim_vecY','orderVec')

% Separation of conditions and applying probe map
DataX.raw=Data.raw(orderVec,vecX);
DataY.raw=Data.raw(orderVec,vecY);
clear Data

% Apply notch filter per channel: removing 50Hz, 100Hz and 150Hz --> power line and its harmonics
DataX.notch=[];
DataY.notch=[];
for i=1:numChannels
    
    DataX.notch(i,:)= detrend(donotch(double(DataX.raw(i,:)),Fs));
    DataY.notch(i,:)= detrend(donotch(double(DataY.raw(i,:)),Fs));
    
end

%save data
save('Datanotch.mat','DataX','DataY''-v7.3');

%% Calculating PSD per animal and per recording
clear all; clc; close all;

animal=[1,2,3,4]; %vector with animal numbers

Fs=30000; %sampling frequency Hz
numChannels=64; %number of channels in the probe
numShanks=8; %number of probe shanks
numCycles=10; %number of cycles per condition
stimTime=15; %Stimulation time in sec
restTime=45; %Resting period in sec

%Condition X
frequency=X; %Frequency in Hertz
timePulseX=1/frequency; %Period for a specific frequency
durX=Fs*timePulseX; %Period for a specific frequency in samples
numpulseX=X;%number of light pulses during stimulation period

for an=animal
    an
    if an==XX
        recording=[1,2,3];
    elseif an==YY || an==ZZ
        recording=[1,2,3,4];
    end
    
    for reco=recording
        reco
        
        clear DataXTMP DataXNEW
        
        cd('path'); %path where raw data is saved
        
        load('workspace.mat', 'stim_vec','stim_vecx')
        load('DataX.mat')
        
        DataXN.notch=DataX.notch;
        clear DataX
        
        % Move to new folder where processed data is saved
        cd('path');
        save('workspace.mat')
        
        %Break data into cycles
        %Calculating timepoints of begining and end of separate cylces (stimulation period centered)
        [PKS,LOCS] = findpeaks(stim_vecx);
        
        clear PKS vecLOC
        i=1;
        
        vecLOC.XHzbegin(i)=LOCS(i)-(restTime/2)*Fs;
        vecLOC.XHzend(i)=LOCS(i)+(stimTime+restTime/2)*Fs;
        difference=vecLOC.XHzend(i)-vecLOC.XHzbegin(i);
        
        for i=1:numCycles
            [PKS,LOCS] = findpeaks(stim_vecX);
            vecLOC.XHzbegin(i)=LOCS(i)-(restTime/2)*Fs;
            vecLOC.XHzend(i)=vecLOC.XHzbegin(i)+difference;
        end
        
        DataXtN.notchCycle=[];
        
        %Sparation of individual cycles
        for s=1:numChannels
            
            for c=1:numCycles
                DataXN.notchDetrendCycle(s,c,:)=DataXN.notch(s,vecLOC.XHzbegin(c):vecLOC.XHzend(c));
            end
        end
        
        %Stimulation vector with stimulation period centered
        stimvecFINALtmp=zeros(1,length(DataXN.notchDetrendCycle));
        stimvecFINALtmp(((restTime/2)*Fs)+1:((restTime/2+stimTime)*Fs)+1)=1;
        save('stimvec.mat','stimvecFINALtmp')
        
        DataXTMP.notchDetrendCycle=DataXN.notchDetrendCycle;
        clear DataXN
        
        %Give each channel a specific weight (different for groups with low frequency regimes and groups witn only high frequency regimes)
        
        %_______________________________________________
        
        %LOW frequencies where there is increased power during stim
        
        DataXTMP.AvCyclePerChannel=squeeze(mean(DataXTMP.notchDetrendCycle,2));
        
        for s=1:numChannels
            fraction(s)=trapz(abs(DataXTMP.AvCyclePerChannel(s,find(stimvecFINALtmp==1)))); %based on the integral of power measured during stimulation
        end
        
        %HIGH frequencies where there is decreased power during stim
        DataXTMP.AvCyclePerChannel=squeeze(mean(DataXTMP.notchDetrendCycle,2));
        
        tpoint=find(stimvecFINALtmp==1);
        fistSecStim=tpoint(1)+2*Fs;
        TwoHalfAfterStim=tpoint(end)+2.5*Fs;
        
        for s=1:numChannels
            
            if prctile(abs(DataXTMP.AvCyclePerChannel(s,[1:tpoint(1)-1,TwoHalfAfterStim:end])),95)>prctile(abs(DataXTMP.AvCyclePerChannel(s,[tpoint(1):fistSecStim,tpoint(end):TwoHalfAfterStim])),95)
                fraction(s)=0;
                multiplier(s)=0;
            else
                multiplier(s)=prctile(abs(DataXTMP.AvCyclePerChannel(s,[tpoint(1):fistSecStim,tpoint(end):TwoHalfAfterStim])),100)-prctile(abs(DataXTMP.AvCyclePerChannel(s,[1:tpoint(1)-1,TwoHalfAfterStim:end])),50);
                fraction(s)=trapz(abs(DataXTMP.AvCyclePerChannel(s,[tpoint(1):fistSecStim,tpoint(end):TwoHalfAfterStim])))*multiplier(s);%based on the integral of the onset and offset peaks (only when these are measured)
                
            end
        end
        %_______________________________________________
        
        %calculated the weight each channel will have from now on
        for s=1:numChannels
            weight(s)=fraction(s)./sum(fraction(:));
        end
        save('workspace.mat','weight','-append')
        
        
        %Weighted combination of different channels
        for s=1:numChannels
            
            for c=1:numCycles
                
                DataXTMP.WnotchDetrendCycle(s,c,:)=DataXTMP.notchDetrendCycle(s,c,:)*weight(s);
                
            end
        end
        DataXTMP.notchAvCycle=squeeze(sum(DataXTMP.WnotchDetrendCycle,1));
        
        %Look at the LFP frequency range
        freqrangeLFP=[0.00001 150];
        
        for i=1:numCycles
            i
            DataXTMP.WCycleLFP(i,:)= bandpassOSL(DataXTMP.notchAvCycle(i,:),freqrangeLFP,Fs);
        end
        
        %Z-score
        DataXTMP.ZWCycleLFP=[];
        for s=1:numCycles
            s
            DataXTMP.ZWCycleLFP(s,:)=(DataXTMP.WCycleLFP(s,:)-mean(DataXTMP.WCycleLFP(s,1:LOCS-1)))...
                ./std(DataXTMP.WCycleLFP(s,1:LOCS-1));
        end
        save('DatacontTMP1.mat','DataXTMP','-v7.3')
        
        %PSD calculation during rest and stimulaton steady state
        windur=5; %Duration of time window in sec
        start4rest=10; %Start of window for the rest period starting at second 10.
        start1=(restTime/2)+2; %Start of window for the first stimulation steady-state period (2sec after stimulus started)
        start2=(restTime/2)+8; %Start of window for the first stimulation steady-state period (8sec after stimulus started)
        
        %Transform intervals into samples
        window2LookRest=[start4rest*Fs:(start4rest+windur)*Fs];
        window2LookStimSST1=[start1*Fs:(start1+windur)*Fs];
        window2LookStimSST2=[start2*Fs:(start2+windur)*Fs];
        
        L=length((window2LookStimSST1));
        winwin=hamming(L,'periodic');
        
        %Calculate power spectrum per cycle
        for c=1:numCycles
            
            A=(DataXTMP.ZWCycleLFP(c,:));
            [FFT_A1,f] =pwelch(A(window2LookRest),winwin,0,[],Fs,'power');
            [FFT_A2,f] =pwelch(A(window2LookStimSST1),winwin,0,[],Fs,'power');
            [FFT_A3,f] =pwelch(A(window2LookStimSST2),winwin,0,[],Fs,'power');
            
            %taking integral around the stimulation frequency (index1 and
            %index2 are 1Hz below and above the target frequency) for the
            %three time intervals
            [val1,idx1]=min(abs(f-index1));
            [val2,idx2]=min(abs(f-index2));
            
            IntegralVal(c,1)=trapz(2*abs(FFT_A2(idx1:idx2)));
            IntegralVal(c,2)=trapz(2*abs(FFT_A3(idx1:idx2)));
            IntegralVal(c,3)=trapz(2*abs(FFT_A1(idx1:idx2)));
            
        end
        
        %Average of all cycles
        A=(mean(DataXTMP.ZWCycleLFP,1));
        [FFT_A1,f] =pwelch(A(window2LookRest),winwin,0,[],Fs,'power');
        [FFT_A2,f] =pwelch(A(window2LookStimSST1),winwin,0,[],Fs,'power');
        [FFT_A3,f] =pwelch(A(window2LookStimSST2),winwin,0,[],Fs,'power');
        
        IntegralValAv(1,1)=trapz(2*abs(FFT_A2(idx1:idx2)));
        IntegralValAv(1,2)=trapz(2*abs(FFT_A3(idx1:idx2)));
        IntegralValAv(1,3)=trapz(2*abs(FFT_A1(idx1:idx2)));
        
        save(['DataXPowerSpec'])
        
        %Spectrogram
        overlp=0.25; %Percentage of overlap between sliding windows
        tempres=0.05; %final temporal resolution
        time = tempres/(1-overlp); %time of window in sec
        window=ceil(time*Fs); %Window in samples
        noverlap=round(overlp*window);
        
        freqvec=[0.1:0.1:50]; %LFP frequency range
        freqvec1=[0.5:0.5:3000]; %MUA frequency range
        
        %Using weighted channel data to create the average run response
        DataXTMP.AvChannotchAvCycle=mean(DataXTMP.notchAvCycle,1);
        
        close all
        %cont light
        clear x X Xnew
        x=double(DataXTMP.AvChannotchAvCycle);
        X = buffer(x,window,noverlap);
        Xnew=X(:,2:size(X,2)-1);
        [s11,f1] =pwelch(Xnew,hamming(window,'periodic'),noverlap,freqvec,Fs); %LFP
        [s12,f2] =pwelch(Xnew,hamming(window,'periodic'),noverlap,freqvec1,Fs); %MUA
        
        t1=((tempres:tempres:size(Xnew,2)*tempres));
        
        %stimvec after binning
        stimvecFINAL=zeros(1,length(s11));
        stimvecFINAL((restTime/2)/tempres+1:((restTime/2)+stimTime)/(tempres))=1;
        [PKS,LOCS] = findpeaks(stimvecFINAL);
        save('stimvec.mat','stimvecFINAL')
        
        DataXNEW.Spec01=s11;
        DataXNEW.Spec02=s12;
        
        %Save data
        save('DataXNEW1.mat','DataXNEW','-v7.3')
        save('spectrogram.mat','time','window','noverlap','Fs','f1','f2','t1','freqvec','freqvec1','overlp');
        
    end
end

%% Averaging animals

close all; clear all clc

cd('path'); %path where one example dataset is salved
load workspace -regexp ^(?!overlp$|t1$|time$|noverlap$|window$|ind1$|ind2$|ind3$|indLFP$|indMUA1$|indMUA2$).
load('spectrogram.mat')

filename='path'; %path where data is saved
numChannels=64;

rats={'Rat1','Rat2','Rat3','Rat4'}; %Animals used

%Getting total number of recordings for the condition being studied
numbtotalFiles=0;
for rat=1:length(rats)
    path2load=[filename,rats{rat},filesep,area{1},filesep,'Record Node 113',filesep,'continuous'];
    files = dir(path2load);
    numbtotalFiles=numbtotalFiles+(length(files)-2);
end

Data.Spec1=zeros(length(freqvec) , length(t1) , numbtotalFiles);
Data.Spec2=zeros(length(freqvec1) , length(t1) , numbtotalFiles);

%load data
clear path2load
r=1;
for rat=1:length(rats)
    
    path2loadinc=[path];
    files = dir(path2loadinc);
    
    for record=1:length(files)-2
        
        load([path2load,'DataXNEW1.mat'])
        
        Data.Spec1(:,:,r)=DataXNEW.Spec01;
        Data.Spec2(:,:,r)=DataXNEW.Spec02;
        
        clear DataXNEW1;
        r=r+1;
        
    end
end

%save data
save('AverageDataX.mat','Data','numbtotalFiles','-v7.3');

animal=[1,2,3,4]; %list of animals in a specific group

freqvec=[0.1:0.1:50]; %LFP range
freqvec1=[0.5:0.5:3000]; %MUA range
freqvec = round(freqvec,2);
freqvec1 = round(freqvec1,2);
tempres=0.05;%final temporal resolution

ind1= find(freqvec==50);
indLFP= find(freqvec1==150);
indMUA1= find(freqvec1==300);
indMUA2= find(freqvec1==3000);

load workspace -regexp ^(?!overlp$|t1$|time$|noverlap$|window$|ind1$|ind2$|ind3$|indLFP$|indMUA1$|indMUA2$).
load('spectrogram.mat')


xvectmp=[t1(1):tempres:t1(end)]; % convert xaxis to seconds
xvec=[-restTime/2+(tempres):tempres:restTime/2+stimTime-(tempres)];
xvec = round(xvec,2);

%Calculate ower in the desired frequency bands
Data.Spec2LFP=squeeze(trapz(abs(Data.Spec2(1:indLFP,:,:)),1)/indLFP);
Data.Spec2MUA=squeeze(trapz(abs(Data.Spec2(indMUA1:indMUA2,:,:)),1)/(indMUA2-indMUA1));

%Z-score
baselineind=find(xvec==0);
for f=1:numbtotalFiles
    f
    Data.Spec2LFPAvZ(:,f)=(Data.Spec2LFP(:,f)-mean(Data.Spec2LFP(1:baselineind,f),1))./(std(Data.Spec2LFP(1:baselineind,f),0,1));
    Data.Spec2MUAAvZ(:,f)=(Data.Spec2MUA(:,f)-mean(Data.Spec2MUA(1:baselineind,f),1))./(std(Data.Spec2MUA(1:baselineind,f),0,1));
    
end

%Mean +/- s.e.m. runs of different animals
Data.Spec2LFPAv=squeeze(mean(Data.Spec2LFPAvZ,2));
Data.Spec2LFPSEM=squeeze(std(Data.Spec2LFPAvZ,0,2)./sqrt(length(animal)));
Data.Spec2MUAAv=squeeze(mean(Data.Spec2MUAAvZ,2));
Data.Spec2MUASEM=squeeze(std(Data.Spec2MUAAvZ,0,2)./sqrt(length(animal)));

%Plot results
close all
figure('units','normalized','outerposition',[0 0 1 1])
subplot(121)
shadedErrorBar(xvec,Data.Spec2LFPAv,Data.Spec2LFPSEM,'lineprops',{'color','b'})
title('cont LFP')
xlim([xvec(1) xvec(end)])
ylim([c d])

subplot(122)
shadedErrorBar(xvec,Data.Spec2MUAAv,Data.Spec2MUASEM,'lineprops',{'color','b'})
title('cont MUA')
xlim([xvec(1) xvec(end)])
ylim([c d])

%Save figure
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'PowerFinalPlots','epsc')
saveas(gcf,'PowerFinalPlots','fig')
saveas(gcf,'PowerFinalPlots','png')

Data.Spec1AvFold=squeeze(mean(abs(Data.Spec1),3));

close all
figure('units','normalized','outerposition',[0 0 1 1])
imagesc(abs(Data.Spec1AvFold(:,:)),[0 2e4])
colormap turbo;

%Save figure
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'FinalSpec','epsc')
saveas(gcf,'FinalSpec','fig')
saveas(gcf,'FinalSpec','png')

%save data
save('DataAvFinal.mat','Data','-v7.3')
save('workspace.mat','ind1','indLFP','indMUA1','indMUA2','xvec','-append')


%% Boxplots

close all
load('stimvec.mat')
load('spectrogram.mat')

tempres=0.05; %final emporal resolution
TR=time*overlp; %temporal resolution
stimdur=stimTime/TR;
xvectmp=[t1(1):tempres:t1(end)];
xvec=[-restTime/2+(tempres):tempres:restTime/2+stimTime-(tempres)]; %vector where 0 is where stimulation starts
xvec = round(xvec,2);

stimvecFINAL=zeros(size(xvec));
baselineind=find(xvec==0);
ind1sec=find(xvec==1); %One second after simulation started
ind15sec=find(xvec==15); %End of stimulation period
cd('path'); %path where data is stored

for f=1:size(Data.Spec2MUAAvZ,2)
    f
    for i=1:stimfreqs
        
        Data.MUAmeansst(f,i)=nanmedian(DataAv.Spec2MUAAvZ(ind1sec+1:ind15sec,f,i),1);%median values during steady-state (one second ...
        %after stimulation started until the end of the stimulation period)
        Data.LFPmeansst(f,i)=nanmedian(DataAv.Spec2LFPAvZ(ind1sec+1:ind15sec,f,i),1);
        
    end
end


%Save data
save('Data.mat','Data','-v7.3')

%% plot
load('stimvec.mat')

close all
figure('units','normalized','outerposition',[0 0 1 1])
subplot(121)
plot(xvec,DataAv.Spec2LFPAvZ)
hold on
plot(xvec,stimvecFINAL*100,'--','Color','r')
legend('0.25 Hz','1 Hz','2 Hz','15 Hz','20 Hz','25 Hz','cont')
ylim([-2 22])
xlim([xvec(1) xvec(end)])
title('Average LFP')

subplot(122)
plot(xvec,DataAv.Spec2MUAAvZ)
hold on
plot(xvec,stimvecFINAL*100,'--','Color','r')
legend('0.25 Hz','1 Hz','2 Hz','15 Hz','20 Hz','25 Hz','cont')
ylim([-2 20])
xlim([xvec(1) xvec(end)])
title('Average MUA')


%Save figure
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'Powerplots','epsc')
saveas(gcf,'Powerplots','fig')
saveas(gcf,'Powerplots','png')

%Box plots
close all
%LFP
nanmean(DataAv.LFPmeansst,1)
nanstd(DataAv.LFPmeansst,0,1)
figure('units','normalized','outerposition',[0 0 1 1])
boxplot(DataAv.LFPmeansst,'Labels',{'025 Hz','1Hz','2Hz','15Hz','20Hz','25Hz','cont'})
hold on
notBoxPlot(DataAv.LFPmeansst,'style','line','interval','tInterval','markMedian',true)
ylim([-2 10])
%save boxplots
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'SST_BoxplotLFP','epsc')
saveas(gcf,'SST_BoxplotLFP','fig')
saveas(gcf,'SST_BoxplotLFP','png')

%MUA
nanmean(DataAv.MUAmeansst,1)
nanstd(DataAv.MUAmeansst,0,1)
figure('units','normalized','outerposition',[0 0 1 1])
boxplot(DataAv.MUAmeansst,'Labels',{'025 Hz','1Hz','2Hz','15Hz','20Hz','25Hz','cont'})
hold on
notBoxPlot(DataAv.MUAmeansst,'style','line','interval','tInterval','markMedian',true)
ylim([-2 10])
%save boxplots
fig2=figure(2);
fig2.Renderer='Painters';
saveas(gcf,'SST_BoxplotMUA','epsc')
saveas(gcf,'SST_BoxplotMUA','fig')
saveas(gcf,'SST_BoxplotMUA','png')

%% Correlation of MUA with fMRI (repeat for LFPs)

load('SteadyState_fMRIData.mat') %load fMRI data for each specific ROI

MUAData=Data.MUAmeansst;
fMRIData=nan(length(MUAData),stimfreqs);
fMRIData(1:length(fMRIData),:)=ROIX_fMRIsteadystateData;


MEANfMRI=[nanmean(fMRIData,1)];
MEANMUA=[nanmean(MUAData,1)];

STDfMRI=[nanstd(fMRIData,0,1)];
STDMUA=[nanstd(MUAData,0,1)];


%Plot the x and Y error bars
close all
figure('units','normalized','outerposition',[0 0 1 1])

for i=2:7 %from 1Hz until continous light (0.25Hz is not included)
    errorbarxy(MEANfMRI(i), MEANMUA(i), STDfMRI(i), STDMUA(i))
    hold on
    pause
end

%Define axis-limits based on the error bars
xlim([a b])
ylim([c d])
%save figure
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'stdXY','epsc')


%Pearson correlation
close all
figure('units','normalized','outerposition',[0 0 1 1])
[RHO,PVAL,H] = corrplot([MEANfMRI(2:7)',MEANMUA(2:7)'],'Type','Pearson') %manually modify x and y axis to correspond to the first plot with x and y error bars
[R,P,RL,RU] = corrcoef([MEANfMRI(2:7)',MEANMUA(2:7)']);

%sava firue
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'Correlation','epsc')

%Resulting image will be an overlap of 'ROIX_corrCloud' and 'ROIX_stdXY'

%%
%%
%% Electrophysiology neurometric curve

clear all; clc; close all;

animal=[1,2,3,4,5]; %animals used in the study

Fs=30000; %sampling frequency Hz
numChannels=64; %number of channels in the probe
numShanks=8; %number of probe shanks
numCycles=10; %number of cycles per condition
stimTime=15; %Stimulation time in sec
restTime=45; %Resting period in sec
numfreqs=X; %number of frequencies tested

intervals=2; %Steady state intervals
windur=5; %durantion in s of the steady-state window
maxrecordings=X;%max number of recordings across animals

%Condition X
frequency=X; %Frequency in Hertz
timePulseX=1/frequency; %Period for a specific frequency
durX=Fs*timePulseX; %Period for a specific frequency in samples
numpulseX=X;%number of light pulses during stimulation period

%initialize vector where data is going to be put
vec4DesicionTrial=nan(numCycles,maxrecordings,length(animal),numfreqs,intervals);

ancount=1;
for an=animal
    
    if an==1 || an==2
        recording=[1,2,3];
        
    elseif  an==3 || an==5
        recording=[1,2,3,4];
        
    elseif an==4
        recording=[1,2,3,4,5];
        
    end
    
    for reco=recording
        reco
        
        cd([path]) %path where PSD data is stored
        load(['DataXPowerSpec.mat'],'peakVal','IntegralVal','peakValAv','IntegralValAv')
        
        vec4DesicionTrial(:,reco,ancount,XX,:)=IntegralVal(:,:,1:2)-IntegralValAv(:,3)';
        
    end
    
    ancount=ancount+1;
end

%Reshape data
reshplotT=reshape(vec4DesicionTrial,[size(vec4DesicionTrial,1)*size(vec4DesicionTrial,2)*size(vec4DesicionTrial,3),numfreqs,intervals,2]);

vec4DesicionTrialRESH=reshape(vec4DesicionTrial,[numCycles*size(vec4DesicionTrial,2),length(animal),numfreqs,intervals]);
DecisionTrial=nan(numCycles*size(vec4DesicionTrial,2),length(animal),numfreqs,intervals);
PercPulsedLight_Trial=nan(length(animal),numfreqs,intervals);

condref=numfreqs-3; %defining a reference condition (this was chosen based on the behavioural results)

%Defining a threshold, using the reference condition for the "continuous" or "pulsed" categorization
for int=1:intervals
    thresIntegralT(int)=round(nanmean(reshplotT(:,condref,int)));
end

%Categorization of each trial as a 'pulsating' or 'continuous light
%'report' according to the defined thresolhd
for int=1:intervals
    ancount=1;
    for an=animal
        for f=1:numfreqs
            for c=1:size(vec4DesicionTrialRESH,1)
                
                if vec4DesicionTrialRESH(c,ancount,f,int)>thresIntegralT(int) && isnan(vec4DesicionTrialRESH(c,ancount,f,int))==0
                    DecisionTrial(c,ancount,f,int)=1;
                elseif vec4DesicionTrialRESH(c,ancount,f,int)<thresIntegralT(int) && isnan(vec4DesicionTrialRESH(c,ancount,f,int))==0
                    DecisionTrial(c,ancount,f,int)=2;
                end
                
            end
            
            %Calculating the % of categorized reports to the pulating port
            PercPulsedLight_Trial(ancount,f,int)=(sum(DecisionTrial(:,ancount,f,int)==1)*100)./(sum(~isnan(DecisionTrial(:,ancount,f,int))));
            
        end
        ancount=ancount+1;
    end
end

%Average and s.e.m. across animals
AvPercPulsedLight_Trial=squeeze(nanmean(PercPulsedLight_Trial,1));

for int=1:intervals
    for i=1:numfreqs
        STDPercPulsedLight_Trial(i,int)=squeeze(nanstd(PercPulsedLight_Trial(:,i,int),0,1)./sqrt(length(animalcond(i))));
    end
end

xvecfreq=[1,2,8,12.5,15,20,25,40,60]';%vector with tested frequencies in Hz. The continuous light is represneted as the 60Hz

close all
for int=1:intervals
    
    %fitting a psychometric curve to the data (used in Wichmann, Felix A., and N. Jeremy Hill., Perception & psychophysics, 2001)
    ymeans=AvPercPulsedLight_Trial(2:end,int)/100; %First frequency regime (0.25Hz) is not considered for the fitting.
    [xData, yData]=prepareCurveData(xvecfreq,ymeans);
    q1=fittype(@(a,b,c,d,x) d + (c-d)./(1+exp(-(2*a)*(x-b))),'independent',{'x'},'dependent',{'y'});
    [f,~]=fit(xData,squeeze(yData),q1,'Lower', [0 xvecfreq(2) 0 0], 'Upper', [Inf xvecfreq(end) 1 1],'startpoint',[1,0,0.05, 0.95]);
    thresh=f.b-(1/(2*f.a))*log(((f.c-f.d)/(0.50-f.d))-1); %Calculation of threshold of 50% choices to pulsating port ('chance level')
    
    %plot results
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(xvecfreq,ymeans,'-o')
    hold on
    plot(f)
    ylim([0 1])
    hold on
    notBoxPlot(PercPulsedLight_Trial(:,2:end,int)/100,xvecfreq,'style','line')
    legend off
end

%save figures
fig1=figure(1);
fig1.Renderer='Painters';
saveas(gcf,'EphysNeurometricCurve_interval1','epsc')
saveas(gcf,'EphysNeurometricCurve_interval1','fig')
saveas(gcf,'EphysNeurometricCurve_interval1','png')

fig2=figure(2);
fig2.Renderer='Painters';
saveas(gcf,'EphysNeurometricCurve_interval2','epsc')
saveas(gcf,'EphysNeurometricCurve_interval2','fig')
saveas(gcf,'EphysNeurometricCurve_interval2','png')

%Bootstrapp with replacement for standard deviation calculation
it=50; %iterations for bootstrapping

for int=1:intervals
    
    d=PercPulsedLight_Trial(:,2:end,int); %First frequency regime (0.25Hz) is not considered for the fitting.
    bootstat  = bootstrp(it,@nanmean,d);
    bootstat(isnan(bootstat))=0;
    
    clear xData yData ymeans1 threshBOOT d
    for a=1:it
        
        %fitting a psychometric curve to the data (used in Wichmann, Felix A., and N. Jeremy Hill., Perception & psychophysics, 2001)
        ymeans1(a,:)=bootstat(a,:)/100;
        [xData, yData(a,:)]=prepareCurveData(xvecfreq,ymeans1(a,:)');
        q1=fittype(@(a,b,c,d,x) d + (c-d)./(1+exp(-(2*a)*(x-b))),'independent',{'x'},'dependent',{'y'});
        [f,~]=fit(xData,squeeze(yData(a,:))',q1,'Lower', [0 xvecfreq(1) 0 0], 'Upper', [Inf xvecfreq(end) 1 1],'startpoint',[1,0,0.05, 0.95]);
        threshBOOT(a)=f.b-(1/(2*f.a))*log(((f.c-f.d)/(0.5-f.d))-1); %Calculation of threshold of 50% choices to pulsating port ('chance level')
    end
    
    finalthreshhighSTDBOOT(int)=std(threshBOOT,0,2);
    
end

%save data
