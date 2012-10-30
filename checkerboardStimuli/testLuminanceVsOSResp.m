% testLuminanceVsOSRep
%
% A script that calculates the cone outer segment response for a variety of
% luminance contexts. The data is taken from Radonjic 2011, and the plots
% parallel those of figure 3. The script is meant to show that there is a
% similar relationship between luminance context and percieved lightness
% and luminance context and  cone outer segment response. However, due to
% factors like the lack of feedback mechanism and relatively small optical
% blurring, the relationship is not as similar as hoped.

CheckIlluminanceData=...
[1 0.03 0.03 0.37 11.33 0.04 2.23 0.10 0.10 0.10 0.10;
2 0.05 0.07 0.48 13.08 0.06 2.53 0.13 0.12 0.16 0.11;
3 0.09 0.09 0.65 15.32 0.10 3.00 0.92 0.13 0.21 0.16;
4 0.11 0.13 0.90 18.04 0.13 3.57 1.35 0.18 0.27 0.23;
5 0.17 0.18 1.19 20.40 0.19 4.03 1.87 0.27 0.40 0.33;
6 0.25 0.26 1.61 23.73 0.25 4.63 2.78 0.35 0.55 0.46;
7 0.35 0.39 2.19 27.57 0.35 5.38 3.95 0.49 0.78 0.70;
8 0.56 0.56 2.94 31.95 0.46 6.36 5.71 0.73 1.08 1.00;
9 0.84 0.82 4.08 36.54 0.61 7.31 8.47 0.88 1.53 1.43;
10 1.25 1.30 5.25 43.05 0.87 8.32 12.19 1.06 2.15 2.00;
11 1.91 1.94 7.24 49.47 1.22 9.57 17.07 1.55 3.06 2.90;
12 2.95 2.93 9.50 57.67 1.67 11.28 24.94 1.77 4.30 4.21;
13 4.52 4.33 12.90 66.55 2.21 12.93 35.54 2.44 6.16 5.95;
14 6.76 6.41 17.23 76.85 3.05 15.03 50.91 3.60 8.85 8.64;
15 10.17 9.46 23.12 89.73 4.09 17.66 72.78 5.20 12.30 12.00;
16 15.27 14.31 31.63 105.44 5.56 20.14 108.56 7.26 17.20 16.84;
17 23.14 21.13 42.30 124.09 7.45 23.31 137.96 10.36 24.51 24.16;
18 34.91 31.49 56.98 146.83 10.14 26.96 151.72 14.52 34.55 34.12;
19 52.00 46.48 76.97 170.27 13.49 30.94 172.93 20.76 48.67 48.2111;
20 78.58 68.51 105.19 199.21 18.50 35.21 180.75 25.61 69.04 68.43;
21 119.77 106.02 148.58 236.83 25.11 42.21 186.54 36.04 98.36 97.55;
22 193.05 167.09 203.34 275.70 33.43 48.17 192.94 51.35 145.06 144.24;
23 295.48 257.38 289.60 325.15 45.12 55.09 228.66 214.59 211.23 210.37;
24 464.32 396.00 400.87 385.64 60.32 63.91 331.87 325.21 332.15 331.08];

% Preallocations
meanRespCol3=zeros(1,24);
testValCol3=zeros(1,24);
meanRespHigh1000=zeros(1,24);
testValHigh1000=zeros(1,24);
meanRespHigh30=zeros(1,24);
testValHigh30=zeros(1,24);
meanRespLow1000=zeros(1,24);
testValLow1000=zeros(1,24);
meanRespLow30=zeros(1,24);
testValLow30=zeros(1,24);
meanRespHighMean=zeros(1,24);
testValHighMean=zeros(1,24);
meanRespLowMean=zeros(1,24);
testValLowMean=zeros(1,24);

% Find indices of center check in the outer segment response matrix by
% running the first calculation outside the loop
    High1000Data=getRespNoPlots('experiment','two','ExpCol',4,'midSquareNum',1);
    [m,n]=size(High1000Data);
    squareSize=m/5-rem(m/5,1);
    edge=rem(m,5)/2;
    centerCheck=High1000Data(edge+2*squareSize:edge+3*squareSize,edge+2*squareSize:edge+3*squareSize);
    meanRespHigh1000(1)=mean(mean(centerCheck));
    testValHigh1000(1)=CheckIlluminanceData(1,4);
    
%% Graph A
% 10,000:1
for k=1:24
    Col3Data=getRespNoPlots('experiment','two','ExpCol',3,'midSquareNum',k);
    % Pulls out all the cones in the middle check of a 5 by 5 checkerboard
    centerCheck=Col3Data(edge+2*squareSize:edge+3*squareSize,edge+2*squareSize:edge+3*squareSize);
    meanRespCol3(k)=mean(mean(centerCheck));
    testValCol3(k)=CheckIlluminanceData(k,3);
end
% High 1,000:1 calculated in graph B
% High 30:1 calculated in graph B
%% Graph B
% High 1000:1
for k=2:24
    High1000Data=getRespNoPlots('experiment','two','ExpCol',4,'midSquareNum',k);
    % Pulls out all the cones in the middle check of a 5 by 5 checkerboard
    centerCheck=High1000Data(edge+2*squareSize:edge+3*squareSize,edge+2*squareSize:edge+3*squareSize);
    meanRespHigh1000(k)=mean(mean(centerCheck));
    testValHigh1000(k)=CheckIlluminanceData(k,4);
end
% Low 1000:1
for k=1:24
    Low1000Data=getRespNoPlots('experiment','two','ExpCol',6,'midSquareNum',k);
    % Pulls out all the cones in the middle check of a 5 by 5 checkerboard
    centerCheck=Low1000Data(edge+2*squareSize:edge+3*squareSize,edge+2*squareSize:edge+3*squareSize);
    meanRespLow1000(k)=mean(mean(centerCheck));
    testValLow1000(k)=CheckIlluminanceData(k,6);
end
%% Graph C
% High 30:1
for k=1:24
    High30Data=getRespNoPlots('experiment','two','ExpCol',5,'midSquareNum',k);
    % Pulls out all the cones in the middle check of a 5 by 5 checkerboard
    centerCheck=High30Data(edge+2*squareSize:edge+3*squareSize,edge+2*squareSize:edge+3*squareSize);
    meanRespHigh30(k)=mean(mean(centerCheck));
    testValHigh30(k)=CheckIlluminanceData(k,5);
end
% Low 30:1
for k=1:24
    Low30Data=getRespNoPlots('experiment','two','ExpCol',7,'midSquareNum',k);
    % Pulls out all the cones in the middle check of a 5 by 5 checkerboard
    centerCheck=Low30Data(edge+2*squareSize:edge+3*squareSize,edge+2*squareSize:edge+3*squareSize);
    meanRespLow30(k)=mean(mean(centerCheck));
    testValLow30(k)=CheckIlluminanceData(k,7);
end
%% Graph D
for k=1:24
    HighMeanData=getRespNoPlots('experiment','two','ExpCol',8,'midSquareNum',k,'midSquareCol',10);
    % Pulls out all the cones in the middle check of a 5 by 5 checkerboard
    centerCheck=HighMeanData(edge+2*squareSize:edge+3*squareSize,edge+2*squareSize:edge+3*squareSize);
    meanRespHighMean(k)=mean(mean(centerCheck));
    testValHighMean(k)=CheckIlluminanceData(k,10);
end
% Low 30:1
for k=1:24
    LowMeanData=getRespNoPlots('experiment','two','ExpCol',9,'midSquareNum',k,'midSquareCol',11);
    % Pulls out all the cones in the middle check of a 5 by 5 checkerboard
    centerCheck=LowMeanData(edge+2*squareSize:edge+3*squareSize,edge+2*squareSize:edge+3*squareSize);
    meanRespLowMean(k)=mean(mean(centerCheck));
    testValLowMean(k)=CheckIlluminanceData(k,11);
end

%% Plots
subplot(2,2,1); hold on
plot(log10(testValCol3), log10(meanRespCol3),'b')
plot(log10(testValHigh1000), log10(meanRespHigh1000),'k')
plot(log10(testValHigh30), log10(meanRespHigh30),'r')
xlabel('Log 10 Test Luminance')
ylabel('Log 10 Cone Outer Segment Response mV')
title('Graph A: Luminance-to-Lightness matching functions')
subplot(2,2,2); hold on
plot(log10(testValHigh1000), (meanRespHigh1000),'b')
plot(log10(testValLow1000),(meanRespLow1000), 'k')
xlabel('Log 10 Test Luminance')
ylabel('Log 10 Cone Outer Segment Response mV')
title('Graph B: Luminance-to-Lightness matching functions for 1000:1 luminance range context')
subplot(2,2,3); hold on
plot(log10(testValHigh30),log10(meanRespHigh30), 'b')
plot(log10(testValLow30),log10(meanRespLow30), 'k')
xlabel('Log 10 Test Luminance')
ylabel('Log 10 Cone Outer Segment Response mV')
title('Graph C: Luminance-to-Lightness matching functions for 30:1 luminance range context')
subplot(2,2,4); hold on
plot(log10(testValHighMean),log10(meanRespHighMean), 'b')
plot(log10(testValLowMean),log10(meanRespLowMean), 'k')
xlabel('Log 10 Test Luminance')
ylabel('Log 10 Cone Outer Segment Response mV')
title('Graph D: Luminance-to-Lightness matching functions for High vs Low mean luminance distribution')

lightnessData;

% Assumes that there is a function relationship between cone response and
% perceived lightness. We start by finding that function from one case,
% then transforming the other cases according to the function. The first function
% turns out to not apply to all cases. So either there is a relationship,
% but it is dependent on luminance/luminance context; or there is not yet a
% viable relationship possibly due to lack of feedback and small amounts of
% optical blurring.

meanRespCol3Interp=interp1(testValCol3,meanRespCol3,10.^A10000(:,1),'linear');

A10000Interp=interp1(10.^A10000(:,2),meanRespCol3Interp, 10.^A10000(:,2),'linear');
A1000Interp=interp1(10.^A10000(:,2),meanRespCol3Interp, 10.^A1000(:,2),'linear');
A30Interp=interp1(10.^A10000(:,2),meanRespCol3Interp, 10.^A30(:,2),'linear');
BHighInterp=interp1(10.^A10000(:,2),meanRespCol3Interp, 10.^BHigh(:,2),'linear');
BLowInterp=interp1(10.^A10000(:,2),meanRespCol3Interp, 10.^BLow(:,2),'linear');
CHighInterp=interp1(10.^A10000(:,2),meanRespCol3Interp, 10.^CHigh(:,2),'linear');
CLowInterp=interp1(10.^A10000(:,2),meanRespCol3Interp, 10.^CLow(:,2),'linear');
DHighInterp=interp1(10.^A10000(:,2),meanRespCol3Interp,10.^DHigh(:,2),'linear');
DLowInterp=interp1(10.^A10000(:,2),meanRespCol3Interp,10.^DLow(:,2),'linear');

figure;

% What the cone response should look like if the model fully explained
% perception
subplot(2,2,1); hold on
plot(A10000(:,1), log10(A10000Interp),'b')
plot(A1000(:,1), log10(A1000Interp),'k')
plot(A30(:,1), log10(A30Interp),'r')
xlabel('Log 10 Test Luminance')
ylabel('Log 10 Cone Outer Segment Response mV')
title('Graph A: Luminance-to-Lightness matching functions')
subplot(2,2,2); hold on
plot(BHigh(:,1), (BHighInterp),'b')
plot(BLow(:,1),(BLowInterp), 'k')
xlabel('Log 10 Test Luminance')
ylabel('Log 10 Cone Outer Segment Response mV')
title('Graph B: Luminance-to-Lightness matching functions for 1000:1 luminance range context')
subplot(2,2,3); hold on
plot(CHigh(:,1),log10(CHighInterp), 'b')
plot(CLow(:,1),log10(CLowInterp), 'k')
xlabel('Log 10 Test Luminance')
ylabel('Log 10 Cone Outer Segment Response mV')
title('Graph C: Luminance-to-Lightness matching functions for 30:1 luminance range context')
subplot(2,2,4); hold on
plot(DHigh(:,1),log10(DHighInterp), 'b')
plot(DLow(:,1),log10(DLowInterp), 'k')
xlabel('Log 10 Test Luminance')
ylabel('Log 10 Cone Outer Segment Response mV')
title('Graph D: Luminance-to-Lightness matching functions for High vs Low mean luminance distribution')


