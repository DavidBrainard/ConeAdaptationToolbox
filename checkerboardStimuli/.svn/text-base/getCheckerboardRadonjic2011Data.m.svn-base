function Checkerboard=getCheckerboardRadonjic2011Data(varargin)
% function Checkerboard=getCheckerboardRadonjic2011Data(varargin)
%
% experiment - 'one' or 'two'; Chooses data and arrangement from specified experiment
% ExpCol - Column number of experiment/condition data to be used 
% midSquareNum - matrix row number of the square that the test square should match 
% midSquareCol - Column of the experiment that the middle square's value
%       will be pulled from in getCheckerboardRadonjic2011Data. Usually the same
%       as the ExpCol
%
% Takes data directly from supplemental material in Radonjic 2011 and turns
% it into a checkerboard matrix of illuminance values that are picked from
% one of the columns of CheckIlluminanceData based on which
% experiment/condition is desired.
% 
% These are true luminance values in cd/m^2. 
%
% 8/7/12 ekf wrote it

%% Set defaults
experiment='two';
ExpCol=4;
midSquareNum=5;
midSquareCol=ExpCol;

%% Change values based on input
if ~isempty(varargin)
    if rem(length(varargin),2)~=0
        error('Arguments must be (pair, val) pairs');
    end
    MSCDefined=0;
    for ii=1:2:(length(varargin)-1)
        parm=varargin{ii};
        val=varargin{ii+1};
        switch parm
            case 'experiment'
                experiment = val;

            case 'ExpCol'
                ExpCol = val;

            case 'midSquareNum'
                midSquareNum = val;
                
            case 'midSquareCol'
                midSquareCol = val;
                MSCDefined = 1;
        end
    end
    if ~MSCDefined
        midSquareCol=ExpCol;
    end
end
%% Setup
numRows=5; %As specified in paper
centerSquare=ceil(numRows/2);

% Data copied from supplementary material of Radonjic 2001; first column is
% just numbering
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

% Get arrangement of square illumances according to paper's supplemental
% material
ExpOneArrangement=[ 7  23  4  10  12;
 18  6  22  17  21;
 1  13 NaN  3  14;
 19  15  16  8  20;
 5  9  2  24  11];

ExpTwoArrangement=[9  24  23  4  22;
 21  10  5  17  3;
 7  20 NaN  6  11;
 18  15  19  12  14;
 2  8  1  16  13];


switch experiment
    case 'one'
        % Experiment 1
        for col=1:numRows
            for row=1:numRows
                if row==centerSquare && col==centerSquare
                    Checkerboard(row,col)=CheckIlluminanceData(midSquareNum,midSquareCol); 
                else
                    Checkerboard(row,col)=CheckIlluminanceData(ExpOneArrangement(row,col),ExpCol); 
                end
            end
        end
    case 'two'
        % Experiment 2
        for col=1:numRows
            for row=1:numRows
                if row==centerSquare && col==centerSquare
                    Checkerboard(row,col)=CheckIlluminanceData(midSquareNum,midSquareCol); 
                else
                    Checkerboard(row,col)=CheckIlluminanceData(ExpTwoArrangement(row,col),ExpCol); 
                end
            end
        end
end




