

%   Objective: Producing 3D Shear wave motion data grid


clear;
close all;
clc;


% cd('N:\Faculty-of-Engineering\Research\I3S\Ultrasound\PhD\Safeer\Experiments_Data_&_codes\RF_data_and_Matlab\SWEI\UARPII\20_05_2016_CPWI_SWEI_EXPERIMENT');
% addpath(genpath('N:\Faculty-of-Engineering\Research\I3S\Ultrasound\PhD\Safeer\Experiments_Data_&_codes\RF_data_and_Matlab\SWEI\UARPII\20_05_2016_CPWI_SWEI_EXPERIMENT'));

cd('D:\UARPII Data\UARPII\UARP - DPB - Experiments\28_06_2016_CPWI_SWEI_EXPERIMENT');
addpath(genpath('D:\UARPII Data\UARPII\UARP - DPB - Experiments\28_06_2016_CPWI_SWEI_EXPERIMENT'));



%% Experiment parameters !!!

Exp.ImagingDepth = 50E-3;
Exp.FrameRate    = 10E3;
Exp.TimeZero     = 110; 
Exp.angles       = [-4 0 4];


%%  NEW Phantom

% fileName = 'UARP_TEMP_DATA_HETROII_NEW13.mat';
% ExpName  = 'PHANTOM-NEW-1';
% Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);
% 
% 
% fileName = 'UARP_TEMP_DATA_HETROII_NEWW13.mat';
% ExpName  = 'PHANTOM-NEW-2';
% Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);


%% Homo


fileName = 'UARP_TEMP_DATA_HOMO1.mat';
ExpName  = 'AP1-SF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO2.mat';
ExpName  = 'AP2-SF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO3.mat';
ExpName  = 'AP3-SF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO4.mat';
ExpName  = 'AP4-SF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO5.mat';
ExpName  = 'AP1-DF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO6.mat';
ExpName  = 'AP2-DF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO7.mat';
ExpName  = 'AP3-DF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO8.mat';
ExpName  = 'AP4-DF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO9.mat';
ExpName  = 'AP5-SF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO10.mat';
ExpName  = 'AP6-SF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO11.mat';
ExpName  = 'AP5-DF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO12.mat';
ExpName  = 'AP6-DF';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO13.mat';
ExpName  = 'U-CUSE';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO14.mat';
ExpName  = 'F-CUSE';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO15.mat';
ExpName  = 'SSI-LEFT';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HOMO16.mat';
ExpName  = 'SSI-RIGHT';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);



%% HETROI


fileName = 'UARP_TEMP_DATA_HETROI1.mat';
ExpName  = 'AP1-SF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI2.mat';
ExpName  = 'AP2-SF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI3.mat';
ExpName  = 'AP3-SF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI4.mat';
ExpName  = 'AP4-SF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI5.mat';
ExpName  = 'AP1-DF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI6.mat';
ExpName  = 'AP2-DF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI7.mat';
ExpName  = 'AP3-DF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI8.mat';
ExpName  = 'AP4-DF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI9.mat';
ExpName  = 'AP5-SF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI10.mat';
ExpName  = 'AP6-SF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI11.mat';
ExpName  = 'AP5-DF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI12.mat';
ExpName  = 'AP6-DF-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);



fileName = 'UARP_TEMP_DATA_HETROI14.mat';
ExpName  = 'F-CUSE-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI15.mat';
ExpName  = 'SSI-LEFT-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROI16.mat';
ExpName  = 'SSI-RIGHT-HETROI';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

%% HETROII


fileName = 'UARP_TEMP_DATA_HETROII1.mat';
ExpName  = 'AP1-SF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII2.mat';
ExpName  = 'AP2-SF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII3.mat';
ExpName  = 'AP3-SF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII4.mat';
ExpName  = 'AP4-SF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII5.mat';
ExpName  = 'AP1-DF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII6.mat';
ExpName  = 'AP2-DF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII7.mat';
ExpName  = 'AP3-DF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII8.mat';
ExpName  = 'AP4-DF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII9.mat';
ExpName  = 'AP5-SF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII10.mat';
ExpName  = 'AP6-SF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII11.mat';
ExpName  = 'AP5-DF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII12.mat';
ExpName  = 'AP6-DF-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);



fileName = 'UARP_TEMP_DATA_HETROII14.mat';
ExpName  = 'F-CUSE-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII15.mat';
ExpName  = 'SSI-LEFT-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

fileName = 'UARP_TEMP_DATA_HETROII16.mat';
ExpName  = 'SSI-RIGHT-HETROII';
Result = SWEI_AX_DISP_ESTIMATOR(fileName, ExpName, Exp);

