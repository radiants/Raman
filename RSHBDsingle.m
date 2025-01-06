% Hydrogen Bond Decoupling for Raman Spectroscopy (HBDRS)
% Smothing with Wavelet method and Baseline Correction with Quadratic Polynomial Fitting
% Promoted by Xian-ze MENG, Fa-He Cao*, 28 Dec 2024, Sun Yat-sen University

clear
close all
% global SpectraType
% Load Raman Spectra file in .csv or .xlsx types
FileName = uigetfile('*.csv;*.xlsx');    % Click this button to load data

tic
RamanData = readtable(FileName); % disp(FileName) % from disk in xlsx or csv format
RamanshiftOri = table2array(RamanData(:,1)); % x axis of Raman Intensity % change Wavenumber formate to double
IntensityOri = table2array(RamanData(:,2)); % y axis of Raman Intensity % change Intensity formate to double
IntensityHeightOri = height(RamanshiftOri);

SpectraType = "Wavenumber"; % Chose the Spectrum type
switch SpectraType
    case "Wavenumber" % from 617.399 to 675.813 nm (2600.011275 to 3999.999822 cm-1);
        [~,edge1] = min(abs(RamanshiftOri - 617.399)); % find the data series1 (2600 cm-1);
        [~,edgeA] = min(abs(RamanshiftOri - 625.118)); % find the data seriesA (2800 cm-1);
        [~,edgeB] = min(abs(RamanshiftOri - 666.800)); % find the data seriesB (3800 cm-1);
        [~,edge2] = min(abs(RamanshiftOri - 675.813)); % find the data series2 (4000 cm-1);
        RamanshiftCut = RamanshiftOri(edge1:edge2);
        IntensityCut = IntensityOri(edge1:edge2); % new Raman Intensity;
        FittingCut = RamanshiftOri([edge1:edgeA, edgeB:edge2]);
        FittingInt = IntensityOri([edge1:edgeA, edgeB:edge2]);
        RamanEdge = RamanshiftOri(edgeA:edgeB);
        IntensityEdge = IntensityOri(edgeA:edgeB);
        NumCut = length(IntensityCut);
        NumEdge = length(IntensityEdge);
    case "RamanShift" % from 2800 to 3800 cm-1
        [~,edge1] = min(abs(RamanshiftOri - 2600)); % find the data series1 (2600 cm-1);
        [~,edgeA] = min(abs(RamanshiftOri - 2800)); % find the data seriesA (2800 cm-1);
        [~,edgeB] = min(abs(RamanshiftOri - 3800)); % find the data seriesB (3800 cm-1);
        [~,edge2] = min(abs(RamanshiftOri - 4000)); % find the data series2 (4000 cm-1);
        RamanshiftCut = RamanshiftOri(edge1:edge2);
        IntensityCut = IntensityOri(edge1:edge2); % new Raman Intensity;
        FittingCut = RamanshiftOri([edge1:edgeA, edgeB:edge2]);
        FittingInt = IntensityOri([edge1:edgeA, edgeB:edge2]);
        RamanEdge = RamanshiftOri(edgeA:edgeB);
        IntensityEdge = IntensityOri(edgeA:edgeB);
        NumCut = length(IntensityCut);
        NumEdge = length(IntensityEdge);
        % otherwise
        % ...
        % return
end

figure; hold on
FittingP2 = polyfit(FittingCut,FittingInt,2);
polyInt2 = polyval(FittingP2,RamanEdge);
NewIntensity = IntensityEdge - polyInt2; % Intensity from 2800 to 3800 cm-1 after baseline correction.
plot(RamanEdge,NewIntensity);
plot(RamanEdge,polyInt2,'.');
plot(RamanshiftOri,IntensityOri);
xlabel('Wavenumber(nm)','fontsize',15);
% xlabel('Raman shift(cm^{-1})','fontsize',15)
ylabel('Intensity(a.u.)','fontsize',15);
hold off

mint=min(RamanEdge);
maxt=max(RamanEdge);
dt=maxt-mint;

Width1=dt/10;
Width2=dt/5;
Width3=dt/5;
Width4=dt/10;
Width5=dt/10;

Height1=0.01;
Height2=0.05;
Height3=0.05;
Height4=0.01;
Height5=0.01;

switch SpectraType
    case "Wavenumber"
        Position1=633.593;
        Position2=641.149;
        Position3=650.745;
        Position4=656.815;
        Position5=659.587;
    case "RamanShift"
        Position1=3014;
        Position2=3220;
        Position3=3430;
        Position4=3572;
        Position5=3635;
end

Peak1=[Height1 Position1 Width1];
Peak2=[Height2 Position2 Width2];
Peak3=[Height3 Position3 Width3];
Peak4=[Height4 Position4 Width4];
Peak5=[Height5 Position5 Width5];

% Starting point "start" has the form [position1 width1 position2 width2]
start=[Position1 Width1 Position2 Width2 Position3 Width3 Position4 Width4 Position5 Width5];  % automatically calculcated start (not very good)
options = optimset('TolX',1);  % Determines how close the model must fit the data
% figure(2); hold on
figure;
parameter=fminsearch(@(lambda)(fitgauss2animatedError(lambda,RamanEdge',NewIntensity')),start)
HBTI=(parameter(6)+parameter(8)+parameter(10))/(parameter(2)+parameter(4))
toc