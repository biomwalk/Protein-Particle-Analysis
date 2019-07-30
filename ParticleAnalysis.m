%% Particle Analysis 

%%% Code for Analyzing Protein, Silicon Oil, and Free Fatty Acid
%%% Particles in solutions. Reads image data from flow cytometer
%%% selects certain columns for display in various ways.
%% Data Imports

[ProtNum, Prottxt, Protraw] = xlsread('181024 Prot 100 data_export.xlsx','A1:AI10000');   % Protein Base Data
[FFANum, FFAtxt, FFAraw] = xlsread('181024 Myr 100 data_export.xlsx','A1:AI10000');       % FFA Base Data
[OilNum, Oiltxt, Oilraw] = xlsread('181024 Oil 100 data_export.xlsx','A1:AI10000');       % Silicon Oil Base Data
[ProtOil5050Num, ProtOil5050txt, ProtOil5050raw] = xlsread('181024 Prot Oil 5050 data_export.xlsx','A1:AI10000');
[FFAProt5050Num, FFAProt5050txt, FFAProt5050raw] = xlsread('181024 Myr Prot 5050 data_export.xlsx','A1:AI10000');
[FFAOil5050Num, FFAOil5050txt, FFAOil5050raw] = xlsread('181024 Myr Oil 5050 data_export.xlsx','A1:AI10000');
[ProtOil9010Num, ProtOil9010txt, ProtOil9010raw] = xlsread('181024 Prot Oil 9010 data_export.xlsx','A1:AI10000');
[ProtFFA9010Num, ProtFFA9010txt, ProtFFA9010raw] = xlsread('181024 Prot Myr 9010 data_export.xlsx','A1:AI10000');
[FFAOil9010Num, FFAOil9010txt, FFAOil9010raw] = xlsread('181024 Myr Oil 9010 data_export.xlsx','A1:AI10000');
[OilProt9010Num, OilProt9010txt, OilProt9010raw] = xlsread('181024 Oil Prot 9010 data_export.xlsx','A1:AI10000');
[FFAProt9010Num, FFAProt9010txt, FFAProt9010raw] = xlsread('181024 Myr Prot 9010 data_export.xlsx','A1:AI10000');
[OilFFA9010Num, OilFFA9010txt, OilFFA9010raw] = xlsread('181024 Oil Myr 9010 data_export.xlsx','A1:AI10000');
[Mix33Num, Mix33txt, Mix33raw] = xlsread('181024 Prot Myr Oil 333333 data_export.xlsx','A1:AI10000');


%% Analysis Tools

SampleNum = FFAOil9010Num;        % Select Dataset for analysis

EndParticle = 500;  % Number of Particles to be analyzed
StartParticle = 1;     % Starting value from excel sheet to be analyzed
NumParticles = EndParticle - StartParticle + 1; % Total Number of Particles to be analyzed
Columns = [4, 7, 12, 18, 21, 25];  % Vals Corresponding to morphologies of interest


%Sorts data files by ESD (descending)
w = sortrows(FFANum,-18);
x = sortrows(ProtNum,-18);
y = sortrows(OilNum,-18);
z = sortrows(SampleNum,-18);


%Creates values matrix for grid plot
W = w(StartParticle:EndParticle, Columns);
X = x(StartParticle:EndParticle, Columns);
Y = y(StartParticle:EndParticle, Columns);
Z = z(StartParticle:EndParticle, Columns);
Data = [W;X;Y;Z];

%Creates values matrix for parallel coords plot

Columns2 = [1:5 6:10 11:15 16:20 21:25 30:35];
A = FFANum(StartParticle:EndParticle, Columns2);
B = ProtNum(StartParticle:EndParticle, Columns2);
C = OilNum(StartParticle:EndParticle, Columns2);
D = SampleNum(StartParticle:EndParticle, Columns2);
Data2 = [A;B;C;D];

yNames = ['AR ', 'CirHu ', 'ESD ', 'GAR ', 'Int ', 'SI'];

id = [];   % Initializes Identity Matrix
for i = 0:3
    id = [id; zeros(NumParticles, 1) + i];  % Assigns Identity values to identity matrix
end

% Implements grid plot matrix.
figure (1)
gplotmatrix(Data,[],id,['b' 'c' 'g' 'r' 'm' 'y'],[],[],false);
text([.08 .24 .43 .66 .83], repmat(-.1,1,5), yNames, 'FontSize',7);
text(repmat(-.12,1,5), [.86 .62 .41 .25 .02], yNames, 'FontSize',7, 'Rotation',90);
title('Plot Matrix of AR, CircHu, ESD, GAR, Int, SI');

labels = FFAtxt(1,Columns2);

% Implements Parallel coordinates plot using Identity matrix values as
% unique conditions for comparison
figure (2)
mixVmyr = ismember(id,[0 1 2 3]);
parallelcoords(Data2(mixVmyr,:), 'group',id(mixVmyr), ...
               'labels',labels,'standardize','on', 'quantile', .25)
title('Coordinate plot of third few flowcam morphologies');

legend('FFA','Prot','Oil','Sample')

    