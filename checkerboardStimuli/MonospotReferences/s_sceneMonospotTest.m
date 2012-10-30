% s_sceneMonospotTest
%
% Excercise ISET, as well as test local routines.
%
% 7/14/11  dhb, gt  Started down this road
% 7/22/11  dhb      Now working for initial monospot.
% 7/25/11  bw       A little more ISET background added

%% Initialize.  Close ISET windows effectively and open ISET again.  
% The command ieMainW will bring up the main window if you want it (probably you don't).
s_initIset;  

%% Create an example scene, using an ISET preset
% We can add this to sceneCreate when you are ready.  Comments within
% sceneMonospot.
scene = sceneMonospot;
vcAddAndSelectObject(scene); sceneWindow;

%% Define an optical system that models the human eye
oi = oiCreate('human');
oi = oiCompute(scene,oi);
vcAddAndSelectObject(oi);
oiWindow;

%% Define a sensor (here, the retina)
sensor = sensorCreate('human');

% We have sensor electrical noise, by default
% Here, we set the noise flag to only photon noise
% 0 is no noise, 1 photon noise, 2 electrical and photon noise
sensor = sensorSet(sensor,'noise flag',1);  % Only photon noise

% Adjust the size of the sensor to match the extent of the scene
sensor = sensorSetSizeToFOV(sensor,sceneGet(scene,'fov'));

% Cautious about exposure time.  Default is auto exposure which produced
% 100 sec of integration time.
expTime = 0.2;
sensor = sensorSet(sensor,'exp time',expTime);   % 200 ms

% Compute sensor image and display
sensor = sensorCompute(sensor,oi);
vcAddAndSelectObject(sensor);
sensorImageWindow;

% Plot a horizontal line about half way through the sensor rows (y)
sz = sensorGet(sensor,'size');
sensorPlotLine(sensor,'h','electrons','space',[1, sz(1)/2])
