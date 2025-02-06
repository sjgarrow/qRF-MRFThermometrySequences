addpath(genpath('C:\Users\garro\Documents\ResearchGeneral\pulseq-master'))
addpath(genpath('C:\Users\garro\Documents\SGarrowMatlabFiles\vds'))

%% vds spiral code
% Define FOV and resolution
fov = 0.256; % m
Nx = 128; Ny = Nx;
% deltak = 1 / fov;
% kWidth = Nx * deltak;
thickness = 5e-3; % EX slice thickness - meters
% TE = 1.8e-3;
TR = 10e-3;
alpha = 10; %flip angle, degrees

%%get the system definition - 'sys' variable
MaxGrad = 22; % REDFLAG set for Vida
MaxSlew = 100; % was 130
sys = mr.opts('MaxGrad', MaxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', MaxSlew, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6,'B0',3);

% system parameters, with safety factors
%	smax = maximum slew rate G/cm/s
%	gmax = maximum gradient G/cm (limited by Gmax or FOV)
%	T = sampling period (s) for gradient AND acquisition.
%	N = number of interleaves.
%	Fcoeff = FOV coefficients with respect to r - see above.
%	rmax= value of k-space radius at which to stop (cm^-1).
%		rmax = 1/(2*resolution);

% Convert max slew rate from T/m/s to G/cm/s -- mr opts converts max slew
% rate to be in Hz/m, so I'm not using the sys var
smax = MaxSlew * 10000 / 100 * 0.9; % 0.9 safety factor 
% Convert max grad from mT/m to G/cm
gmax = MaxGrad / 10;
% T = .000004;	 % Seconds
T = sys.gradRasterTime / 10; % 1e-6 s
nInterleaves = 12;		 % Interleaves
% Fcoeff = [25.6 25.6];	% FOV decreases linearly from 24 to 12cm.
Fcoeff = [fov*100];	% FOV is 256e-3m, * 100 would make it cm
% res = 1;
res = (fov / Nx) * 100 ; % cm --> 2 mm x 2 mm resolution
rmax = 1/(2 * res);		% cm^(-1) --> corresponds to 2mm resolution.

disp('Calculating Gradient');
[k, g, ~, time, r, theta] = vds(smax,gmax,T,nInterleaves,Fcoeff,rmax);
% Tramp = 0.0005; 
Tramp = 0.00025; 
gx_ramp = real(g(end)) * (Tramp / T:-1:1) / (Tramp / T);
gy_ramp = imag(g(end)) * (Tramp / T:-1:1) / (Tramp / T);
g = [g (gx_ramp + 1i * gy_ramp)];


% WAG: I get the same waveforms with this (good!): 
% [km, gm, ~, time] = vdsmex(nInterleaves, Fcoeff, res * 10, gmax, smax, T, 10000000);
% max(abs([real(k)'; imag(k)'] - col(km(1:end-1,:))))
% max(abs([real(g)'; imag(g)'] - col(gm(1:end-1,:))))
ktrunc = k(1:10:end);
kspiral = [real(ktrunc); imag(ktrunc)];
gspiral = g(1:10:end);
gspiral = [real(gspiral); imag(gspiral)];
% TO DO - rotate k space for 12 trajectories
kSpirals = zeros(2,size(kspiral,2),nInterleaves);
% gSpirals = zeros(2,size(kspiral,2),nInterleaves);

%%
	% OUTPUTS:
%	--------
%	k = k-space trajectory (kx+iky) in cm-1.
%	g = gradient waveform (Gx+iGy) in G/cm.
%	s = derivative of g (Sx+iSy) in G/cm/s.
%	time = time points corresponding to above (s).
%	r = k-space radius vs time (used to design spiral)
%	theta = atan2(ky,kx) = k-space angle vs time.

[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',0.5e-3,...
    'SliceThickness',thickness,'timeBwProduct',2,'system',sys); %,'use','excitation');
% [rf, gz] = mr.makeGaussPulse(alpha*pi/180,'system',sys,...
%     'timeBwProduct',2,'SliceThickness',thickness,'Duration',1.2e-3,...
%     'maxGrad', sys.maxGrad, 'maxSlew', sys.maxSlew);
gzReph = mr.makeTrapezoid('z', sys, 'Area', -gz.area / 2);


gxSpiral = mr.convert(squeeze(gspiral(1, :)) * 10, 'mT/m', 'Hz/m'); 
gySpiral = mr.convert(squeeze(gspiral(2, :)) * 10, 'mT/m', 'Hz/m');

gx = mr.makeArbitraryGrad('x', gxSpiral, 'first', 0, ...
    'Delay', 7.5000e-04);
gy = mr.makeArbitraryGrad('y', gySpiral, 'first', 0, ...
    'Delay', 7.5000e-04);

% make a rewinder for whichever channel needs the most rewinding
% copy that rewinder and scale it for the other channel
if abs(gx.area) > abs(gy.area)
    gx_rew = mr.makeTrapezoid('x', sys, 'Area', -gx.area); 
    gy_rew = gx_rew; 
    gy_rew.area = -gy.area;
    gy_rew.flatArea = gx_rew.flatArea * (-gy.area) / (-gx.area);
    gy_rew.amplitude = gx_rew.amplitude * (-gy.area) / (-gx.area);
    gy_rew.channel = 'y';
else
    gy_rew = mr.makeTrapezoid('y', sys, 'Area', -gy.area); 
    gx_rew = gy_rew; 
    gx_rew.area = -gx.area;
    gx_rew.flatArea = gy_rew.flatArea * (-gx.area) / (-gy.area);
    gx_rew.amplitude = gy_rew.amplitude * (-gx.area) / (-gy.area);
    gx_rew.channel = 'x';
end

% Define other gradients and ADC events - round down to nearest multiple
% of 1000 samples
adc = mr.makeAdc(floor(mr.calcDuration(gx) / 2e-6 / 1000) * 1000, ...
    'Dwell', 2e-6, 'Delay', 7.5000e-04); 


% Create a new sequence object
seq = mr.Sequence(sys);

delayTR = ceil((TR - max([mr.calcDuration(gx), mr.calcDuration(gy), mr.calcDuration(adc)])  ...
                   - max([mr.calcDuration(gx_rew), mr.calcDuration(gy_rew)]) ...
                   - mr.calcDuration(gz) - 2 * mr.calcDuration(gzReph)) / seq.gradRasterTime / 2) * seq.gradRasterTime;

adc.delay = delayTR;
gx.delay = delayTR;
gy.delay = delayTR;
gx.waveform = gx.waveform.';
gy.waveform = gy.waveform.';

delayTR = mr.makeDelay(delayTR);
delayTRpost = delayTR;
delayTRpost.delay = delayTRpost.delay - seq.gradRasterTime + 7.5e-04;

rf_inc = 4.2;


% add delay TE back in

% Define sequence blocks
for s = 0 : 4
    rf_phase = wrapTo180(rf_inc * s.^2);
    rf_phase
    rf.phaseOffset = rf_phase / 180 * pi;
    adc.phaseOffset = rf_phase / 180 * pi;
    seq.addBlock(rf, gz);
    
    ind = mod(s, nInterleaves);
    theta = ind * 2 * pi / nInterleaves;
    seq.addBlock(gzReph);

    % seq.addBlock(delayTR);
    seq.addBlock(mr.rotate('z', theta, gx, gy, adc));
    seq.addBlock(mr.rotate('z', theta, gx_rew, gy_rew));

    seq.addBlock(delayTRpost);
    seq.addBlock(gzReph)

end

[ok, error_report] = seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

seq.setDefinition('FOV', [fov fov thickness]);
seq.setDefinition('Name', 'gre');

% seq.write('Therm12ShotSpiral.seq')       % Write to pulseq file

%seq.install('siemens');

seq.plot('timeRange', [0 2]*TR);


%%
% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kx-axis
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
%%
rep = seq.testReport;
fprintf([rep{:}]);