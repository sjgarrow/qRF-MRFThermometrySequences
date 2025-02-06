addpath(genpath('C:\Users\garro\Documents\ResearchGeneral\pulseq-master'))

% set system limits
%%get the system definition - 'sys' variable
MaxGrad = 22; % REDFLAG set for Vida
MaxSlew = 100; % was 130
sys = mr.opts('MaxGrad', MaxGrad, 'GradUnit', 'mT/m', ...
    'MaxSlew', MaxSlew, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, ...
    'adcDeadTime', 10e-6,'B0',3);

seq=mr.Sequence(sys);           % Create a new sequence object
fov=256e-3; Nx=128; Ny=128;     % Define FOV and resolution
alpha = 19;                       % flip angle
sliceThickness = 5e-3;            % slice
TR = 17e-3;                       % repetition time TR
TE = 12e-3;                        % echo time TE

%%
% more in-depth parameters
rfSpoilingInc=117;              % RF spoiling increment
roDuration= 3.2e-3;              % ADC duration

% Create fat-sat pulse
% (in Siemens interpreter from January 2019 duration is limited to 8.192 ms, and although product EPI uses 10.24 ms, 8 ms seems to be sufficient)
% B0=2.89; % 1.5 2.89 3.0
% sat_ppm=-3.45;
% sat_freq=sat_ppm*1e-6*B0*lims.gamma;
% rf_fs = mr.makeGaussPulse(110*pi/180,'system',lims,'Duration',8e-3,...
%     'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
% gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

% Create alpha-degree slice selection pulse and gradient
% [rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',3e-3,...
%     'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);

[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',0.5e-3,...
    'SliceThickness',sliceThickness,'timeBwProduct',2,'system',sys); %,'use','excitation');

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys);
% adc = mr.makeAdc(floor(mr.calcDuration(gx) / 2e-6 / 1000) * 1000, ...
%     'Dwell', 2e-6,'Delay',1e-3); 
% adc = mr.makeAdc(floor(mr.calcDuration(gx) * 5 / 1000) * 1000, ...
%     'Dwell', gx.flatTime, 'Delay', gx.riseTime); 
% adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys, ...
%    'Dwell',Nx * sys.adcRasterTime);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
% adc = mr.makeAdc(1000,'Delay',gx.riseTime,'system',sys, ...
%    'Dwell',Nx * sys.adcRasterTime);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',1e-3,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',1e-3,'system',sys);
phaseAreas = ((0:Ny-1)-Ny/2)*deltak;
gyPre = mr.makeTrapezoid('y','Area',max(abs(phaseAreas)),'Duration',mr.calcDuration(gxPre),'system',sys);
peScales=phaseAreas/gyPre.area;


% gradient spoiling
gxSpoil=mr.makeTrapezoid('x','Area',2*Nx*deltak,'system',sys);
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys);

% Calculate timing
delayTE=ceil((TE - mr.calcDuration(gxPre) - gz.fallTime - gz.flatTime/2 ...
    - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
delayTR=ceil((TR - mr.calcDuration(gz) - mr.calcDuration(gxPre) ...
    - mr.calcDuration(gx) - delayTE)/seq.gradRasterTime)*seq.gradRasterTime;
assert(all(delayTE>=0));
assert(all(delayTR>=mr.calcDuration(gxSpoil,gzSpoil)));

rf_phase=0;
rf_inc=0;
% rf_phasevec = zeros(1,Ny);
% Loop over phase encodes and define sequence blocks
nDyn = 55; % number of dynamic images to get
for ii = 1:nDyn
    for i=1:Ny
        for c=1:length(TE)
            %seq.addBlock(rf_fs,gz_fs); % fat-sat
            % rf_inc=mod(rfSpoilingInc, 360.0);
            rf_phase=mod(0.5 * 117 * ((i-1)^2 + (i-1) +1), 360.0);

            % rf_phasevec(i) = rf_phase;
            rf.phaseOffset=rf_phase/180*pi;
            adc.phaseOffset=rf_phase/180*pi;
            %
            seq.addBlock(rf,gz);
            seq.addBlock(gxPre,mr.scaleGrad(gyPre,peScales(i)),gzReph);
            seq.addBlock(mr.makeDelay(delayTE(c)));
            seq.addBlock(gx,adc);
            %gyPre.amplitude=-gyPre.amplitude;
            seq.addBlock(mr.makeDelay(delayTR(c)),gxSpoil,mr.scaleGrad(gyPre,-peScales(i)),gzSpoil)
        end
    end
end
% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end
% Timing check passed successfully
% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'gre');

% seq.write('GREtherm19FAQuadSpoiling.seq')       % Write to pulseq file

%seq.install('siemens');
% plot sequence and k-space diagrams
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
% Current plot held
% Current plot held
% 