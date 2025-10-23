%{
Increased the input amplitude for stimulus encoding to simulate stronger
mem. encoding
%}

function [ang_cu_mat, ang_uc_mat, metrics] = simulation_inc_Memenco(ori_deg, dist_params, rng_seed)

% Stimulus angles in degrees b/w -90 and +90
stim_ori_cu = deg2rad(ori_deg(1));  stim_ori_uc = deg2rad(ori_deg(2)); dist_ori = deg2rad(ori_deg(3));
dist_on_ms = dist_params(1); 
dist_side = dist_params(2);

if nargin < 3
    rng_seed = 'default';
end

rng(rng_seed); % For reproducibility

%%%%% PARAMETERS

N=512; 		% number of "neurons" in the rate model
npop=8; 	% number of cues presented

totalTime=3050;	% total time of the simulation in ms
dt=2; 		% integration step in ms

tauE=20;	% time constant of rate equation for excitatory neurons
tauI=10;	% time constant of rate equation for inhibitory neurons

kappa=1.5;	% parameter defining concentration of e-to-e connectivity
GEE=6;	 	% strength of excitation to excitatory neurons
GEI=4;		% strength of excitation to inhibitory neurons
GIE=3.4;	% strength of inhibition to excitatory neurons
GII=0.85;	% strength of inhibition to inhibitory neurons

I0E=0.2;	% external bias current to excitatory neurons
I0I=0.5;	% external bias current to inhibitory neurons

sigE=3;	%1	% standard deviation of additive noise in rate equation of e-cells
sigI=5;	%3	% standard deviation of additive noise in rate equation of i-cells
sigE1=3;	% standard deviation of additive noise in rate equation of e-cells PFC
sigI1=5;    % standard deviation of additive noise in rate equation of i-cells PFC

stimon = 500;	% time when external stimulus is applied in ms
stimoff = 1000;	% time when external stimulus ceases in ms
time_cueON_ms = 1100; % time when the retrocue appears
stim = 200;  	% strength of external stimulus

mutual_inhibition = 0.02;%0.02; 
synaptic_plasticity = 0; %0.075;
GEE1 = 3; % strength of bottom-up excitatory connections from VC to PFC 
GEE2 = 0.1; % strength of top-down excitatory connections from PFC to VC
%%%%% PRELIMINARY CALCULATIONS

rE_cu=zeros(N,1); rE_uc=zeros(N,1); 
rI_cu=zeros(N,1); rI_uc=zeros(N,1);
rE_cu1=zeros(N,1); rE_uc1=zeros(N,1); 
rI_cu1=zeros(N,1); rI_uc1=zeros(N,1);
nsteps=floor(totalTime/dt);
delayPop=zeros(N,1);

% E-to-E connectivity
theta = [0:N-1]/N*pi;
v = exp(kappa*cos(theta*2));
v = v/sum(v);
WE = gallery('circul',v);

% stimulus parameters
theta=theta-pi/2;
v_cu = exp(1.8*kappa*cos(2*(theta - stim_ori_cu)));  v_uc = exp(1.8*kappa*cos(2*(theta - stim_ori_uc)));      v_d = exp(1.5*kappa*cos(2*(theta - dist_ori)));
v_cu = v_cu/sum(v_cu);                             v_uc = v_uc/sum(v_uc);                                 v_d = v_d/sum(v_d);
stimulus_cu = 5*stim*v_cu';                          stimulus_uc = 5*stim*v_uc';                              stimulus_dist = 5*stim*v_d'; 
stimon = floor(stimon/dt);                                                                                dist_on = floor(dist_on_ms/dt);    
stimoff = floor(stimoff/dt);                                                                              dist_off = floor((dist_on_ms+100)/dt);    
time_cueON = floor(time_cueON_ms/dt);

% Synaptic plasticity
v_cu = v_cu - min(v_cu);                              v_uc = v_uc - min(v_uc);
v_cu = v_cu/max(v_cu);                                v_uc = v_uc/max(v_uc);
WE_cu1 = WE.*(1 + synaptic_plasticity*(v_cu+v_cu'));  WE_uc1 = WE.*(1 + synaptic_plasticity*(v_uc+v_uc'));
WE_cu = WE;                                           WE_uc = WE;

% input-output function for all cells, as used previously (Brunel, Cereb Cortex 13:1151, 2003)
f = inline('x.*x.*(x>0).*(x<1)+sqrt(4*x-3).*(x>=1)');

% population vector decoder given the rates r for neurons with selectivity th
decode = inline('atan2(sum(r.*sin(th)),sum(r.*cos(th)))','r','th');

nbl=floor(N/npop);

%%%% SIMULATION LOOP

for i=1:nsteps,

  % additive noise for each population
  noiseE_cu = sigE*randn(N,1);   noiseE_uc = sigE*randn(N,1);
  noiseI_cu = sigI*randn(N,1);   noiseI_uc = sigI*randn(N,1);
  noiseE_cu1 = sigE1*randn(N,1);   noiseE_uc1 = sigE1*randn(N,1);
  noiseI_cu1 = sigI1*randn(N,1);   noiseI_uc1 = sigI1*randn(N,1);
  
  % Reduce the noise on the cued side after cued onset 
  if (i>time_cueON) % After cue onset
        noiseE_cu = 0.4*sigE*randn(N,1); 
        noiseI_cu = 0.4*sigI*randn(N,1);
        noiseE_cu1 = 0.4*sigE1*randn(N,1); 
        noiseI_cu1 = 0.4*sigI1*randn(N,1);
  end
  
  % current input to each population
  IE_cu=GEE*WE_cu*rE_cu+(I0E-GIE*mean(rI_cu))*ones(N,1);    IE_uc=GEE*WE_uc*rE_uc+(I0E-GIE*mean(rI_uc))*ones(N,1);
  II_cu=(GEI*mean(rE_cu)-GII*mean(rI_cu)+I0I)*ones(N,1);    II_uc=(GEI*mean(rE_uc)-GII*mean(rI_uc)+I0I)*ones(N,1);
  IE_cu1=GEE*WE_cu*rE_cu1+(I0E-GIE*mean(rI_cu1))*ones(N,1); IE_uc1=GEE*WE_uc*rE_uc1+(I0E-GIE*mean(rI_uc1))*ones(N,1);
  II_cu1=(GEI*mean(rE_cu1)-GII*mean(rI_cu1)+I0I)*ones(N,1); II_uc1=(GEI*mean(rE_uc1)-GII*mean(rI_uc1)+I0I)*ones(N,1);
  
  % external task-dependent inputs
  if i>stimon & i<stimoff, 
    IE_cu=IE_cu + stimulus_cu;  IE_uc=IE_uc + stimulus_uc;
    WE_cu = WE_cu1;             WE_uc = WE_uc1; % Synaptic plasticity -- only once
  
    % Bottom-up inputs from VC to PFC
    IE_cu1=IE_cu1 + GEE1*rE_cu;  IE_uc1=IE_uc1 + GEE1*rE_uc; % Bottom-up connection from VC to PFC
  end
  
  % Distractor
  if (i>dist_on & i<dist_off) & (dist_side ~=3),
      if dist_side ==1
          IE_cu=IE_cu + stimulus_dist;
          IE_uc = IE_uc - 4*8*mutual_inhibition*rE_cu; % Cross-hemispheric inhibition
      elseif dist_side ==2
          IE_uc=IE_uc + stimulus_dist;
          IE_cu = IE_cu - 4*8*mutual_inhibition*rE_uc; % Cross-hemispheric inhibition
      end
  end
  
  % Delay period 
  if (i>stimoff)
      
      % Top-down inputs from the PFC to VC
      if (i>time_cueON) % After cue onset
          IE_cu = IE_cu + 1.2*GEE2*rE_cu1;  IE_uc = IE_uc + GEE2*rE_uc1;
      else % Before cue onset
          IE_cu = IE_cu + GEE2*rE_cu1;      IE_uc = IE_uc + GEE2*rE_uc1;
      end
    
    % Inhibition across hemisphere
    IE_cu = IE_cu - mutual_inhibition*rE_uc;  IE_uc = IE_uc - mutual_inhibition*rE_cu;
  end
  
  % integration with time-step dt: Newton method
  rE_cu = rE_cu + (f(IE_cu) - rE_cu + noiseE_cu)*dt/tauE;          rE_uc = rE_uc + (f(IE_uc) - rE_uc + noiseE_uc)*dt/tauE;
  rI_cu = rI_cu + (f(II_cu) - rI_cu + noiseI_cu)*dt/tauI;          rI_uc = rI_uc + (f(II_uc) - rI_uc + noiseI_uc)*dt/tauI;
  rE_cu1 = rE_cu1 + (f(IE_cu1) - rE_cu1 + noiseE_cu1)*dt/tauE;     rE_uc1 = rE_uc1 + (f(IE_uc1) - rE_uc1 + noiseE_uc1)*dt/tauE;
  rI_cu1 = rI_cu1 + (f(II_cu1) - rI_cu1 + noiseI_cu1)*dt/tauI;     rI_uc1 = rI_uc1 + (f(II_uc1) - rI_uc1 + noiseI_uc1)*dt/tauI;

  metrics.IE_cu_mat(i,:) = IE_cu;     metrics.IE_uc_mat(i,:) = IE_uc;
  metrics.II_cu_mat(i,:) = II_cu;     metrics.II_uc_mat(i,:) = II_uc;
  
  metrics.rE_cu_mat(i,:) = rE_cu;     metrics.rE_uc_mat(i,:) = rE_uc;
  metrics.rI_cu_mat(i,:) = rI_cu;     metrics.rI_uc_mat(i,:) = rI_uc;
  metrics.rE_cu_mat1(i,:) = rE_cu1;   metrics.rE_uc_mat1(i,:) = rE_uc1;
  metrics.rI_cu_mat1(i,:) = rI_cu1;   metrics.rI_uc_mat1(i,:) = rI_uc1;
  
  % get decoded angle from network activity  
  ang_cu=decode(rE_cu,2*theta')/2;      ang_uc=decode(rE_uc,2*theta')/2;
  ang_cu1=decode(rE_cu1,2*theta')/2;    ang_uc1=decode(rE_uc1,2*theta')/2;
  
  ang_cu_mat(i,:) = [rad2deg(ang_cu), rad2deg(ang_cu1)]; ang_uc_mat(i,:) = [rad2deg(ang_uc), rad2deg(ang_uc1)]; 

end

end