%
% COAG_DRIVER is the main command-line driver for calculating particle size
% spectra. 
%
% USAGE:
%   Type coag_driver at the command line and follow the instructions
%
% HISTORY:
%   23-04-09: First Cut
%
% Adrian Burd, University of Georgia, 2009
%
%

%% Set up and get user options
%

p = SetUpCoag;

%% Calculate the sectionally integrated coagulation kernels
%

disp('Calculating kernels')

%b_brown = CalcBetas(p);
%b_brown.b1 = b_brown.b1*p.conBr*p.day_to_sec;
%b_brown.b2 = b_brown.b2*p.conBr*p.day_to_sec;
%b_brown.b3 = b_brown.b3*p.conBr*p.day_to_sec;
%b_brown.b4 = b_brown.b4*p.conBr*p.day_to_sec;
%b_brown.b5 = b_brown.b5*p.conBr*p.day_to_sec;

p.kernel='KernelFracSh';
b_shear = CalcBetas(p);
b_shear.b1 = b_shear.b1*p.gamma*p.day_to_sec;
b_shear.b2 = b_shear.b2*p.gamma*p.day_to_sec;
b_shear.b3 = b_shear.b3*p.gamma*p.day_to_sec;
b_shear.b4 = b_shear.b4*p.gamma*p.day_to_sec;
b_shear.b5 = b_shear.b5*p.gamma*p.day_to_sec;
b_shear.b25 = b_shear.b25*p.gamma*p.day_to_sec;

p.kernel='KernelFracDS';
b_ds    = CalcBetas(p);
b_ds.b1 = b_ds.b1*p.setcon*p.day_to_sec;
b_ds.b2 = b_ds.b2*p.setcon*p.day_to_sec;
b_ds.b3 = b_ds.b3*p.setcon*p.day_to_sec;
b_ds.b4 = b_ds.b4*p.setcon*p.day_to_sec;
b_ds.b5 = b_ds.b5*p.setcon*p.day_to_sec;
b_ds.b25 = b_ds.b25*p.setcon*p.day_to_sec;

% Pack up the betas and store them in a new structure that will get passed
% to the derivative and jacobian calculation routines

% p2.b1 = b_brown.b1 + b_shear.b1 + b_ds.b1;
% p2.b2 = b_brown.b2 + b_shear.b2 + b_ds.b2;
% p2.b3 = b_brown.b3 + b_shear.b3 + b_ds.b3;
% p2.b4 = b_brown.b4 + b_shear.b4 + b_ds.b4;
% p2.b5 = b_brown.b5 + b_shear.b5 + b_ds.b5;


p2.b1 =  b_shear.b1 + b_ds.b1;
p2.b2 =  b_shear.b2 + b_ds.b2;
p2.b3 =  b_shear.b3 + b_ds.b3;
p2.b4 =  b_shear.b4 + b_ds.b4;
p2.b5 =  b_shear.b5 + b_ds.b5;

p2.b25 = p2.b2 - p2.b3 - p2.b4 - p2.b5;

%% Calculate the linear terms in the population balance equation
%

growth    = CalcGrowth(p);
sink_loss = CalcSinkingLoss(p);

p2.linear   = growth - sink_loss;

%% Initial Size Spectrum
%  Caclculate the initial size spectrum used for both estimation of the
%  steady state or evolving solution - put both in later versions

spec_init = CalcInitialSpec(p, p2);

%% Integrate
%  Set up for integrating over time

disp('Solving ODEs')

calcomp = 1:p.n_sections;
abs_tol = 1.0e-18;        % Absolute Tolerance baseline
rel_tol = 3.0e-14;        % Relative tolerance

at = (abs_tol * 1.5 .^ (-(calcomp-1)));

t_span = p.t_init : p.delta_t : p.t_final - 1;

ode_options = odeset('RelTol', rel_tol, 'Refine', 0, 'AbsTol', at, 'Jacobian', @CalcCoagJac);
%ode_options = odeset('Jacobian', @CalcCoagJac);

[t_out, y] = ode15s(@CalcCoagDeriv, t_span, spec_init, ode_options, p2); 

%% Output
%

outflag = CoagOutput(p, p2, t_out, y);
