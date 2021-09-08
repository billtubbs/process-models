function CSTR5_S_function(block)
%
%   Level-2 MATLAB S-Function for CSTR reactor simulation in Simulink
%   This is based on the MSFUNTMPL_BASIC Template from the documentation
%   found here:
%   https://www.mathworks.com/help/simulink/sfg/writing-level-2-matlab-s-functions.html

%%
%% The setup method is used to set up the basic attributes of the
%% S-function such as ports, parameters, etc. Do not add any other
%% calls to the main body of the function.
%%
setup(block);

%endfunction

%% Function: setup ===================================================
%% Abstract:
%%   Set up the basic characteristics of the S-function block such as:
%%   - Input ports
%%   - Output ports
%%   - Dialog parameters
%%   - Options
%%
%%   Required         : Yes
%%   C MEX counterpart: mdlInitializeSizes
%%
function setup(block)

% Register number of ports
block.NumInputPorts  = 3;
block.NumOutputPorts = 2;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
block.InputPort(1).Dimensions  = 1;
block.InputPort(1).DatatypeID  = 0;  % double
block.InputPort(1).Complexity  = 'Real';
block.InputPort(1).DirectFeedthrough = true;
block.InputPort(2).Dimensions  = 1;
block.InputPort(2).DatatypeID  = 0;  % double
block.InputPort(2).Complexity  = 'Real';
block.InputPort(2).DirectFeedthrough = true;
block.InputPort(3).Dimensions  = 5;
block.InputPort(3).DatatypeID  = 0;  % double
block.InputPort(3).Complexity  = 'Real';
block.InputPort(3).DirectFeedthrough = true;

% Override output port properties
block.OutputPort(1).Dimensions  = 1;
block.OutputPort(1).DatatypeID  = 0; % double
block.OutputPort(1).Complexity  = 'Real';
block.OutputPort(2).Dimensions  = 1;
block.OutputPort(2).DatatypeID  = 0; % double
block.OutputPort(2).Complexity  = 'Real';

% Define number of continuous states
block.NumContStates = 5;  % TODO: Check this is correct place for this

% Register parameters
block.NumDialogPrms     = 0;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
%% The MATLAB S-function uses an internal registry for all
%% block methods. You should register all relevant methods
%% (optional and required) as illustrated below. You may choose
%% any suitable name for the methods and implement these methods
%% as local functions within the same file. See comments
%% provided for each function for more information.
%% -----------------------------------------------------------------

block.RegBlockMethod('InitializeConditions', @InitializeConditions);
block.RegBlockMethod('Outputs', @Outputs);     % Required
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('SetInputPortSamplingMode', @SetInputPortSamplingMode);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup

%%
%% InitializeConditions:
%%   Functionality    : Called at the start of simulation and if it is 
%%                      present in an enabled subsystem configured to reset 
%%                      states, it will be called when the enabled subsystem
%%                      restarts execution to reset the states.
%%   Required         : No
%%   C MEX counterpart: mdlInitializeConditions
%%
function InitializeConditions(block)

% Initial values of states
block.ContStates.Data(1) = 5.63;
block.ContStates.Data(2) = 0.6388;
block.ContStates.Data(3) = 331.25;
block.ContStates.Data(4) = 0.001659;
block.ContStates.Data(5) = 44.26;

%end InitializeConditions

%%
%% Outputs:
%%   Functionality    : Called to generate block outputs in
%%                      simulation step
%%   Required         : Yes
%%   C MEX counterpart: mdlOutputs
%%
function Outputs(block)

x = [block.ContStates.Data(1);
     block.ContStates.Data(2);
     block.ContStates.Data(3);
     block.ContStates.Data(4);
     block.ContStates.Data(5)];

y = cstr5_measurements(x);

block.OutputPort(1).Data = y(1);
block.OutputPort(2).Data = y(2);

%end Outputs

%%
%% Derivatives:
%%   Functionality    : Called to update derivatives of
%%                      continuous states during simulation step
%%   Required         : No
%%   C MEX counterpart: mdlDerivatives
%%
function Derivatives(block)

% Kinetic parameters
params.Z.Tc = 3.8223e10;  % kmol/m^3.h
params.Z.Td = 3.14157e11;  % kmol/m^3.h
params.Z.I = 3.7920e18;  % h^-1
params.Z.P = 1.7700e9;  % kmol/m^3.h
params.Z.fm = 1.0067e15;  % kmol/m^3.h
params.E.Tc = 2.9442e3;  % kJ/kmol
params.E.Td = 2.9442e3;  % kJ/kmol
params.E.I = 1.2877e5;  % kJ/kmol
params.E.P = 1.8283e4;  % kJ/kmol
params.E.fm = 7.4478e4;  % kJ/kmol

% Physical parameters
params.F = 1.0;  % m^3.h^-1
params.V = 0.1;  % m^3
params.rho = 866;  % kg.m^-3
params.Tj = 295;  % K
params.Tin = 350;  % K
params.A = 2.0;  % m^2
params.f_star = 0.58;
params.FI = 0.08;  % m^3.h^-1
params.Fm = 0.92;  % m^3.h^-1
params.Mm = 100.12;  % kg.kmol^-1
params.mDeltaHp = 57800;  % kJ.kmol^-1
params.cp =  2.0;  % kJ.kg^-1
params.U = 800;  % kJ.h^-1.K-1.m^-2
params.R = 8.314;  % kJ.mol^-1.K^-1

% State vector
x = [block.ContStates.Data(1);
     block.ContStates.Data(2);
     block.ContStates.Data(3);
     block.ContStates.Data(4);
     block.ContStates.Data(5)];

% Input vector
p = [block.InputPort(1).Data;
     block.InputPort(2).Data];

% State disturbances
w = block.InputPort(3).Data;

% Continuous-time differential equations
dx = cstr5_dynamics_CT(x, p, w, params);

block.Derivatives.Data(1) = dx(1);
block.Derivatives.Data(2) = dx(2);
block.Derivatives.Data(3) = dx(3);
block.Derivatives.Data(4) = dx(4);
block.Derivatives.Data(5) = dx(5);

%end Derivatives


%% SetInpPortFrameData function

function SetInputPortSamplingMode(block, idx, fd)
  block.InputPort(idx).SamplingMode = fd;
  block.OutputPort(1).SamplingMode  = fd;
  block.OutputPort(2).SamplingMode  = fd;

%end SetInpPortFrameData


%%
%% Terminate:
%%   Functionality    : Called at the end of simulation for cleanup
%%   Required         : Yes
%%   C MEX counterpart: mdlTerminate
%%
function Terminate(block)

%end Terminate

