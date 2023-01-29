function draw_crane(y,params)
% draw_crane(y,params)
% Renders a drawing of the gantry-crane system
% in the current figure.
%
% Arguments:
%   y : (4, 1) double
%     State vector.
%       y(1) : cart horizontal position
%       y(2) : cart velocity (not used)
%       y(3) : pendulum angle, anti-clockwise from vertical down
%       y(4) : pendulum angular velocity (not used)
%   params : struct
%       Parameter values which are used to determine
%       the rendered dimensions of the cart-pole:
%           params.L : cable length
%           params.mc : cart mass
%           params.mp : load mass.
%

% This code is based on a similar script for a 
% cart-pole system by Steve Brunton from his Control 
% Bootcamp course.

% Current state
x = y(1);    % cart horizontal position
th = y(3);   % cable angle, anti-clockwise from 
             % vertical down position

% Parameters
m_c = params.m_c;
L = params.L;
r = params.r;

% Dimensions
W = 1*sqrt(m_c/5);  % cart width
H = .5*sqrt(m_c/5); % cart height
wr = .2; % wheel radius

% Positions
y = wr/2 + H/2; % cart vertical position
w1x = x - .9*W/2;
w1y = 0;
w2x = x + .9*W/2 - wr;
w2y = 0;

% Pole end position
px = x + L*sin(th);
py = y + L*cos(th);

% Colours
edge_color = [0 0 0];
axes_color = [0.6 0.6 0.6];
cart_face_color = [1 0.1 0.1];
pole_face_color = [.3 0.3 1];

% Draw line for cart track
plot([-10 10],[0 0],'Color',edge_color,'LineWidth',2)
hold on

% Cart
rectangle( ...
    'Position',[x-W/2,y-H/2,W,H], ...
    'Curvature',.1, ...
    'FaceColor',cart_face_color, ...
    'EdgeColor',edge_color ...
)
% Pivot
rectangle( ...
    'Position',[x-0.05,y-0.05,0.1,0.1], ...
    'Curvature',[1 1], ...
    'FaceColor',edge_color ...
)

% Wheels
rectangle( ...
    'Position',[w1x,w1y,wr,wr], ...
    'Curvature',[1 1], ...
    'FaceColor',edge_color, ...
    'EdgeColor',edge_color ...
)
rectangle( ...
    'Position',[w2x,w2y,wr,wr], ...
    'Curvature',[1 1], ...
    'FaceColor',edge_color, ...
    'EdgeColor',edge_color ...
)

% Load
rectangle( ...
    'Position',[px-r/2,py-r/2,r,r], ...
    'FaceColor',pole_face_color, ...
    'EdgeColor',edge_color ...
)
% Pole
plot([x px],[y py],'k','LineWidth',2)
% Hook
rectangle( ...
    'Position',[px-0.05,py-0.05,0.1,0.1], ...
    'Curvature',[1 1], ...
    'FaceColor',edge_color ...
)

% This ensures the x and y axis have same scaling
axis equal

% set(gca,'YTick',[])
% set(gca,'XTick',[])
xlim([-5 5]);
ylim([-5.5 1]);
set(gca,'XColor',axes_color,'YColor',axes_color)
set(gcf,'Position',[50 500 800 600])

% box off
drawnow
hold off