function draw_cartpole(y,params)
% Renders a drawing of an inverted pendulum on a cart
% in the current figure.
%
% Arguments:
%   y : (4, 1) double
%     State vector.
%       y(1) : cart horizontal position
%       y(2) : cart velocity (not used)
%       y(3) : pendulum angle, clockwise from vertical up
%       y(4) : pendulum angular velocity (not used)
%   params : struct
%       Parameter values which are used to determine
%       the rendered dimensions of the cart-pole:
%           params.l : pendulum length
%           params.mc : cart mass
%           params.mp : pole mass.
%

% This code is based on a similar script by Steve Brunton
% from his Control Bootcamp course.

% Current state
x = y(1);    % cart horizontal position
th = y(3);   % pendulum angle, clockwise from 
             % vertical up position

% Parameters
mp = params.mp;
mc = params.mc;
l = params.l;

% Dimensions
W = 1*sqrt(mc/5);  % cart width
H = .5*sqrt(mc/5); % cart height
wr = .2; % wheel radius
mr = .3*sqrt(mp); % mass radius

% Positions
y = wr/2 + H/2; % cart vertical position
w1x = x - .9*W/2;
w1y = 0;
w2x = x + .9*W/2 - wr;
w2y = 0;

% Pole end position
px = x + l*sin(th);
py = y + l*cos(th);

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

% Pole
plot([x px],[y py],'k','LineWidth',2)
rectangle( ...
    'Position',[px-mr/2,py-mr/2,mr,mr], ...
    'Curvature',[1 1], ...
    'FaceColor',pole_face_color, ...
    'EdgeColor',edge_color ...
)

% This ensures the x and y axis have same scaling
axis equal

% set(gca,'YTick',[])
% set(gca,'XTick',[])
xlim([-5 5]);
ylim([-2 2.5]);
set(gca,'XColor',axes_color,'YColor',axes_color)
set(gcf,'Position',[50 500 800 400])

% box off
drawnow
hold off