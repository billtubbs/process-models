function drawcartpend_bw(y,params)
% Renders a drawing of an inverted pendulum on a cart
% in the current figure.
%
% Arguments:
%   y : (4, 1) double
%     State vector.
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

% positions
y = wr/2 + H/2; % cart vertical position
w1x = x - .9*W/2;
w1y = 0;
w2x = x + .9*W/2 - wr;
w2y = 0;

% Pole end position
px = x + l*sin(th);
py = y + l*cos(th);

% Draw line for cart track
plot([-10 10],[0 0],'w','LineWidth',2)
hold on

% Cart
rectangle( ...
    'Position',[x-W/2,y-H/2,W,H], ...
    'Curvature',.1, ...
    'FaceColor',[1 0.1 0.1], ...
    'EdgeColor',[1 1 1] ...
)

% Wheels
rectangle( ...
    'Position',[w1x,w1y,wr,wr], ...
    'Curvature',[1 1], ...
    'FaceColor',[1 1 1], ...
    'EdgeColor',[1 1 1] ...
)
rectangle( ...
    'Position',[w2x,w2y,wr,wr], ...
    'Curvature',[1 1], ...
    'FaceColor',[1 1 1], ...
    'EdgeColor',[1 1 1] ...
)

% Pole
plot([x px],[y py],'w','LineWidth',2)
rectangle( ...
    'Position',[px-mr/2,py-mr/2,mr,mr], ...
    'Curvature',[1 1], ...
    'FaceColor',[.3 0.3 1], ...
    'EdgeColor',[1 1 1] ...
)

% This ensures the x and y axis have same scaling
axis equal

% set(gca,'YTick',[])
% set(gca,'XTick',[])
xlim([-5 5]);
ylim([-2 2.5]);
axes_color = [0.5 0.5 0.5];
set(gca,'Color','k','XColor',axes_color,'YColor',axes_color)
set(gcf,'Position',[50 500 800 400])
set(gcf,'Color','k')
set(gcf,'InvertHardcopy','off')

% box off
drawnow
hold off