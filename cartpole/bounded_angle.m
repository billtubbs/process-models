function theta = bounded_angle(theta)
% If theta is out of bounds (-pi, pi) bring back
% within these bounds.
%

    if (theta > pi)
        theta = rem(theta - pi, 2*pi) - pi;       
    elseif (theta < -pi)
        theta = rem(theta + pi, 2*pi) + pi;
    end

end