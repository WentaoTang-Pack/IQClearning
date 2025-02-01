function g = util_bandpass(omega1, omega2, butter)
if nargin<3
    butter = 1;
end
if butter == 2
    if omega1 == 0 && omega2 == inf
        g = tf(1, 1);
    elseif omega1 == 0
        g = tf(1, [1/omega2^2 sqrt(2)/omega2 1]);
    elseif omega2 == inf
        g = tf([1/omega1^2 0 0], [1/omega1^2 sqrt(2)/omega1 1]);
    else
        g = tf([1/omega1^2 0 0], [1/omega1^2 sqrt(2)/omega1 1]) * tf(1, [1/omega2^2 sqrt(2)/omega2 1]) * (omega1^2/omega2^2 + 1);
    end
else
    if omega1 == 0 && omega2 == inf
        g = tf(1, 1);
    elseif omega1 == 0
        g = tf(1, [1/omega2 1]);
    elseif omega2 == inf
        g = tf([1/omega1 0], [1/omega1, 1]);
    else
        g = tf([1/omega1 0], [1/omega1, 1]) * tf(1, [1/omega2 1]) * (omega1/omega2 + 1);
    end
end
end