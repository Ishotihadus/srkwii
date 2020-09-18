function hz = mel2hz(mel)
% mel2hz - convert from hertz to mel scale
%   hz = mel2hz(mel)

hz = 700 * (10 .^ (mel / 2595) - 1);
