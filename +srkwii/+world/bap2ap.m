function ap = bap2ap(bap, coarse_axis, frequency_axis)
% bap2ap - convert band aperiodicity into full aperiodicity

coarse_axis = coarse_axis(:);
frequency_axis = frequency_axis(:);
bap = max(min(bap, 0), -60);
ap = cell2mat(...
  cellfun(@(b) 10 .^ (interp1(coarse_axis, b, frequency_axis, 'linear') / 20), num2cell(bap, 1), 'un', 0));
