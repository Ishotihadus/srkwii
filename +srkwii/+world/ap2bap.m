function bap = ap2bap(ap, frequency_axis, coarse_axis)
% ap2bap - convert full aperiodicity into band aperiodicity

coarse_axis = coarse_axis(:);
frequency_axis = frequency_axis(:);
ap = max(min(ap, 1), 1e-3);
bap = cell2mat(...
  cellfun(@(b) interp1(frequency_axis, 20 * log10(b), coarse_axis, 'linear'), num2cell(ap, 1), 'un', 0));
