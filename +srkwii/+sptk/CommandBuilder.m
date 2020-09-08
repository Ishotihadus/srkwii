classdef CommandBuilder < matlab.mixin.Copyable
  properties (GetAccess = private, SetAccess = private)
    Path
    Options = struct('opt', {}, 'type', {}, 'val', {})
    Arguments = {}
    Output = '/dev/null'
  end

  methods
    function obj = CommandBuilder(command)
      obj.Path = checkcommandpath(srkwii.sptk.getsptkpath(command));
    end

    function AddOption(obj, opt, opttype, optval)
      n = numel(obj.Options) + 1;
      obj.Options(n).opt = opt;
      obj.Options(n).type = opttype;
      obj.Options(n).val = optval;
    end

    function AddArgument(obj, arg)
      n = numel(obj.Arguments) + 1;
      obj.Arguments{n} = arg;
    end

    function SetOutput(obj, output)
      obj.Output = output;
    end

    function Run(obj)
      command = obj.escape(obj.Path);
      for n = 1:numel(obj.Options)
        opt = obj.Options(n);
        switch opt.type
        case {'bool', 'boolean'}
          if opt.val
            command = sprintf('%s -%s', command, obj.escape(opt.opt)); end
        case {'int', 'integer'}
          command = sprintf('%s -%s %d', command, obj.escape(opt.opt), opt.val);
        case {'float', 'double', 'real'}
          command = sprintf('%s -%s %g', command, obj.escape(opt.opt), opt.val);
        case {'char', 'string'}
          command = sprintf('%s -%s %s', command, obj.escape(opt.opt), obj.escape(opt.val));
        otherwise
          error('unknown argument type `%s` for -%s', obj.type, obj.opt)
        end
      end
      for n = 1:numel(obj.Arguments)
        command = sprintf('%s %s', command, obj.escape(obj.Arguments{n})); end
      if ~isempty(obj.Output)
        command = sprintf('%s > %s', command, obj.escape(obj.Output)); end
      [status, ~] = system(command, '-echo');
      if status ~= 0
        error('command exited with status code %d: %s', status, command); end
    end

    function output = Exec1by1(obj, input, precision)
      if nargin <= 2
        precision = getsptknative; end

      infile = tempname;
      outfile = tempname;
      srkwii.io.writebin(infile, input, precision);
      obj.AddArgument(infile);
      obj.SetOutput(outfile);
      obj.Run;
      output = srkwii.io.readbin(outfile, precision);
      delete(infile);
      delete(outfile);
    end

    function output = Exec1by1Parallel(obj, input, batchnumel, precision)
      if nargin <= 3
        precision = getsptknative; end

      numworkers = ceil(numel(input) / batchnumel);
      if numworkers <= 1
        output = obj.Exec1by1(input, precision);
        return
      end

      builders = arrayfun(@(~) copy(obj), 1:numworkers, 'un', 0)';
      distsz = repmat(batchnumel, numworkers, 1);
      distsz(end) = numel(input) - batchnumel * (numworkers - 1);
      input_split = mat2cell(input(:), distsz, 1);
      output = cell(numworkers, 1);
      parfor widx = 1:numworkers
        output{widx, 1} = builders{widx, 1}.Exec1by1(input_split{widx, 1}, precision); end
      output = cell2mat(output);
    end
  end

  methods (Access = private, Static = true)
    function str = escape(str)
      if isa(str, 'string')
        str = convertStringsToChars(str); end
      str = ['''' strrep(str, '''', '''\\''''') ''''];
    end
  end
end
