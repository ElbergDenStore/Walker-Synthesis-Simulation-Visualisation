classdef CancelToken < handle
% CANCELTOKEN  Cancellation handle for constellation_simulator.
%
%   token = CancelToken(queue, total_sats)
%
% Drains a parallel.pool.PollableDataQueue for cancellation messages from a
% coordinator and decides whether the current simulation should abort.
%
% Pass `@() token.check()` as the `should_cancel` argument to
% constellation_simulator.  After the simulation returns, read
% `token.ConsumedThreshold` to update the worker's local skip threshold.
%
% Recognised queue messages
%   struct('skip_above_sats', N)    -> cancel if total_sats > N
%   struct('cancel_detailed', true) -> always cancel
%
% Example (gridsearch worker)
%   token = CancelToken(q_in, Cfg.Total_sats);
%   m     = constellation_simulator(Cfg, false, false, true, @() token.check());
%   if m.cancelled
%       skip_above_sats = min(skip_above_sats, token.ConsumedThreshold);
%   end

    properties
        Queue                       % parallel.pool.PollableDataQueue
        TotalSats          double   % current Cfg.Total_sats
        ConsumedThreshold  double = Inf
        Cancelled          logical = false
    end

    methods
        function obj = CancelToken(queue, total_sats)
            obj.Queue     = queue;
            obj.TotalSats = total_sats;
        end

        function tf = check(obj)
        % Drain pending messages, update state, return true if simulation
        % should stop.  Cheap to call (just a non-blocking poll).
            while true
                [msg, got] = poll(obj.Queue, 0);
                if ~got, break; end
                if isfield(msg, 'skip_above_sats')
                    obj.ConsumedThreshold = min(obj.ConsumedThreshold, msg.skip_above_sats);
                end
                if isfield(msg, 'cancel_detailed') && msg.cancel_detailed
                    obj.Cancelled = true;
                end
            end
            if obj.TotalSats > obj.ConsumedThreshold
                obj.Cancelled = true;
            end
            tf = obj.Cancelled;
        end
    end
end
