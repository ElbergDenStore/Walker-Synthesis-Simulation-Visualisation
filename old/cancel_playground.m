function cancel_playground()
%CANCEL_PLAYGROUND  Empirically tests parfeval cancellation strategies.
%
%   cancel_playground()
%
% Tests three concrete patterns and prints timings so you can see what
% actually works.  Run this before touching gridsearch.
%
%   Test A  — worker-created PollableDataQueue, client tries to send back
%   Test B  — client-created PollableDataQueue, passed to worker as arg
%   Test C  — per-task parfeval; cancel(queued_futures) vs cancel(running)

    pool = gcp('nocreate');
    if isempty(pool), pool = parpool(4); end
    fprintf('\nPool: %d workers (Processes pool)\n', pool.NumWorkers);

    test_A_worker_created_queue(pool);
    test_B_client_created_queue(pool);
    test_C_cancel_queued_vs_running(pool);

    fprintf('\n--- Done ---\n');
end

% =========================================================================

function test_A_worker_created_queue(pool)
% HYPOTHESIS: worker creates PollableDataQueue, sends handle to client;
%             client calls send() on it to signal the worker.
%             PollableDataQueue docs say "client polls, workers send" —
%             the reverse direction may silently fail.
    fprintf('\n=== Test A: WORKER-created queue (client sends back) ===\n');

    q_results = parallel.pool.PollableDataQueue;
    f = parfeval(pool, @worker_A, 0, q_results, 8);

    % Collect worker's queue handle
    q_worker = [];
    t_reg = tic;
    while isempty(q_worker) && toc(t_reg) < 5
        [msg, got] = poll(q_results, 0.1);
        if got && isfield(msg, 'registered')
            q_worker = msg.queue;
        end
    end
    if isempty(q_worker)
        fprintf('  FAIL: worker never sent its queue handle\n');
        cancel(f); return;
    end
    fprintf('  Worker registered in %.2f s\n', toc(t_reg));

    % Let worker run 2 iterations then signal stop
    pause(2.2);
    t_send = tic;
    send(q_worker, struct('stop', true));
    fprintf('  send() returned in %.4f s\n', toc(t_send));

    % Did the worker actually stop?
    [resp, got] = poll(q_results, 3);
    if got && isfield(resp, 'stopped')
        fprintf('  Worker stopped at iter %d — WORKS (%.2f s)\n', resp.iter, resp.elapsed);
    else
        fprintf('  Worker did NOT stop within 3 s — BROKEN\n');
    end
    cancel(f);
    wait(f);
end

function worker_A(q_out, n_iters)
    my_q = parallel.pool.PollableDataQueue;
    send(q_out, struct('registered', true, 'queue', my_q));
    t0 = tic;
    for i = 1:n_iters
        [msg, got] = poll(my_q, 0);
        if got && isfield(msg, 'stop')
            send(q_out, struct('stopped', true, 'iter', i, 'elapsed', toc(t0)));
            return;
        end
        pause(1);
    end
end

% =========================================================================

function test_B_client_created_queue(pool)
% CORRECT direction: client creates PollableDataQueue, passes to worker
% as a function argument.  Worker polls it.  Client sends to it.
    fprintf('\n=== Test B: CLIENT-created queue (passed to worker) ===\n');

    q_results = parallel.pool.PollableDataQueue;
    q_control = parallel.pool.PollableDataQueue;   % created on client

    f = parfeval(pool, @worker_B, 0, q_results, q_control, 8);

    % Let worker run 2 iterations then signal stop
    pause(2.2);
    t_send = tic;
    send(q_control, struct('stop', true));
    fprintf('  send() returned in %.4f s\n', toc(t_send));

    [resp, got] = poll(q_results, 3);
    if got && isfield(resp, 'stopped')
        fprintf('  Worker stopped at iter %d — WORKS (%.2f s)\n', resp.iter, resp.elapsed);
    else
        fprintf('  Worker did NOT stop within 3 s — BROKEN\n');
    end
    cancel(f);
    wait(f);
end

function worker_B(q_out, q_control, n_iters)
    t0 = tic;
    for i = 1:n_iters
        [msg, got] = poll(q_control, 0);
        if got && isfield(msg, 'stop')
            send(q_out, struct('stopped', true, 'iter', i, 'elapsed', toc(t0)));
            return;
        end
        pause(1);
    end
end

% =========================================================================

function test_C_cancel_queued_vs_running(pool)
% CORE QUESTION: how fast is cancel() on QUEUED vs RUNNING futures?
%
% If queued-cancel is instant → per-task architecture gives free cancellation.
% If running-cancel blocks   → avoid waiting; fire-and-forget.
    fprintf('\n=== Test C: cancel() timing — queued vs running ===\n');

    N = 60;
    nw = pool.NumWorkers;
    fprintf('  Submitting %d tasks to %d workers...\n', N, nw);

    futs(N) = parallel.FevalFuture;
    for i = 1:N
        futs(i) = parfeval(pool, @slow_task, 0, 30);  % each "runs" for 30 s
    end
    pause(1);   % let pool start `nw` of them

    states    = {futs.State};
    n_running = sum(strcmp(states, 'running'));
    n_queued  = sum(strcmp(states, 'queued'));
    fprintf('  %d running, %d queued\n', n_running, n_queued);

    % --- cancel queued ---
    q_mask = strcmp({futs.State}, 'queued');
    t = tic;
    cancel(futs(q_mask));
    t_q = toc(t);
    fprintf('  cancel(%d queued): %.3f s\n', sum(q_mask), t_q);

    % --- cancel running, do NOT wait ---
    r_mask = strcmp({futs.State}, 'running');
    t = tic;
    cancel(futs(r_mask));
    t_r_noWait = toc(t);
    fprintf('  cancel(%d running) without wait: %.3f s\n', sum(r_mask), t_r_noWait);

    % --- now wait to see how long before they actually finish ---
    t = tic;
    wait(futs(r_mask));
    t_actual = toc(t);
    fprintf('  wait(running after cancel): %.3f s (workers run until next safe point)\n', t_actual);

    fprintf('\n  CONCLUSION:\n');
    fprintf('  - cancel(queued) : %.3f s  → instant, use freely\n', t_q);
    fprintf('  - cancel(running): %.3f s (call) + %.3f s (actual stop)\n', t_r_noWait, t_actual);
    if t_r_noWait < 0.1
        fprintf('  - cancel() call is NON-BLOCKING → safe to fire-and-forget\n');
    else
        fprintf('  - cancel() call BLOCKS → avoid on running futures\n');
    end
end

function slow_task(duration_s)
    % Simulates an "atomic" computation (like states()) that cannot be
    % interrupted from outside.  In practice this is replaced by the
    % satellite toolbox's states() call.
    pause(duration_s);
end
