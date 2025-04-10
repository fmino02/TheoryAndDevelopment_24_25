function status = progress_output(t, y, flag)
    % Declare persistent variables to retain across calls
    persistent h t_start t_final

    status = 0;  % Required return value (0 = continue, 1 = stop)

    switch flag
        case 'init'
            % Initialization: create waitbar
            t_start = t(1);
            t_final = evalin('base', 'tfinal');  % Read final time from base workspace
            h = waitbar(0, sprintf('Simulation Progress: t = %.2f', t_start));
        
        case ''
            % Regular update
            if ~isempty(t) && ishandle(h)
                frac_done = (t(end) - t_start) / (t_final - t_start);
                waitbar(frac_done, h, sprintf('Simulation Progress: t = %.2f', t(end)));
            end
        
        case 'done'
            % Close waitbar when finished
            if ishandle(h)
                waitbar(1, h, 'Simulation complete!');
                pause(0.5);  % So you see the "100%" before closing
                close(h);
            end
            clear h t_start t_final
    end
end
