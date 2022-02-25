function phase = create_movie_phase(action_potential,v_gate)

a = find( action_potential > v_gate );
if ~isempty(a)
    b = find( abs(diff(a)) > 1);
    p = [a(1); a(b); a(b+1); a(end)];
    p = sort(p,'ascend');

    debug_plot = 0;
    if debug_plot == 1
        figure;
        plot(action_potential,'b');
        hold on;
        scatter(p,action_potential(p),100,'.r');
        hold off;
    end

    % phase interval
    phase_interval = zeros(length(p)/2,2);
    for n = 1:length(p)/2
        phase_interval(n,:) = p([(n-1)*2+1 (n-1)*2+2]);
    end

    % phase
    l_median = ceil( median(phase_interval(:,2)-phase_interval(:,1)) );
    phase = zeros(size(action_potential));
    for n = 1:size(phase_interval,1)    
        m = phase_interval(n,:);
        l = length(m(1):m(2)-1);

        if m(1) == 1
            phase(m(1):m(2)-1) = linspace(m(1)/l_median,1,l);
        end

        if m(1) ~= 1 && m(2) ~= length(phase)        
            phase(m(1):m(2)-1) = linspace(0,1,l);
        end

        if m(2) == length(phase)
            phase(m(1):m(2)-1) = linspace(0,l/l_median,l);
        end
    end
elseif isempty(a)
    phase = zeros(size(action_potential));
end

debug_plot = 0;
if debug_plot == 1
    figure;
    plot(action_potential,'b');
    hold on;
    plot(phase,'r');
    hold off;
    axis tight;
end

end