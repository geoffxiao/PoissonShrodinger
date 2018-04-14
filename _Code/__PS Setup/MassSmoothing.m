function out = MassSmoothing(masses, smoothing_distance, smoothing_factor, interfaces, x)
    % Mass smoothing at the interfaces
    % Given abrupt interfaces, smooth the interfaces using tanh(x)

    % smoothing_distance; % how much the smoothing (index-wise) extends into a layer for

    interfaces_smoothing = [0; interfaces; numel(x)];
    out = [];
    
    for i = 1 : numel(interfaces) + 1
    
        out = [out; repmat(masses(i), interfaces_smoothing(i + 1) - interfaces_smoothing(i), 1)];
        
    end
 
    for i = 1 : numel(interfaces)

        % sharp interface
        mass_here = repmat(masses(i), interfaces_smoothing(i + 1) - interfaces_smoothing(i), 1);
        mass_next = repmat(masses(i + 1), interfaces_smoothing(i + 2) - interfaces_smoothing(i + 1), 1);
        mass_combined = [mass_here; mass_next];
        indices = 1 : numel(mass_combined);
        
        % smooth it now
        center = interfaces(i);
        amplitude = (masses(i + 1) - masses(i)) / 2;
        to_eval = abs( indices - center) < smoothing_distance;
        out( to_eval ) = ...
            amplitude * tanh( (indices(to_eval) - center) * smoothing_factor) + ( masses(i) + masses(i + 1)) * 0.5;
        
    end

end