function [d] = Build_d(m,PastPeriods,T_full,Shocks_u,past_obj)

    x_full = zeros(m.nx,PastPeriods+T_full);
    shks_full = zeros(m.ns,PastPeriods+T_full);
    shks_full(:,1:size(Shocks_u,2)) = Shocks_u;

    if ~isstruct(past_obj)
        t0_inherited = past_obj;
        if length(t0_inherited) ~= m.nx
            error('Wrong size for vector of inherited values from period 0');
        end
        for t = 1:PastPeriods
            if t == 1
                x_full(:,t) = m.PHI*shks_full(:,t) + t0_inherited;
            else
                x_full(:,t) = m.PHI*shks_full(:,t) + m.B*x_full(:,t-1);
            end
        end
        d.Past.Deviations = x_full(:,1:PastPeriods);
        d.Past.Shocks_u   = shks_full(:,1:PastPeriods);
        d.Past.Trends     = repmat(m.SS,1,PastPeriods);
    else
        d_Past = past_obj;
        if size(d_Past.Past.Deviations,1) <= m.nx
            d.Past = d_Past.Past;
        else
            error('Trouble!');
        end

        if size(d_Past.Past.Deviations,1) < m.nx
            d.Past.Deviations = [d.Past.Deviations ; zeros(m.nx-size(d_Past.Past.Deviations,1),PastPeriods)];
            d.Past.Trends = repmat(m.SS,1,PastPeriods);
        end
        
        x_full(:,1:PastPeriods) = d.Past.Deviations;
    end
    
    for t = (PastPeriods+1):(PastPeriods+T_full)
        x_full(:,t) = m.PHI*shks_full(:,t) + m.B*x_full(:,t-1);
    end
    
    d.Forecast.Deviations = x_full(:,(PastPeriods+1):end);
    d.Forecast.Shocks_u   = shks_full(:,(PastPeriods+1):end);
    d.Forecast.Trends     = repmat(m.SS,1,T_full);
end