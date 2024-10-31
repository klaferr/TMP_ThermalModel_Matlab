au  = 1.4959787061e11;           % Size of one astronomical unit in meters
gc  = 6.67259e-11;               % Gravitational constant
sm  = 1.9891e30;                 % Solar mass
yearLength = 2*pi/sqrt(gc*sm/(s.a*au)^3);   % length of Mars year in seconds using Kepler's 3rd law
EarthYearLength = 2*pi*sqrt(au^3/(gc*sm));   % Length of one Earth year in seconds

TotalSurfSublimation = 0; %G_FFsubl;
lagGrowthAnnConv = TotalSurfSublimation*s.iceDust /((1-s.iceDust)*(1-s.regolithPorosity));

%% DIFFUSIVE LOSS THROUGH LAG
% Save final lag thickness to use as starting point for next run
% Does this work well though if saves as z_ei which means no pf? check
% this later....
lagThickness = z_ei;
old_lagThickness = z_ei;

meanTsurf = mean(Tsurf); % mean surface temperature
atmosDensitySF12 = atmosPressSF12 * PtoRhoFactor / meanTsurf;

% Diffusive case
if z_ei > 0
    
    % Calculate T's and rho_v at the excess ice sheet boundary (remembering
    % now that z_ei_index is for the last layer of the lag deposit
    % and z_ei_index+1 is first layer of ice sheet.
    % layerboundary(iceTableIndex) is the depth at the ice table (top of
    % the ice layer)
    iceTableTemps = (k(z_ei_index+1).*dz(z_ei_index).*Temps(z_ei_index+1,:) + k(z_ei_index).*dz(z_ei_index+1).*Temps(z_ei_index,:))/(k(z_ei_index+1)*dz(z_ei_index) + k(z_ei_index)*dz(z_ei_index+1));
    iceTableTemps = iceTableTemps';
    iceTable_Pv = 611 .* exp( (-51058/8.31).*(1./iceTableTemps - 1/273.16) ); % compute vapor pressure at ice table
    iceTable_rhov = iceTable_Pv .* (0.01801528 / 6.022140857e23) ./ (1.38064852e-23 .* iceTableTemps); % compute vapor densities at ice table
    
    iceTable_rhov_mean = mean(iceTable_rhov); % get mean for that year of ice table vapor densities
    iceTableTemps_mean = mean(iceTableTemps); % mean temperature at the ice table over the year
    
    if PFICE == 0
        % This is the case when pf ice isn't stable so the ice sheet
        % can retreat and if that's the case, grow the lag deposit
        
        % Calculate new lag thickness for next run
        % This formulation is from the integration in Schorghofer 2010 that
        % has a 1/lagThickness in the equation, so this breaks down if the
        % lag is very thin, if the timestep is too large, or if the vapor
        % density difference is large or changing between steps.
        if s.growLag == 'y'
            % Updating to include lag stripping term
            %lagThickness = sqrt(lagThickness^2 + (((2*porosity_param*s.diffusionCoef*(iceTable_rhov_mean - atmosDensitySF12)/densityIce)-(s.removeLag/EarthYearLength))*(s.Laskar_step*1000*EarthYearLength)))  %how much occurs over the 1000 years of laskar solution (earth years...could also use yearLength/1.88)
            lagThickness = sqrt(lagThickness^2 + (((2*porosity_param*s.diffusionCoef*(iceTable_rhov_mean - atmosDensitySF12)/densityIce))*(s.Laskar_step*1000*EarthYearLength)))  %how much occurs over the 1000 years of laskar solution (earth years...could also use yearLength/1.88)

            aveLagThickness = .5*(lagThickness + old_lagThickness);
            
            % Vapor density flux
            annualVaporFlux = s.diffusionCoef * (iceTable_rhov_mean - atmosDensitySF12) /  aveLagThickness;
            
        elseif s.growLag == 'n'
            annualVaporFlux = s.diffusionCoef * (iceTable_rhov_mean - atmosDensitySF12) /  lagThickness
        end
        % Annual growth rate of lag in diffusive case
        lagGrowthAnnDiff = (lagThickness - old_lagThickness)/(s.Laskar_step*1000);
        
        % Change in ice table, accounts for excess ice porosity and dirt content
        DiffusiveIceLoss = annualVaporFlux * EarthYearLength / (densityIce * (1-s.icePorosity-s.iceDust)) % annual retreat of ice (in Earth years)
        
        if DiffusiveIceLoss < TotalSurfSublimation
            % if diffusive loss is smaller than surface sublimation than
            % this will be the ice loss for the timestep, otherwise will
            % stay at the total surf sublimation value set above
            dz_iceLoss = DiffusiveIceLoss;
        else
            dz_iceLoss = TotalSurfSublimation;
        end
        
    else
        % pore filling ice is stable - assume no retreat from ice sheet
        DiffusiveIceLoss = 0; 
        dz_iceLoss = DiffusiveIceLoss;
        lagGrowthAnnDiff = 0;
        % lagThickness remains unchanged
    end
    
else
% No lag

    iceTable_rhov_mean = -1;
    iceTableTemps_mean = -1;
    DiffusiveIceLoss = -1; % -1 means there is no lag to even calculate diffusive loss
    lagGrowthAnnDiff = -1;
    dz_iceLoss = TotalSurfSublimation;
    
    % Grow lag from surface sublimation instead
    if s.growLag == 'y'
        lagThickness = lagGrowthAnnConv*(s.Laskar_step*1000);
        lagGrowthAnnDiff = lagThickness;
    end
end

dz_iceLossTimestep(counterOrbital) = dz_iceLoss * (s.Laskar_step*1000);  % each timestep in Laskar is 1000 years so assume same conditions for that long to get retreat over each time step, but Laskar is in Earth years so divide by 1.88 because Mars years are longer...but just used Earth year time in seconds so its all good
dz_CumulateIceLoss = cumsum(dz_iceLossTimestep); % sum up retreat over time
delta_lagThickness = lagThickness - old_lagThickness;

fprintf('\nFinal sublimation from diffusive loss per year is %10.9f mm\n', DiffusiveIceLoss*1000);
fprintf('\nTotal sublimation for this timestep %10.9f m\n', dz_iceLossTimestep(counterOrbital));

if s.removeLag > 0
    fprintf('\nLag thickness is %20.9f and removal will be %20.9f so new z_ei is %20.9f\n', lagThickness, (s.removeLag * (s.Laskar_step*1000)), lagThickness - (s.removeLag * (s.Laskar_step*1000)));
	z_ei = lagThickness - (s.removeLag * (s.Laskar_step*1000)); % Set new z_ei but remove some lag  
    lagThickness = z_ei;
    amtLagRemoved = (s.removeLag * (s.Laskar_step*1000));
    
    % Make sure lag doesn't get below h0 value
    if h0 > 0  % Is there an h0 value
        if z_ei < h0
            % Set to h0 if there is one and z_ei drops below that
            fprintf('\nz_ei of %20.9f is less than h0 of %20.9f\n', z_ei, h0);
            z_ei = h0;
            lagThickness = z_ei; 
        end
    else
        if z_ei < s.thinnestLayer % Somewhat arbitrary but want to make sure the code doesn't try to use a super small lag in the cases where there is no minimum h0 value
            % Just set to 0, this should be the cases where convective loss
            % would be pretty small anyways...let's see how it goes
            fprintf('\nz_ei of %20.9f is less than 0.0005 so setting to 0\n', z_ei);
            z_ei = 0;
            lagThickness = z_ei;
        end
    end
else
    z_ei = lagThickness; % Set new z_ei
    fprintf('\nSetting z_ei to lag thickness without removing an lag: %20.9f\n', z_ei);
    amtLagRemoved = 0;
end

z_pf_old = z_pf;
if PFICE == 0
    % If pore-filling ice isn't stable, set it equal to z_ei to come in at
    % the next timestep as if without pore filling ice
    z_pf = z_ei;
else
    % check if z_pf now bigger than the new z_ei after lag has been
    % stripped and if so, also start z_pf at z_ei, otherwise can come in
    % with the pf that it's already at
    if z_pf > z_ei
        z_pf = z_ei
    end
end


% Calculate horizontal retreat if sloped - REMEMBER TO USE SIN in DEGREES
if s.slope > 0
    dz_iceLoss_horiz = dz_iceLoss/sind(s.slope);
else
    dz_iceLoss_horiz = dz_iceLoss;
end
dz_iceLossTimestep_horiz(counterOrbital) = dz_iceLoss_horiz * (s.Laskar_step*1000); 
dz_CumulateIceLoss_horiz = cumsum(dz_iceLossTimestep_horiz);
    
