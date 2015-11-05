function outflag = CoagOutput(p, p2, t_out, spec)
%
% This function does some simple graphical output
%

outflag = 0;

%% Book-keeping
%  Here we calculate some things for the graphical output

total_spec = sum(spec');

axis_limits = [0, p.n_sections, min(total_spec), max(total_spec)];

n_times = length(t_out);

%% Calculate additional spectra
%  Calculate the number spectrum, mass spectrum and flux with respect to
%  particle diameters.

nspec_v    = zeros(n_times, p.n_sections);
massspec_v = nspec_v;
fluxsect   = nspec_v;
fluxspec   = nspec_v;

r_i = p.amfrac *p.av_vol.^p.bmfrac;
r_v = (0.75/pi*p.av_vol).^(1.0/3.0);

set_vel = SettlingVelocity(r_i, r_v, p.setcon);
set_vel = set_vel/100*p.day_to_sec;

diam_i = 2.0*p.r_to_rg*r_i;
diam_v = 2.0*r_v;

diam_i = diam_i';
diam_v = diam_v';

diam_i_mat = diam_i(ones(n_times, 1), :);
diam_v_mat = diam_v(ones(n_times, 1), :);

for jindx = 1 : n_times
    
    yout = spec(jindx,:);
    nspec_v(jindx,:)   = yout./(1.5*p.v_lower')./p.dwidth';
    masspec_v(jindx,:) = yout./p.dwidth';
    fluxsect(jindx,:)  = yout.*set_vel'*1e6;
    fluxspec(jindx,:)  = masspec_v(jindx,:).*set_vel'*1e6;
    
end
total_flux = sum(fluxsect,2);
total_mass = sum(spec, 2);

diaratio = (p.fr_dim/3)*diam_v_mat./diam_i_mat;
nspec_i = nspec_v.*diaratio;
masspec_i = masspec_v.*diaratio;
fluxspec_i = fluxspec.*diaratio;

%keyboard

%% Simple 2D plots

fig_h1 = figure(1);

sp_h_221 = subplot(2,2,1)
loglog(diam_i, nspec_i);
xlabel('Particle diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
axis tight

sp_h_223 = subplot(2,2,3);
plot(diam_i_mat(:,2:end)', fluxspec_i(:,2:end)')
xlabel('Particle image diameter [cm]')
ylabel('Volume flux spectra [cm^2 m^{-2} d^{-1}]')
axis tight

sp_h_247 = subplot(2,4,7);
plot(t_out, fluxsect, t_out, total_flux, '*--')
xlabel('Time [d]')
ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]')

sp_h_248 = subplot(2,4,8);
plot(t_out, total_flux./total_mass/1e6)
xlabel('Time [d]')
ylabel('Average v [m d^{-1}]')

sp_h_243 = subplot(2,4,3)
semilogy(t_out, spec, t_out, total_mass, '*--')
xlabel('Time [d]')
ylabel('Sectional concentration [vol/vol/sect]')



