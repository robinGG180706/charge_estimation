function [R1, C1, tau, fit] = id_RC_from_one_step(I, V, k0, Nrelax, dt, R0, Npre)
% ID_RC_FROM_ONE_STEP  Identification de l'étage RC sur UN step de courant
%
% Step entre k0-1 (I_prev) et k0 (I_new), courant ensuite ~constant.
%
% Modèle après le step :
%   V(t) = V_inf + A * exp( -t / tau )
%
% On ajuste V_inf, A et tau par moindres carrés NON linéaires
% directement sur V, puis on déduit R1 à partir des deux plateaux
% (avant / après step).
%
% INPUTS :
%   I, V   : vecteurs courant / tension (même longueur)
%   k0     : indice du 1er échantillon APRÈS le step
%   Nrelax : nb de points à utiliser après k0
%   dt     : pas d'échantillonnage [s]
%   R0     : (optionnel) résistance ohmique [Ω] (0 si inconnu)
%   Npre   : (optionnel) nb de points du plateau avant step (défaut 10)
%
% OUTPUTS :
%   R1, C1, tau : paramètres de l'étage RC
%   fit         : structure pour vérifier le fit (V_inf, t_rel, V_model, k0)

    if nargin < 6 || isempty(R0)
        R0 = 0;
    end
    if nargin < 7 || isempty(Npre)
        Npre = 10;
    end

    I = I(:);
    V = V(:);
    N = length(I);
    if length(V) ~= N
        error('I et V doivent avoir la même longueur');
    end

    k_end = k0 + Nrelax - 1;
    if k0 < 2 || k_end > N
        error('Fenêtre de relaxation hors limites ou k0 < 2');
    end

    % --- 1) indices de relaxation et plateau avant step ---
    idx_relax = (k0:k_end).';
    Nrelax    = length(idx_relax);

    pre_start = max(1, k0 - Npre);
    pre_end   = k0 - 1;
    if pre_end <= pre_start
        error('Pas assez de points avant le step pour le plateau précédent');
    end
    idx_pre = (pre_start:pre_end).';

    % --- 2) données pour le fit ---
    t_rel   = (0:Nrelax-1).' * dt;      % temps relatif après step
    V_relax = V(idx_relax);

    % estimation grossière de V_inf (moyenne de la fin de la fenêtre)
    n_tail   = max(5, round(0.3 * Nrelax));
    idx_tail = (Nrelax-n_tail+1):Nrelax;
    V_inf0   = mean(V_relax(idx_tail));
    A0       = V_relax(1) - V_inf0;
    tau0     = max(dt, Nrelax*dt/3);    % estimation grossière de tau

    % --- 3) fit non linéaire V(t) = p1 + p2*exp(-t / (p3^2)) ---
    % p3^2 garantit tau > 0
    model = @(p,t) p(1) + p(2).*exp(-t./(p(3).^2));

    errfun = @(p) sum( (V_relax - model(p,t_rel)).^2 );

    p0 = [V_inf0, A0, sqrt(tau0)];
    opts = optimset('Display','off');
    p_opt = fminsearch(errfun, p0, opts);

    V_inf = p_opt(1);
    A     = p_opt(2);
    tau   = p_opt(3)^2;

    % --- 4) estimation des plateaux et de R1 ---
    V_before = mean( V(idx_pre) );
    I_prev   = mean( I(idx_pre) );

    % plateau après step : on considère que V_inf est la tension finale,
    % et que le courant après step est le courant moyen sur les premiers points
    n_head   = min(5, Nrelax);
    I_new    = mean( I(idx_relax(1:n_head)) );

    dI    = I_new - I_prev;
    if abs(dI) < 1e-6
        error('ΔI trop petit pour estimer R1');
    end

    % variation de tension entre les plateaux (steady-state)
    DeltaV_ss = (V_inf - V_before) - R0*dI;
    R1 = DeltaV_ss / dI;
    C1 = tau / R1;

    % forcer les signes physiques
    if C1 < 0
        R1 = -R1;
        C1 = -C1;
    end

    % --- 5) calcul de la courbe modélisée pour affichage ---
    V_model = model(p_opt, t_rel);

    fit = struct();
    fit.V_inf   = V_inf;
    fit.t_rel   = t_rel;
    fit.V_model = V_model;
    fit.k0      = k0;
end



function results = id_RC_from_pulses(I, V, dt, R0, dI_min, Nrelax_max)

    if nargin < 4 || isempty(R0)
        R0 = 0;
    end
    if nargin < 5 || isempty(dI_min)
        dI_min = 1;      % seuil de step (A)
    end
    if nargin < 6 || isempty(Nrelax_max)
        Nrelax_max = 50; % nb max de points de relax
    end

    I = I(:);
    V = V(:);
    N = length(I);
    if length(V) ~= N
        error('I et V doivent avoir la même longueur');
    end

    % --- 1) détection des steps sur le courant ---
    dI = diff(I);
    idx_steps = find(abs(dI) > dI_min) + 1;   % premier point après step

    if isempty(idx_steps)
        error('Aucune marche de courant détectée (baisser dI_min ?)');
    end

    R1_all  = [];
    C1_all  = [];
    tau_all = [];

    % *** CORRECTION ICI : initialiser un tableau de structs AVEC les bons champs ***
    fits = struct('V_inf', {}, ...
                  'DeltaV_sample', {}, ...
                  'k_relax', {}, ...
                  'V_model', {}, ...
                  'k0', {});

    for s = 1:length(idx_steps)
        k0 = idx_steps(s);

        % fenêtre maxi jusqu'au prochain step ou fin de série
        if s < length(idx_steps)
            k_next_step = idx_steps(s+1);
            Nrelax = min(Nrelax_max, k_next_step - k0);
        else
            Nrelax = min(Nrelax_max, N - k0);
        end

        if Nrelax < 10
            % trop court pour un fit propre
            continue;
        end

        try
            [R1, C1, tau, fit] = id_RC_from_one_step(I, V, k0, Nrelax, dt, R0);

            R1_all(end+1,1)  = R1;
            C1_all(end+1,1)  = C1;
    tau_all(end+1,1) = tau;

    fits(end+1) = fit;
catch ME
    warning('Step à k0=%d ignoré : %s', k0, ME.message);
end
   
    end

    if isempty(R1_all) || isempty(fits)
        error('Aucun step exploitable trouvé');
    end

    results = struct();
    results.R1_all   = R1_all;
    results.C1_all   = C1_all;
    results.tau_all  = tau_all;
    results.R1_mean  = mean(R1_all);
    results.tau_mean = mean(tau_all);
    results.C1_mean  = results.tau_mean / results.R1_mean;
    results.steps    = fits;
end

results = id_RC_from_pulses(current_vector, voltage_vector, dt, R0, dI_min, Nrelax_max);

fprintf('=== Moyenne sur tous les steps détectés ===\n');
fprintf('R1 moyen  = %.6g ohms\n', results.R1_mean);
fprintf('C1 moyen  = %.6g F\n',   results.C1_mean);
fprintf('tau moyen = %.6g s\n',   results.tau_mean);

if ~isempty(results.steps)
    fit1 = results.steps(1);
    idx_relax = fit1.k0 + fit1.k_relax;

    figure;
    plot(idx_relax*dt, voltage_vector(idx_relax), 'o', 'DisplayName', 'V mesurée');
    hold on;
    plot(idx_relax*dt, fit1.V_model, '--', 'DisplayName', 'V modèle RC');
    xlabel('Temps [s]');
    ylabel('Tension [V]');
    grid on;
    legend('Location','best');
    title('Relaxation tension sur un step de courant');
end