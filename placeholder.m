for k = 1:N
    t = (k-1)*dt;

    % 1. MinSnap Referanslarını Çek
    x_ref = xref(:,k);
    v_ref = vref(:,k);
    a_ref = aref(:,k);
    j_ref = jref(:,k);
    s_ref = sref(:,k);

    % 2. YAW PLANNER ÇAĞRISI (Bunu eklemeyi unutma!)
    st_ref.v = v_ref;
    st_ref.a = a_ref;
    
    % Bir önceki adımdaki psi değerini hatırla (Süreklilik için)
    if k > 1
        params.last_psi = log.omega_ref(3, k-1) * dt + log.q(4, k-1)*2; % Yaklaşık veya logdan çek
        % Daha doğrusu: psi_ref'i loglayıp oradan çekmektir.
        % Şimdilik basitçe şunu yapalım:
        params.last_psi = psi_ref; % Bir önceki döngüden kalan değer
    else
        params.last_psi = 0;
    end

    % Tangent modunda çalıştır
    [psi_ref, psi_dot_ref, psi_ddot_ref] = yaw_planner(t, 'tangent', st_ref, params);

    %% Dynamics ... (kodun geri kalanı aynı)
