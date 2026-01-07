function [psi, psi_dot, psi_ddot] = yaw_planner(t, mode, state_ref, params)
%YAW_PLANNER Yaw referanslarını hesaplar.
%   mode: 'constant', 'tangent', 'coordinated'

    switch lower(mode)
        case "constant"
            psi = 0; psi_dot = 0; psi_ddot = 0;
            
        case {"tangent", "coordinated"} 
            % 'tangent' ve 'coordinated' burada aynı matematiksel temeli kullanır
            % Amaç: Gövde X eksenini Hız vektörüne kilitlemek.
            
            vx = state_ref.v(1);
            vy = state_ref.v(2);
            vz = state_ref.v(3);
            
            ax = state_ref.a(1);
            ay = state_ref.a(2);
            az = state_ref.a(3);
            
            v_sq = vx^2 + vy^2 + vz^2; % 3D Hız büyüklüğü
            
            % Hız eşiği: Çok yavaşsa (Hover) yön değiştirmeye çalışma
            if v_sq < 0.1 
                if isfield(params, 'last_psi')
                    psi = params.last_psi;
                else
                    psi = 0;
                end
                psi_dot = 0; 
                psi_ddot = 0;
            else
                %% KOORDİNELİ YAW HESABI (Doğal Takip)
                % Adım 1: İstenen Thrust Yönü (Body Z - b3)
                % a_ref = T_vec + g -> T_vec = a_ref - g
                % Drone Z ekseni Thrust'ın tersinedir.
                g_vec = [0;0;9.81];
                t_vec = [ax; ay; az] - g_vec; 
                
                % Eğer serbest düşüşteyse (t_vec=0), Z eksenini aşağı kabul et
                if norm(t_vec) < 1e-3
                    b3 = [0;0;1];
                else
                    b3 = -t_vec / norm(t_vec);
                end
                
                % Adım 2: İstenen Gidiş Yönü (Velocity Unit Vector)
                vel_dir = [vx; vy; vz] / sqrt(v_sq);
                
                % Adım 3: Body Y Ekseni (b2)
                % Hız vektörü ile Thrust vektörünün oluşturduğu düzleme diktir.
                % Bu, drone'un "yatış" (bank) eksenidir.
                y_temp = cross(b3, vel_dir);
                
                % Singularity Check: Eğer hız vektörü ile Z ekseni paralelse (Tam dik dalış/tırmanış)
                if norm(y_temp) < 1e-3
                     % Standart bir Y ekseni seç (Örn: Doğu)
                     b2 = [0;1;0];
                else
                     b2 = y_temp / norm(y_temp);
                end
                
                % Adım 4: Body X Ekseni (b1) - Doğal Burun Yönü
                % Z ve Y belliyse, X bunların cross product'ıdır.
                b1 = cross(b2, b3);
                
                % Adım 5: Rotasyon Matrisinden Yaw Çıkarma
                % R_des = [b1, b2, b3]
                % Yaw (Psi) = atan2(R(2,1), R(1,1))
                target_psi = atan2(b1(2), b1(1));
                
                % AÇI SÜREKLİLİĞİ (Unwrap)
                if isfield(params, 'last_psi')
                    diff = target_psi - params.last_psi;
                    while diff > pi,  diff = diff - 2*pi; end
                    while diff < -pi, diff = diff + 2*pi; end
                    psi = params.last_psi + diff;
                else
                    psi = target_psi;
                end

                % TÜREVLER (Feedforward)
                % Coordinated turn için tam analitik türev çok karmaşıktır.
                % Ancak "tangent" türev yaklaşımı (yatay düzlem) burada da %95 oranında iş görür.
                % Çünkü Yaw, yerel Z ekseni etrafındaki dönüştür.
                
                v_hor_sq = vx^2 + vy^2;
                if v_hor_sq > 0.1
                     psi_dot = (vx*ay - vy*ax) / v_hor_sq;
                else
                     psi_dot = 0;
                end
                psi_ddot = 0; 
            end
            
        case "poi"
            % (Önceki POI kodu buraya gelecek...)
            psi = 0; psi_dot=0; psi_ddot=0; 
            
        otherwise
            psi = 0; psi_dot = 0; psi_ddot = 0;
    end
end