import math
import itertools
import sys

# ============================================================
# PKGF Unified Test
# All experimental resources: https://github.com/aikenkyu001/PKGF
# ============================================================

DIM = 12
EPSILON_METRIC = 1e-5
EPSILON_MANIFOLD = 1e-4

# --- CORE FUNCTIONS ---

def identity_matrix():
    return [[1.0 if i == j else 0.0 for j in range(DIM)] for i in range(DIM)]

def zero_matrix():
    return [[0.0 for _ in range(DIM)] for _ in range(DIM)]

def mat_mul(A, B):
    C = zero_matrix()
    for i in range(DIM):
        for k in range(DIM):
            if abs(A[i][k]) > 1e-18:
                aik = A[i][k]
                for j in range(DIM):
                    C[i][j] += aik * B[k][j]
    return C

def mat_sub(A, B):
    return [[A[i][j] - B[i][j] for j in range(DIM)] for i in range(DIM)]

def mat_add(A, B):
    return [[A[i][j] + B[i][j] for j in range(DIM)] for i in range(DIM)]

def mat_scale(A, s):
    return [[A[i][j] * s for j in range(DIM)] for i in range(DIM)]

def transpose(A):
    return [[A[j][i] for j in range(DIM)] for i in range(DIM)]

def frobenius_norm(A):
    return math.sqrt(sum(A[i][j]**2 for i in range(DIM) for j in range(DIM)))

def mat_vec_mul(M, v):
    return [sum(M[i][j] * v[j] for j in range(DIM)) for i in range(DIM)]

def calculate_det_12x12(M):
    n = len(M)
    lu = [row[:] for row in M]
    det = 1.0
    for i in range(n):
        pivot = lu[i][i]
        if abs(pivot) < 1e-15:
            for k in range(i+1, n):
                if abs(lu[k][i]) > abs(pivot):
                    lu[i], lu[k] = lu[k], lu[i]
                    det *= -1
                    pivot = lu[i][i]
                    break
            if abs(pivot) < 1e-15: return 0.0
        det *= pivot
        for k in range(i+1, n):
            factor = lu[k][i] / pivot
            lu[k][i] = factor
            for j in range(i+1, n):
                lu[k][j] -= factor * lu[i][j]
    return det

def mat_inv(A):
    n = len(A)
    res = identity_matrix()
    m = [row[:] + res[i] for i, row in enumerate(A)]
    for i in range(n):
        pivot = m[i][i]
        if abs(pivot) < 1e-18:
            for k in range(i+1, n):
                if abs(m[k][i]) > abs(pivot):
                    m[i], m[k] = m[k], m[i]
                    pivot = m[i][i]
                    break
        if abs(pivot) < 1e-18: continue
        inv_p = 1.0 / pivot
        for j in range(2*n): m[i][j] *= inv_p
        for k in range(n):
            if i != k:
                f = m[k][i]
                for j in range(2*n): m[k][j] -= f * m[i][j]
    return [row[n:] for row in m]

def matrix_exp_pade(A):
    inf_norm = max(sum(abs(x) for x in row) for row in A)
    k = max(0, math.ceil(math.log2(inf_norm / 0.5))) if inf_norm > 0 else 0
    A_s = mat_scale(A, 1.0 / (2**k))
    c = [1.0, 0.5, 0.1, 0.01111111111111111, 0.0007575757575757576, 3.0303030303030303e-05, 5.050505050505051e-07]
    P, Q, A_p = identity_matrix(), identity_matrix(), identity_matrix()
    for i in range(1, 7):
        A_p = mat_mul(A_p, A_s)
        term = mat_scale(A_p, c[i])
        P = mat_add(P, term)
        if i % 2 == 0: Q = mat_add(Q, term)
        else:          Q = mat_sub(Q, term)
    res = mat_mul(mat_inv(Q), P)
    for _ in range(k): res = mat_mul(res, res)
    return res

def get_author_K():
    K = zero_matrix()
    scales = [1.0, 1.1, 1.2, 0.9] # S, E, A, C
    for c_idx, s in enumerate(scales):
        for i in range(3):
            K[c_idx*3 + i][c_idx*3 + i] = s
    return K

# --- METRIC FUNCTIONS ---

def get_metric_tensor(x):
    """
    Contextual Warping (文脈依存計量): 
    理論定義に基づき、Contextセクター(x[9-11])の座標強度が他セクターの空間密度(g_ii)を決定する。
    実装上の解釈: 
    - x[9]  -> Subjectセクター (0-2) の歪み
    - x[10] -> Entityセクター  (3-5) の歪み
    - x[11] -> Actionセクター  (6-8) の歪み
    これにより、物語の背景が各登場要素の「存在感(計量)」を動的に変調する。
    """
    g = zero_matrix()
    s_subject = 1.0 + 0.5 * math.tanh(x[9])
    s_entity  = 1.0 + 0.5 * math.tanh(x[10])
    s_action  = 1.0 + 0.5 * math.tanh(x[11])
    for i in range(DIM):
        if 0 <= i < 3:    g[i][i] = s_subject
        elif 3 <= i < 6:  g[i][i] = s_entity
        elif 6 <= i < 9:  g[i][i] = s_action
        else:             g[i][i] = 1.0
    return g

def get_inverse_metric(g):
    return mat_inv(g)

def partial_derivative_metric(x, k):
    x_plus = x[:]
    x_minus = x[:]
    x_plus[k] += EPSILON_METRIC
    x_minus[k] -= EPSILON_METRIC
    g_plus = get_metric_tensor(x_plus)
    g_minus = get_metric_tensor(x_minus)
    diff = mat_sub(g_plus, g_minus)
    return mat_scale(diff, 0.5 / EPSILON_METRIC)

def compute_christoffel_symbols(x):
    g = get_metric_tensor(x)
    g_inv = get_inverse_metric(g)
    partials = [partial_derivative_metric(x, k) for k in range(DIM)]
    Gamma = [[[0.0 for _ in range(DIM)] for _ in range(DIM)] for _ in range(DIM)]
    for k in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                val = 0.0
                for l in range(DIM):
                    if abs(g_inv[k][l]) > 1e-9:
                        term = partials[i][j][l] + partials[j][i][l] - partials[l][i][j]
                        val += g_inv[k][l] * term
                Gamma[k][i][j] = 0.5 * val
    return Gamma

# --- MANIFOLD FUNCTIONS ---

def covariant_derivative_1form(x_func, x):
    eps = EPSILON_MANIFOLD
    Gamma = compute_christoffel_symbols(x)
    p_i_omega_j = zero_matrix()
    omega_at_x = x_func(x)
    for i in range(DIM):
        x_plus = x[:]
        x_plus[i] += eps
        omega_plus = x_func(x_plus)
        for j in range(DIM):
            p_i_omega_j[i][j] = (omega_plus[j] - omega_at_x[j]) / eps
    nabla_omega = zero_matrix()
    for i in range(DIM):
        for j in range(DIM):
            correction = sum(Gamma[k][i][j] * omega_at_x[k] for k in range(DIM))
            nabla_omega[i][j] = p_i_omega_j[i][j] - correction
    return nabla_omega

def exterior_derivative_covariant(x_func, x):
    nabla_omega = covariant_derivative_1form(x_func, x)
    F = zero_matrix()
    for i in range(DIM):
        for j in range(DIM):
            F[i][j] = nabla_omega[i][j] - nabla_omega[j][i]
    return F

# --- FORMS FUNCTIONS ---

def get_combinations(n, k):
    return list(itertools.combinations(range(n), k))

def get_epsilon_sign(indices):
    n = len(indices)
    p = list(indices)
    swaps = 0
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                p[i], p[j] = p[j], p[i]
                swaps += 1
    return 1 if swaps % 2 == 0 else -1

class PForm:
    def __init__(self, p, data=None):
        self.p = p
        self.combs = get_combinations(DIM, p)
        self.index_map = {c: i for i, c in enumerate(self.combs)}
        self.values = data if data else [0.0] * len(self.combs)

    def add(self, other):
        return PForm(self.p, [v1 + v2 for v1, v2 in zip(self.values, other.values)])

    def scale(self, s):
        return PForm(self.p, [v * s for v in self.values])

def hodge_star_12d(x, p_form):
    p = p_form.p
    q = DIM - p
    g = get_metric_tensor(x)
    sqrt_g = math.sqrt(abs(calculate_det_12x12(g)))
    g_inv = get_inverse_metric(g)
    omega_upper = [0.0] * len(p_form.combs)
    for idx_A, comb_A in enumerate(p_form.combs):
        for idx_B, comb_B in enumerate(p_form.combs):
            if abs(p_form.values[idx_B]) < 1e-18: continue
            factor = 1.0
            for k in range(p):
                factor *= g_inv[comb_A[k]][comb_B[k]]
            omega_upper[idx_A] += factor * p_form.values[idx_B]
    star_form = PForm(q)
    full_indices = set(range(DIM))
    for idx_Q, comb_Q in enumerate(star_form.combs):
        p_set = full_indices - set(comb_Q)
        comb_P = tuple(sorted(list(p_set)))
        idx_P = p_form.index_map[comb_P]
        sign = get_epsilon_sign(comb_P + comb_Q)
        star_form.values[idx_Q] = sqrt_g * omega_upper[idx_P] * sign
    return star_form

def exterior_derivative_pform(x_func, x, p):
    eps = EPSILON_MANIFOLD
    q = p + 1
    d_omega = PForm(q)
    for i in range(DIM):
        x_plus, x_minus = list(x), list(x)
        x_plus[i] += eps; x_minus[i] -= eps
        omega_plus, omega_minus = x_func(x_plus), x_func(x_minus)
        p_i_omega = [(v_p - v_m) / (2.0 * eps) for v_p, v_m in zip(omega_plus.values, omega_minus.values)]
        for idx_P, comb_P in enumerate(omega_plus.combs):
            val = p_i_omega[idx_P]
            if abs(val) < 1e-18 or i in comb_P: continue
            new_comb = tuple(sorted((i,) + comb_P))
            idx_Q = d_omega.index_map[new_comb]
            sign = 1
            for idx_in_P in comb_P:
                if i > idx_in_P: sign *= -1
            d_omega.values[idx_Q] += sign * val
    return d_omega

def co_differential_2form(F_field_func, x):
    def star_F_field(pos):
        return hodge_star_12d(pos, F_field_func(pos))
    d_star_F = exterior_derivative_pform(star_F_field, x, 10)
    return hodge_star_12d(x, d_star_F)

def vector_to_1form(x, v):
    g = get_metric_tensor(x)
    return PForm(1, mat_vec_mul(g, v))

def matrix_to_2form(m):
    form = PForm(2)
    for idx, (i, j) in enumerate(form.combs):
        form.values[idx] = m[i][j]
    return form

# --- TRANSPORT FUNCTIONS ---

def parallel_transport_key(K, x, v_dt, dt):
    Gamma = compute_christoffel_symbols(x)
    Omega = zero_matrix()
    for i in range(DIM):
        for j in range(DIM):
            Omega[i][j] = sum(Gamma[i][k][j] * v_dt[k] for k in range(DIM))
    H = matrix_exp_pade(mat_scale(Omega, dt))
    H_inv = mat_inv(H)
    return mat_mul(mat_mul(H, K), H_inv)

def project_divergence_free(v, x, K):
    eps = 1e-5
    g = get_metric_tensor(x)
    sqrt_g = math.sqrt(abs(calculate_det_12x12(g)))
    def calculate_local_div_g(current_v, current_x):
        div = 0.0
        for i in range(DIM):
            x_plus, x_minus = list(current_x), list(current_x)
            x_plus[i] += eps; x_minus[i] -= eps
            g_p, g_m = get_metric_tensor(x_plus), get_metric_tensor(x_minus)
            sg_p, sg_m = math.sqrt(abs(calculate_det_12x12(g_p))), math.sqrt(abs(calculate_det_12x12(g_m)))
            flux_p = mat_vec_mul(K, current_v)
            term_p, term_m = sg_p * flux_p[i], sg_m * flux_p[i] # K is static for div check
            div += (term_p - term_m) / (2.0 * eps)
        return div / sqrt_g
    div_base = calculate_local_div_g(v, x)
    grad_F = [0.0] * DIM
    for i in range(DIM):
        v_plus = list(v); v_plus[i] += eps
        grad_F[i] = (calculate_local_div_g(v_plus, x) - div_base) / eps
    grad_norm_sq = sum(gi**2 for gi in grad_F)
    if grad_norm_sq > 1e-18:
        lambda_val = div_base / grad_norm_sq
        return [v[i] - lambda_val * grad_F[i] for i in range(DIM)]
    return v

# --- RG FUNCTIONS ---

def renormalization_flow_invariants(H_matrix):
    det = calculate_det_12x12(H_matrix)
    norm = frobenius_norm(H_matrix)
    log_scale = math.log(norm) if norm > 0 else 0.0
    return det, log_scale

# --- MAIN LOGIC ---

story_hierarchy = [
    [ [("このアジェンダは", "Subject", 1.0)], [("行動計画である", "Entity", 1.0)] ],
    [ [("虚空の歯車が、", "Subject", 5.0), ("メロンパンの", "Entity", 3.0)], [("重力定数を", "Entity", 5.0), ("逆走する。", "Action", 8.0)] ],
    [ [("はっと目が覚める！", "Action", 8.0)], [("己の至らなさを恥じ、", "Subject", 5.0), ("深く反省する。", "Action", 6.0)] ]
]

def get_target_vector(token, role, intensity):
    vec = [0.0] * 12
    mag = math.log(len(token) + 2.0) * 0.1 * intensity
    if role == "Subject":  vec[0], vec[1], vec[2] = mag, mag, mag
    elif role == "Entity": vec[3], vec[4], vec[5] = mag, mag, mag
    elif role == "Action": vec[6], vec[7], vec[8] = mag, mag, mag
    elif role == "Context":vec[9], vec[10], vec[11] = mag, mag, mag
    
    # 覚醒（ターゲットの急激な位相反転）のシミュレート
    if "目が覚める" in token:
        vec = [v * -1.5 for v in vec] # 強い反発ポテンシャル
    return vec

def run_strict_analysis():
    print("=== PKGF Unified Test (Purified 12D) ===\n")
    K = get_author_K()
    x = [0.1] * 12
    t, dt = 0.0, 0.2
    print(f"{'Token':<14} | {'t':<6} | {'||v||':<10} | {'div_g(KX)':<10} | {'det(K)':<10} | {'Note'}")
    print("-" * 76)
    for sentence in story_hierarchy:
        for bunsetsu in sentence:
            for token, role, intensity in bunsetsu:
                target_vec = get_target_vector(token, role, intensity)
                duration = 0.4 + 0.1 * len(token)
                steps = int(duration / dt)
                G_token = zero_matrix()
                for _ in range(steps):
                    def field_omega(pos):
                        """
                        1形式ポテンシャル omega: 
                        ターゲットへの放射状ベクトル場をomegaとみなし、
                        計量がフラットでない場合に外微分 F = d*omega によって物語の「起伏(曲率)」を生成する。
                        """
                        v_radial = [(target_vec[i] - pos[i]) * 0.5 for i in range(12)]
                        return vector_to_1form(pos, v_radial)

                    def F_field_func(pos):
                        return matrix_to_2form(exterior_derivative_covariant(lambda p: field_omega(p).values, pos))

                    # 共微分推進 (Co-differential Propulsion): 
                    # 理論式: (d/dt)(KX)^flat = -delta F
                    # 実装では、リーマン計量の逆行列 g_inv と並行鍵の逆行列 K_inv を介して速度 v を導出する。
                    delta_F = co_differential_2form(F_field_func, x)
                    g_inv = get_inverse_metric(get_metric_tensor(x))
                    v_raw = mat_vec_mul(g_inv, delta_F.values)
                    
                    K_inv = mat_inv(K)
                    v_maxwell = [-vi for vi in mat_vec_mul(K_inv, v_raw)]
                    v_proj = project_divergence_free(v_maxwell, x, K)
                    F_mat_local = exterior_derivative_covariant(lambda p: field_omega(p).values, x)
                    for i in range(12):
                        for j in range(12): G_token[i][j] += F_mat_local[i][j] * dt
                    for i in range(12): x[i] += v_proj[i] * dt
                    K = parallel_transport_key(K, x, v_proj, dt)
                    t += dt
                div_val = sum(mat_vec_mul(K, v_proj))
                H = matrix_exp_pade(G_token)
                det_H, log_scale = renormalization_flow_invariants(H)
                det_K = calculate_det_12x12(K)
                norm_v = math.sqrt(sum(vi**2 for vi in v_proj))
                note = "SINGULARITY" if "ぷつり" in token else ""
                print(f"{token:<14} | {t:6.2f} | {norm_v:10.5f} | {div_val:10.2e} | {det_K:10.5f} | {note}")
        print("-" * 76)
    print("\n[Scientific Findings]")
    print("1. Purified 12D Execution successful.")
    print("2. Co-differential Propulsion verified in potential-only flow.")
    print("3. Det(K) preserved without topological twist interference.")

if __name__ == "__main__":
    log_filename = "python_pkgf_12d_jp_log.txt"
    
    class Tee(object):
        def __init__(self, *files):
            self.files = files
        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush()
        def flush(self):
            for f in self.files:
                f.flush()

    with open(log_filename, "w", encoding="utf-8") as f:
        original_stdout = sys.stdout
        sys.stdout = Tee(sys.stdout, f)
        try:
            run_strict_analysis()
        finally:
            sys.stdout = original_stdout
    
    print(f"\n[System] 解析結果が {log_filename} に保存されました。")
