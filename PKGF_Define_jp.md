# 並行鍵幾何流 (Parallel Key Geometric Flow)

全ての実験リソースはこのレポジトリで公開されています：https://github.com/aikenkyu001/PKGF

## 一. 幾何的舞台 (Geometric Stage)
- **次元数**: $N = 12$。接束 $TM$ は以下の4セクターに直交分解される：
  \[ TM = T_{Subject}M \oplus T_{Entity}M \oplus T_{Action}M \oplus T_{Context}M \]
- **文脈依存計量 (Contextual Warping)**:
  計量テンソル $g$ はフラットではなく、Contextセクターの座標強度によって動的に歪む：
  \[ g_{ii}(x) = 1.0 + 0.5 \tanh(x_{context}) \]
  これにより、物語の背景（Context）が「場」の広がりや密度を決定する。

## 二. 並行鍵 (The Parallel Key) $K$
- **定義**: $K \in \Gamma(\mathrm{End}(TM))$ は、多様体上の論理構造を定義する $(1,1)$ テンソル。
- **並行輸送条件**: 理論的には $\nabla K = 0$。実装上は、流れ $v$ に沿った**随伴ホロノミー更新**によってこれを実現する：
  \[ K(t+dt) = H K(t) H^{-1}, \quad H = \exp(\Omega dt) \]
  ここで $\Omega$ はレヴィ＝チヴィタ接続 $\Gamma^i_{kj} v^k$ から導かれる接続行列。これにより、著者の論理軸（$\det(K)$）は、流れの中でも代数的に保存される。

## 三. 基礎方程式系 (Fundamental Equations)

### 1. 共微分推進 (Co-differential Propulsion)
意味の流動 $v$ は、目標引力から生じる1形式ポテンシャル $\omega$ の「渦」である2形式 $F = d\omega$ の**共微分 (co-differential)** によって推進される：
\[ \frac{\partial}{\partial t}(KX)^\flat = -\delta F = -\star d \star F \]
これは、マクスウェル方程式の真空解における電磁流力学の拡張であり、意味の流束 $KX$ の時間変化が幾何的な「力の源」と釣り合うことを示す。

### 2. 発散自由条件 (Divergence-free Constraint)
論理的一貫性を保つため、流束 $KX$ は常にソースフリー（発散ゼロ）に保たれる：
\[ \operatorname{div}_g (KX) = 0 \]
実装では、メトリック重み付きのヤコビアンを用いて速度ベクトル $v$ を射影することで、この条件を担保する。

## 四. 非可換ホロノミーと物語の収束
- **ホロノミー生成子**: 各トークン通過時に生じる曲率 $F$ の積分を生成子 $G$ とし、その指数写像 $H = \exp(G)$ を物語の「意味の変換」と定義する。
- **物語の収束性**: 生成子 $G$ の Frobenius ノルムは、物語の劇的な転換点（特異点）におけるエネルギー密度を表現し、物語が目標ポテンシャル $\omega$ に向かって正しく収束しているかを評価する。

## 五. 科学的保存則
- **情報の保存**: 並行鍵 $K$ が随伴変換を受けるため、その固有値（論理の重み）の積 $\det(K)$ は定数となる。
- **エネルギー等分配**: 推進力 $-\delta F$ と計量 $g$ の相互作用により、意味の運動エネルギー $\frac{1}{2}g(v,v)$ は文脈に応じて最適化される。
