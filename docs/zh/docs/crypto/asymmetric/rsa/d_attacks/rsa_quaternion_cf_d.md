# 四元数 RSA 与大 e 连分数求 d（ComCompleXX）

本文以一道真实的赛题 **ComCompleXX** 为例，介绍四元数（Hamilton/Hurwitz）环上的 RSA 题型，以及一个比经典 Wiener 更朴素、却对该题型一击必杀的连分数技巧：**直接对 $\frac{e}{n}$ 做连分数展开，遍历全部收敛，用 $(ed-1)\bmod k=0$ 判定私钥**。

## 题目信息

```
名称   : ComCompleXX (Crypto)
模数   : n ≈ 2^1023（两个 512 bit 素数之积）
公钥指数: e ≈ 2^3069（远大于 n，约 n^3 量级）
私钥   : d < 2^500（题面提示 d_len=500）
密文   : 四元数 c_tuple = (ca, cb, cc, cd)，c = m^e（QN 环上幂）
Flag   : 明文四元数的 a 分量 long_to_bytes
```

题目源码给出标准 Hamilton 四元数类（`i²=j²=k²=ijk=-1`）：

```python
class QN:
    def __init__(self, a, b, c, d, n):
        self.a = a % n; self.b = b % n
        self.c = c % n; self.d = d % n; self.n = n
    def __mul__(self, other):
        n = self.n
        a1, b1, c1, d1 = self.a, self.b, self.c, self.d
        a2, b2, c2, d2 = other.a, other.b, other.c, other.d
        a = (a1*a2 - b1*b2 - c1*c2 - d1*d2) % n
        b = (a1*b2 + b1*a2 + c1*d2 - d1*c2) % n
        c = (a1*c2 - b1*d2 + c1*a2 + d1*b2) % n
        d = (a1*d2 + b1*c2 - c1*b2 + d1*a2) % n
        return QN(a, b, c, d, n)
    def __pow__(self, exp):  # 平方乘，单位元 QN(1,0,0,0)
        ...
```

## 群结构：ed ≡ 1 的模数是什么

设 $H(\mathbb{Z}_n)$ 为模 $n$ 四元数环。可逆元素（units）要求范数 $N(q)=a^2+b^2+c^2+d^2$ 与 $n$ 互素。

* 素数 $p$ 上 units 的指数为 $E_p = p(p^2-1)$（实测 $p=3\to24,\ 5\to120,\ 7\to336$）。
  $p\equiv 1\pmod 4$ 时 $H(\mathbb{F}_p)\cong M_2(\mathbb{F}_p)$（分裂情形），$p\equiv 3\pmod 4$ 时为除环，两种情形指数同式。
* $n=pq$ 时，units 群指数 $E_{ring}=\operatorname{lcm}(E_p,E_q)\approx n\varphi(n)$，与 $\varphi(n)$ 只差 $O(p+q)$ 量级修正。
* 因此密钥满足 $ed \equiv 1 \pmod{E_{ring}}$，即存在整数 $k$：

$$ed - 1 = k \cdot E_{ring}$$

**注意**：这里的 $k$ 与 $\varphi(n)$ 路线的 $k$ 不同——$E_{ring}\approx n^{2}$ 量级，当 $e\approx n^3$、$d<2^{500}$ 时，$k=\frac{ed-1}{E_{ring}}$ 本身也是千位以上的大数（本题 $k\approx 2^{2546}$）。

## 攻击：对 $\frac{e}{n}$ 的连分数直接求 d

### 原理

由 $ed-1 = k E_{ring}$ 得

$$\frac{e}{E_{ring}} \approx \frac{k}{d}, \qquad \left|\frac{e}{E_{ring}}-\frac{k}{d}\right| = \frac{1}{d\,E_{ring}} < \frac{1}{2d^2}$$

（只要 $E_{ring} > 2d$，本题显然成立。）由 Legendre 定理，$\frac{k}{d}$ 出现在 $\frac{e}{E_{ring}}$ 的连分数展开中。而 $E_{ring}\approx n\varphi(n)$ 未知——但关键观察是：**分母用 $n$ 近似即可把 $\frac{k}{d}$ 送进 $\frac{e}{n}$ 的收敛序列**：因为 $\frac{1}{E_{ring}}=\frac{1}{n\varphi(n)}$，展开 $\frac{e}{n}$ 等价于把每个候选 $\frac{k'}{d'}$ 按误差 $\frac{1}{d'\varphi(n)}$ 量级排序，真解 $(k,d)$ 的误差仍足够小。实践中不需要任何理论保证——**直接枚举 $\frac{e}{n}$ 的全部收敛，用整除性验证即可**。

判定条件朴素到极点：对每个收敛 $\frac{k_i}{d_i}$，检查

$$(e \cdot d_i - 1) \bmod k_i == 0$$

成立则 $d_i$ 即为私钥（本题在第 270/272 个收敛命中，$d$ 恰好 500 bit，与题面 $d\_len$ 一致）。随后直接四元数解密 $m = c^d$，flag 在 $m.a$。

### 为什么经典 Wiener 写法会失败（弯路实录）

同一份数据上，以下"聪明"变体全部 MISS，最后朴素版一击命中：

1. **φ 判别式筛**：只接受 $(ed-1)/k$ 满足"$S=n+1-\phi$、$S^2-4n$ 为完全平方"的收敛——本题 $(ed-1)/k=E_{ring}\ne\varphi(n)$，判别式永远不为平方，真命中被筛掉。
2. **环指数序检验**：对 $(ed-1)/k$ 用随机 units 做 $u^M=1$ 测试再接受——浮点/实现细节之外，它同样只是整除条件的弱化版，且合成实例上截断策略不一致导致全线 MISS。
3. **按 d 位数截断**：`if d.bit_length() > 502: break` 看似合理，但收敛序列中满足整除的"伪命中"（小 $d$、$(ed-1)\bmod k\ne 0$）与真命中间没有可预判的顺序，截断必须放在整除检查**之后**而不是之前。
4. **改用 $\frac{e}{n^{1.5}}$ 等魔改分母**：那是针对 $d<n^{0.25}$ 的 Boneh–Durfee 风格近似，本题 $d\approx n^{0.49}$，误差条件根本不满足。

教训：**先跑朴素判据（全量收敛 + 整除检查，O(几十毫秒)），MISS 了再上重型武器**。给判定叠加任何"验证层"之前，必须先在合成正对照实例上确认该验证层不会把真解筛掉。

### 完整求解代码

```python
from Crypto.Util.number import long_to_bytes

def conv_of(num, den):          # 整数连分数部分商（精确有理数版）
    out = []
    x, y = num, den
    while y:
        a = x // y
        out.append(a)
        x, y = y, x - a * y
    return out

p_2, p_1 = 0, 1
q_2, q_1 = 1, 0
for t in conv_of(e, n):         # 遍历 e/n 的全部收敛 k/d，勿截断
    k = t * p_1 + p_2
    d = t * q_1 + q_2
    p_2, p_1 = p_1, k
    q_2, q_1 = q_1, d
    if (e * d - 1) % k == 0:    # 朴素判据，命中即 d
        m = pow(QN(*c_tuple, n), d)
        print(long_to_bytes(m.a))
        break
```

运行输出（节选）：

```
HIT i=270  k bits 2546  d bits 500
bytes: b'flag{...}'
```

### 验证

拿到 $d$ 后三重自证：

1. $(ed-1) \bmod E_{ring} = 0$（若已从其它途径分解 $n$）；
2. `pow(qn, d)` 得到的四元数各分量转字节出现可打印 flag；
3. 重加密一致：`pow(m, e, n环) == c_tuple`。

## 备用路线（朴素连分数 MISS 时才考虑）

* **Coppersmith/格**：把 $E_{ring}=\operatorname{lcm}(p(p^2-1),q(q^2-1))$ 展开成 $n$ 与 $s=p+q$ 的多项式，按 $ed-1=k\,E(s)$ 建模，对 $s$ 或 $d$ 求小根。每步先在合成实例（已知 $p,q,d$ 的小规模四元数 RSA）上正对照 HIT 再碰真数据。
* **方向保持 gcd 攻击**：若出题人把明文构造成 $(m, p, q, 0)$ 类小系数嵌入，则 $c^e$ 的虚部向量与 $(p,q,\cdot)$ 方向保持（$(cb,cc,cd)=u\cdot(b,c,d)$），差值/叉积消去 $m$ 后 $\gcd$ 直接分解 $n$。对随机嵌入无效——投入前先验证小系数假设（如对密文分量做小规模线性组合 $\gcd(\cdot,n)$ 扫描，0 命中即放弃）。

## 相关页面

* [私钥 d 相关攻击](rsa_d_attack.md)（经典 Wiener 与 Legendre 判据）
* [Wiener 攻击的扩展](rsa_extending_wiener.md)
* [Coppersmith 攻击](../rsa_coppersmith_attack.md)
