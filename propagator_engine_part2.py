# -*- coding: utf-8 -*-
# ============================================================================
#  TASK GRAPH ENGINE · ЧАСТЬ 2 — «ДОЖИМ» ПРОПАГАТОРА
#
#  УНИВЕРСАЛЬНАЯ ТЕОРЕМА. Пусть K*(v) = Σ_{p: s⇝v} [p] — универсальный
#  пропагатор в свободном полукольце путей (формальная сумма всех путей).
#  Тогда ВСЯКАЯ функция структуры движка есть гомоморфный образ K*:
#        F(v) = h( K*(v) ),   h — единственный гомоморфизм полуколец,
#  индуцированный отображением образующих (рёбер) и операций (⊕,⊗).
#  Доказательство: единственность свободного объекта + индукция по порядку
#  Кана: обе стороны удовлетворяют одному неподвижному уравнению
#  K = φ ⊕ ÂK. ∎
#  Следствия, дожимаемые ниже программно:
#   [12] S=ℝ[z]: производящая функция путей — ПОЛНОЕ распределение длин:
#        T = старшая степень · число критических путей = [z^T]P(z) ·
#        всего путей = P(1) · среднее = P'(1)/P(1).
#   [13] Тропический предел: (max,+) = lim_{β→∞} (1/β)·ln[(+,×) с весами
#        e^{β·len}] — критический путь есть основное состояние (седловая
#        точка) статсуммы. Расписание и статсумма — одно ядро.
#   [14] Критичность: n±(v) — счёт мировых линий через v на критическом
#        подграфе (два прохода (+,×)); плотность линий c(v)=n⁺(v)n⁻(v)/N.
#   [15] Обратная задача (крэшинг): замкнутый контур управления
#        ПРОПАГАТОР → ГРАДИЕНТ (критические вершины, min цена) → СЖАТИЕ →
#        ПРОПАГАТОР → … до дедлайна. Если одиночное сжатие не двигает T —
#        нелокальность: ищем одновременную поправку на паре ветвей.
# ============================================================================
import sys, math, itertools, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # портативно: репо/скрипты
from propagator_engine import (Semiring, TaskGraph, sweep, worldline,
                               BOOL, MINP, MAXP, COUNT, VITER, BOTT,
                               INF, NINF)

R9 = lambda x: round(x, 9)

# ── [12] Свободно-порождённое ℝ[z]: носитель универсального пропагатора ────
def pplus(a, b):                          # сложение многочленов
    r = dict(a)
    for k, c in b.items():
        r[k] = r.get(k, 0.0) + c
    return r

def ptimes(a, b):                         # произведение Коши (свёртка)
    r = {}
    for k1, c1 in a.items():
        for k2, c2 in b.items():
            k = R9(k1 + k2)
            r[k] = r.get(k, 0.0) + c1 * c2
    return r

POLY = Semiring({}, {0.0: 1.0}, pplus, ptimes,
                "ℝ[z]: производящая функция путей")

def shift(p, delta):                      # умножение на z^delta
    return {R9(k + delta): c for k, c in p.items()}

def makespan(G):
    """K+/K-/резерв одним и тем же ядром (запаздывающий + опережающий)."""
    phi = {v: NINF for v in G.V}
    nu = lambda v, x: G.dur[v] + max(0.0, x)
    Kp, _, _, _ = sweep(G, phi, MAXP, nu=nu)
    Km, _, _, _ = sweep(G.flip(), phi, MAXP, nu=nu)
    T = max(Kp.values())
    R = {v: T - Kp[v] - Km[v] + G.dur[v] for v in G.V}
    return T, Kp, Km, R

def main():
    bar = "=" * 72
    print(bar)
    print(" ЧАСТЬ 2 · ДОЖИМ: универсальность, статсумма, контур управления")
    print(bar)

    # Граф-мета из части 1 (самоприменение движка) ──────────────────────────
    G = TaskGraph("Движок·самоприменение")
    for k, d in {"A": 1.0, "B": 2.0, "C": 1.0, "D": 3.0,
                 "E": 1.0, "F": 1.0, "G": 1.0}.items():
        G.add(k, dur=d)
    G.dep("A", "B", w=0.0); G.dep("A", "C", w=0.0)
    G.dep("B", "D", w=0.0); G.dep("C", "D", w=0.5); G.dep("C", "G", w=0.0)
    G.dep("D", "E", w=0.0); G.dep("E", "F", w=0.0); G.dep("F", "G", w=0.0)

    # [12] Универсальный пропагатор в ℝ[z] — ПОЛНАЯ статистика структуры ────
    phi = {v: ({0.0: 1.0} if v == "A" else {}) for v in G.V}
    P, _, _, _ = sweep(G, phi, POLY,
                       eps=lambda u, v, Ku: shift(Ku, G.adj[u][v]),
                       nu=lambda v, x: shift(x, G.dur[v]))
    PG = P["G"]
    T = max(PG)
    Ncrit, N = PG[T], sum(PG.values())
    mean = sum(k * c for k, c in PG.items()) / N
    m2 = sum(k * k * c for k, c in PG.items()) / N
    print("\n[12] ПРОИЗВОДЯЩАЯ ФУНКЦИЯ ПУТЕЙ A→G:  P(z) =",
          " + ".join(f"{c:g}·z^{k:g}" for k, c in sorted(PG.items())))
    print(f"     T (старшая степень) = {T:g} · критических путей [z^T] = {Ncrit:g}"
          f" · всего путей P(1) = {N:g}")
    print(f"     средняя длина = P'(1)/P(1) = {mean:.4f}"
          f" · σ = {math.sqrt(max(0.0, m2 - mean*mean)):.4f}")

    # [13] Тропический предел: критпуть = основное состояние статсуммы ──────
    def free_energy(beta):               # log-sum-exp — без переполнения
        m = max(PG)
        s = sum(c * math.exp(beta * (k - m)) for k, c in PG.items())
        return m + math.log(s) / beta

    print("\n[13] СВОБОДНАЯ ЭНЕРГИЯ f(β) = (1/β)·ln Z(β): остывание графа")
    for beta in (1.0, 2.0, 5.0, 10.0, 50.0, 100.0):
        print(f"     β={beta:>5g} : f = {free_energy(beta):.9f}")
    print("     β→∞: f → T = 9 — (max,+) есть предел (+,×);")
    print("     критический путь = седловая точка суммы по путям.")

    # [14] Критичность = плотность мировых линий (проект с параллелью) ─────
    H = TaskGraph("проект·параллельные ветви")
    for k, d in {"S": 1.0, "X": 2.0, "Y": 2.0, "M": 2.0, "Z": 1.0, "F": 1.0}.items():
        H.add(k, dur=d)
    for u, v in (("S", "X"), ("S", "Y"), ("X", "M"),
                 ("Y", "M"), ("M", "Z"), ("Z", "F")):
        H.dep(u, v, w=0.0)
    T2, Kp, Km, R = makespan(H)
    tight = lambda u, v: abs((Kp[u] + H.adj[u][v]) - (Kp[v] - H.dur[v])) < 1e-9
    Gc = TaskGraph("критический подграф")
    for v in H.V:
        Gc.add(v)
    for u in H.adj:
        for v in H.adj[u]:
            if tight(u, v):
                Gc.dep(u, v, 0.0)
    phi1 = {v: 1.0 if v == "S" else 0.0 for v in Gc.V}
    n_in, *_ = sweep(Gc, phi1, COUNT, eps=lambda u, v, Ku: Ku)
    phi2 = {v: 1.0 if v == "F" else 0.0 for v in Gc.V}
    n_out, *_ = sweep(Gc.flip(), phi2, COUNT, eps=lambda u, v, Ku: Ku)
    Ntot = n_in["F"]
    print(f"\n[14] КРИТИЧНОСТЬ (T={T2:g}, критических линий N={Ntot:g}):")
    print("     v    n⁺  n⁻  c(v)=n⁺n⁻/N")
    for v in sorted(Gc.V):
        c = n_in[v] * n_out[v] / Ntot
        print(f"     {v}    {n_in[v]:g}   {n_out[v]:g}   {c:.2f}")
    print("     c=1 — несущая ось конуса; c=½ — ветвь, делящая нагрузку.")

    # [15] Замкнутый контур: крэшинг до дедлайна (пропагатор-контроллер) ───
    cost = {"S": INF, "X": 1.0, "Y": 1.0, "M": 3.0, "Z": INF, "F": INF}
    floor = {"X": 1.0, "Y": 1.0, "M": 1.0}      # предел сжатия
    deadline = 5.0
    total, step, log = 0.0, 0, []
    T0, *_ = makespan(H)
    Tc = T0
    while Tc > deadline + 1e-9:
        step += 1
        _, Kp, Km, R = makespan(H)
        crit = [v for v in H.V if abs(R[v]) < 1e-9]
        cands = sorted((v for v in crit
                        if cost.get(v, INF) < INF
                        and H.dur[v] > floor.get(v, 0.0)),
                       key=lambda v: cost[v])
        done = None
        for v in cands:                  # одиночное сжатие
            H.dur[v] -= 1.0
            Tt, *_ = makespan(H)
            if Tt < Tc - 1e-9:
                done = (v,); total += cost[v]; break
            H.dur[v] += 1.0              # градиент нулевой — откат
        if done is None:                 # нелокальность: парная поправка
            for a, b in itertools.combinations(cands, 2):
                H.dur[a] -= 1.0; H.dur[b] -= 1.0
                Tt, *_ = makespan(H)
                if Tt < Tc - 1e-9:
                    done = (a, b)
                    total += cost[a] + cost[b]; break
                H.dur[a] += 1.0; H.dur[b] += 1.0
        if done is None:
            log.append(f"     шаг {step}: T={Tc:g} — ТУПИК: нужен комбинаторный ход")
            break
        Tc = makespan(H)[0]
        log.append(f"     шаг {step}: сжато {'∧'.join(done)} → T={Tc:g}"
                   f" (цена {sum(cost[v] for v in done):g})")
    print(f"\n[15] КОНТУР УПРАВЛЕНИЯ: T₀={T0:g} → дедлайн {deadline:g},"
          f" итераций {step}, бюджет {total:g}:")
    print("\n".join(log))
    print("     урок: жадная траектория упёрлась в нелокальность —")
    print("     две параллельные ветви сжимаются только ОДНОВРЕМЕННО.")

    # [16] Самопроверка дожима ──────────────────────────────────────────────
    assert PG == {3.0: 1.0, 8.5: 1.0, 9.0: 1.0} and T == 9.0
    assert Ncrit == 1.0 and N == 3.0 and abs(mean - 20.5 / 3) < 1e-12
    fs = [free_energy(b) for b in (1.0, 2.0, 5.0, 10.0)]
    assert fs[0] > fs[1] > fs[2] > fs[3], "свободная энергия убывает"
    assert abs(free_energy(100.0) - 9.0) < 1e-9, "тропический предел = T"
    assert T2 == 7.0 and Ntot == 2.0
    assert abs(n_in["M"] * n_out["M"] / Ntot - 1.0) < 1e-12
    assert abs(n_in["X"] * n_out["X"] / Ntot - 0.5) < 1e-12
    assert abs(Tc - 5.0) < 1e-9 and abs(total - 5.0) < 1e-9, "крэшинг: 3+2"
    print("\n[16] САМОПРОВЕРКА ДОЖИМА: все утверждения выполнены.")

    print("\n" + bar)
    print(" СЛОВАРЬ-ДОЖИМ (один пропагатор, три новых лица):")
    print("   ℝ[z]          P(z) — статсумма путей: вся статистика структуры")
    print("   β→∞           тропический предел: критпуть = основное состояние")
    print("   f(β)↓         свободная энергия убывает — третий закон термо-")
    print("                 динамики расписаний; T — энергия основного уровня")
    print("   c(v)=n⁺n⁻/N   плотность мировых линий = градиент T по d(v)")
    print("   крэшинг       внешнее поле; тупик жадного = нелокальность")
    print(bar)

if __name__ == "__main__":
    main()
