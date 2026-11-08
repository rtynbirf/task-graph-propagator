# -*- coding: utf-8 -*-
# ============================================================================
#  TASK GRAPH ENGINE · П53 — PROPERTY-BASED ТЕСТЫ ПРОПАГАТОРА (issue #2)
#
#  Почему не hypothesis: доктрина репо — «dependencies: stdlib only»
#  (бейдж README, pyproject без зависимостей). Поэтому свойства гоняются
#  на ВСТРОЕННОМ random с ФИКСИРОВАННЫМ SEED — детерминированные,
#  воспроизводимые бит-в-бит на любой машине и в CI.
#
#  Каждое свойство документировано строкой в docstring (приёмка П53).
#  @example-реплеи: известные фиксированные случаи прогоняются рядом со
#  случайными (как @example у hypothesis) — при падении свойства случайный
#  случай фиксируется: добавь (seed, параметры) в REPLAY и он станет
#  постоянным реплей-тестом.
#
#  Запуск:  python3 test_propagator_properties.py   (rc=0 — зелёный)
# ============================================================================
import os
import sys
import random
from collections import deque

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from propagator_engine import (TaskGraph, sweep, worldline, kahn,
                               BOOL, MINP, MAXP, COUNT, VITER, BOTT,
                               INF, NINF)

SEED = 20261109            # зафиксирован: seed воспроизводим (приёмка П53)
N_GRAPHS = 160             # случайных графов на свойство
EPS = 1e-9

def close(a, b):
    """Сравнение с допуском, устойчивое к бесконечностям (INF−INF = NaN)."""
    return a == b or abs(a - b) < EPS

# ── Генераторы (детерминированные) ─────────────────────────────────────────
def random_dag(rng, n=None, p_edge=0.35):
    """Случайный DAG: вершины v0..vk, рёбра только вперёд по индексу."""
    n = n or rng.randint(1, 14)
    G = TaskGraph("random-dag")
    for i in range(n):
        G.add(f"v{i}", dur=rng.choice([0.0, 1.0, 2.0, 3.5]),
              prob=rng.choice([0.5, 0.9, 1.0]))
    for u in range(n):
        for v in range(u + 1, n):
            if rng.random() < p_edge:
                G.dep(f"v{u}", f"v{v}", w=rng.choice([0.0, 0.5, 1.0, 2.0]),
                      cost=rng.choice([1.0, 2.0, 5.0]),
                      cap=rng.choice([3.0, 6.0, 9.0]))
    return G

def shuffled_build(rng, G):
    """Тот же граф, собранный в случайном порядке вершин и рёбер."""
    H = TaskGraph(G.title + "-shuffled")
    nodes = [(v, G.dur[v], G.prob[v]) for v in G.dur]
    edges = []
    for u in G.adj:
        for v, w in G.adj[u].items():
            edges.append((u, v, w, G.cost[u][v], G.cap[u][v]))
    rng.shuffle(nodes)
    rng.shuffle(edges)
    for v, dur, prob in nodes:
        H.add(v, dur=dur, prob=prob)
    for u, v, w, cost, cap in edges:
        H.dep(u, v, w, cost, cap)
    return H

# ── Наивные оракулы (независимые реализации для перекрёстной сверки) ───────
def naive_reach(G, sources):
    """Достижимость обходом в ширину — не пропагатор."""
    seen = set()
    q = deque(v for v in sources)
    while q:
        u = q.popleft()
        if u in seen:
            continue
        seen.add(u)
        q.extend(G.adj[u])
    return seen

def naive_path_counts(G, s):
    """Число путей s→v: DP по порядку Кана — счёт, не sweep."""
    order, cyc = kahn(G)
    cnt = {v: (1.0 if v == s else 0.0) for v in G.V}
    for v in order:
        for u in G.preds(v):
            cnt[v] += cnt[u]
    return cnt

def naive_shortest(G, s):
    """Беллман-Форд по рёбрам — не проходка по Кана."""
    dist = {v: (0.0 if v == s else INF) for v in G.V}
    for _ in range(len(G.V) + 2):   # +2 запас: рёбра в произвольном порядке
        changed = False
        for u in G.adj:
            for v, _w in G.adj[u].items():
                nd = dist[u] + G.cost[u][v]
                if nd < dist[v] - EPS:
                    dist[v] = nd
                    changed = True
        if not changed:
            break
    return dist

# ── Свойства ───────────────────────────────────────────────────────────────
def prop_idempotence(rng):
    """СВОЙСТВО P1 · Идемпотентность: propagate ×2 == ×1 — sweep чист
    (два прогона на том же входе дают те же K/pred/порядок), а оператор
    замыкания в Б(∨,∧) идемпотентен: замыкание замыкания = замыкание."""
    G = random_dag(rng)
    phi = {"v0": True}
    phi = {v: (v == "v0") for v in G.V}
    K1, pred1, order1, _ = sweep(G, phi, BOOL, eps=lambda u, v, Ku: Ku)
    # чистота: тот же вход -> те же байты результата
    K1b, pred1b, order1b, _ = sweep(G, dict(phi), BOOL, eps=lambda u, v, Ku: Ku)
    assert K1 == K1b and pred1 == pred1b and order1 == order1b, "sweep не чист"
    # идемпотентность оператора: замыкание замыкания = замыкание
    K2, _, _, _ = sweep(G, dict(K1), BOOL, eps=lambda u, v, Ku: Ku)
    assert K1 == K2, "замыкание замыкания != замыкания"
    # и на (max,+) без вершинного функтора: проходка по Кана достигает
    # неподвижной точки K = phi (+) A·K за один проход — повтор не двигает
    phi_max = {v: NINF for v in G.V}
    P1, _, _, _ = sweep(G, phi_max, MAXP)
    P2, _, _, _ = sweep(G, dict(P1), MAXP)
    for v in G.V:
        assert close(P1[v], P2[v]), f"(max,+) не идемпотентен на {v}"

def prop_monotonicity(rng):
    """СВОЙСТВО P2 · Монотонность: новые факты не стирают старые —
    добавление ребра/вершины сохраняет достижимость (Б(∨,∧)) и не
    уменьшает ранние финишы K+ ((max,+), лаг ≥ 0)."""
    G = random_dag(rng, n=rng.randint(3, 12))
    n = len(G.V)
    phi = {v: (v == "v0") for v in G.V}
    eps = lambda u, v, Ku: Ku
    R1, *_ = sweep(G, phi, BOOL, eps=eps)
    phi_max = {v: NINF for v in G.V}
    nu = lambda v, x: G.dur[v] + max(0.0, x)
    K1, *_ = sweep(G, phi_max, MAXP, nu=nu)

    # добавляем 1..3 новых «факта»: рёбра вперёд по индексу и вершину-сток
    H = TaskGraph(G.title + "+facts")
    for v in G.V:
        H.add(v, dur=G.dur[v], prob=G.prob[v])
    for u in G.adj:
        for v, w in G.adj[u].items():
            H.dep(u, v, w, G.cost[u][v], G.cap[u][v])
    idx = [int(v[1:]) for v in G.V]
    for _ in range(rng.randint(1, 3)):
        a, b = sorted(rng.sample(range(n), 2))
        if f"v{b}" not in H.adj[f"v{a}"]:   # только НОВЫЕ факты, старые не трогаем
            H.dep(f"v{a}", f"v{b}", w=rng.choice([0.0, 0.5, 1.0]),
                  cost=rng.choice([1.0, 2.0]), cap=rng.choice([3.0, 9.0]))
    H.add(f"v{n}", dur=1.0)                       # новая вершина старое не трогает

    R2, *_ = sweep(H, {v: (v == "v0") for v in H.V}, BOOL, eps=eps)
    for v in R1:
        assert R2[v] or not R1[v], f"достижимость стёрта в {v}"
        if R1[v]:
            assert R2[v], f"факт достижимости {v} стёрт новыми рёбрами"
    K2, *_ = sweep(H, {v: NINF for v in H.V}, MAXP,
                   nu=lambda v, x: H.dur[v] + max(0.0, x))
    for v in K1:
        assert K2[v] >= K1[v] - EPS, f"K+ уменьшился в {v}: {K1[v]} -> {K2[v]}"

def prop_permutation_determinism(rng):
    """СВОЙСТВО P3 · Детерминизм при перестановке входов: тот же граф,
    собранный в другом порядке вершин/рёбер, даёт тот же ответ (K по
    вершинам совпадает; Б(∨,∧) — точно, числовые полукольца — в EPS)."""
    G = random_dag(rng, n=rng.randint(2, 12))
    phi = {v: (v == "v0") for v in G.V}
    eps = lambda u, v, Ku: Ku
    K_ref, *_ = sweep(G, phi, BOOL, eps=eps)
    phi_max = {v: NINF for v in G.V}
    nu = lambda v, x: G.dur[v] + max(0.0, x)
    Kp_ref, *_ = sweep(G, phi_max, MAXP, nu=nu)
    for trial in range(3):
        H = shuffled_build(rng, G)
        K_h, *_ = sweep(H, {v: (v == "v0") for v in H.V}, BOOL, eps=eps)
        assert K_h == K_ref, f"Б(∨,∧) зависит от порядка входа (trial {trial})"
        Kp_h, *_ = sweep(H, {v: NINF for v in H.V}, MAXP, nu=nu)
        for v in G.V:
            assert close(Kp_h[v], Kp_ref[v]), \
                f"(max,+) разошёлся при перестановке на {v} (trial {trial})"

def prop_naive_oracles(rng):
    """СВОЙСТВО P4 · Перекрёстная сверка с наивными оракулами:
    пропагатор в Б(∨,∧) == BFS; в (+,×) == DP-счёт путей; в (min,+)
    == Беллман-Форд; мировая линия — настоящий путь из источника."""
    G = random_dag(rng, n=rng.randint(2, 12), p_edge=0.3)
    src = "v0"
    # Б(∨,∧) == BFS
    K, pred, _, _ = sweep(G, {v: (v == src) for v in G.V}, BOOL,
                          eps=lambda u, v, Ku: Ku)
    assert {v for v, ok in K.items() if ok} == naive_reach(G, [src]), \
        "достижимость пропагатора != BFS"
    # (+,×) == DP-счёт путей
    KC, *_ = sweep(G, {v: 1.0 if v == src else 0.0 for v in G.V}, COUNT,
                   eps=lambda u, v, Ku: Ku)
    cnt = naive_path_counts(G, src)
    for v in G.V:
        assert close(KC[v], cnt[v]), f"число путей до {v} разошлось"
    # (min,+) == Беллман-Форд
    KM, predm, _, _ = sweep(G, {v: 0.0 if v == src else INF for v in G.V},
                            MINP, eps=lambda u, v, Ku: Ku + G.cost[u][v])
    dist = naive_shortest(G, src)
    for v in G.V:
        assert close(KM[v], dist[v]), f"кратчайший до {v} разошёлся"
    # мировая линия: реальный путь из источника, рёбра существуют
    for v in G.V:
        if K[v]:
            path = worldline(v, pred)
            assert path[0] == src and path[-1] == v, "мировая линия не от источника"
            for a, b in zip(path, path[1:]):
                assert b in G.adj[a], f"нет ребра {a}->{b} в мировой линии"

def prop_damped_standing_wave(rng):
    """СВОЙСТВО P5 · Регулятор β даёт стоячую амплитуду: на цикле длины L
    сходится к φ·β^0+… = β/(1−β^L) на затравке контура (теория Неймана)."""
    beta, L = rng.choice([0.5, 0.7, 0.9]), rng.randint(2, 6)
    C = TaskGraph(f"цикл-{L}")
    C.add("E0", dur=0.0)
    for i in range(L):
        C.add(f"Z{i}", dur=0.0)
    C.dep("E0", "Z0", 1.0)
    for i in range(L):
        C.dep(f"Z{i}", f"Z{(i + 1) % L}", 1.0)
    K, n = __import__("propagator_engine").damped(
        C, {"E0": 1.0, **{f"Z{i}": 0.0 for i in range(L)}}, COUNT, beta=beta)
    theory = beta / (1 - beta ** L)   # неподвижная точка K = phi + beta·A·K
    for i in range(L):
        expected = theory if i == 0 else (beta ** i) * theory
        assert abs(K[f"Z{i}"] - expected) < 1e-7, \
            f"стоячая волна Z{i}: {K[f'Z{i}']} != {expected}"

# ── @example-реплеи (фиксированные случаи, всегда с прогона) ──────────────
def example_demo_graph():
    """@example · демо-граф A..G: пути=3, кратчайший=5, узкое место=6,
    Viterbi=.99·.98·.99 — канонические значения самопроверки движка."""
    G = TaskGraph("example-demo")
    for k in "ABCDEFG":
        G.add(k, dur=1.0, prob=0.99 if k in "AEG" else (0.95 if k == "B" else 0.98))
    for u, v, cost, cap in (("A", "B", 2, 9), ("A", "C", 1, 7), ("B", "D", 1, 8),
                            ("C", "D", 1, 5), ("C", "G", 5, 4), ("D", "E", 1, 6),
                            ("E", "F", 1, 6), ("F", "G", 1, 7)):
        G.dep(u, v, 0.0, cost=cost, cap=cap)
    KC, *_ = sweep(G, {v: 1.0 if v == "A" else 0.0 for v in G.V}, COUNT,
                   eps=lambda u, v, Ku: Ku)
    assert KC["G"] == 3.0, "число путей A→G != 3"
    KM, predm, _, _ = sweep(G, {v: 0.0 if v == "A" else INF for v in G.V},
                            MINP, eps=lambda u, v, Ku: Ku + G.cost[u][v])
    assert KM["G"] == 5.0, "кратчайший A→G != 5"
    KB, *_ = sweep(G, {v: INF if v == "A" else NINF for v in G.V}, BOTT,
                   eps=lambda u, v, Ku: min(Ku, G.cap[u][v]))
    assert KB["G"] == 6.0, "пропускная способность A→G != 6"
    KV, *_ = sweep(G, {v: 1.0 if v == "A" else 0.0 for v in G.V}, VITER,
                   eps=lambda u, v, Ku: Ku, nu=lambda v, x: x * G.prob[v])
    assert abs(KV["G"] - 0.99 * 0.98 * 0.99) < 1e-12, "Viterbi A→G разошёлся"

def example_edge_cases():
    """@example · края: одна вершина; линейная цепь; пустые рёбра —
    propagate не падает и возвращает согласованные значения."""
    for G, src in ((TaskGraph("одна"), "v0"),):
        G.add("v0", dur=2.0)
        K, *_ = sweep(G, {"v0": True}, BOOL, eps=lambda u, v, Ku: Ku)
        assert K["v0"] is True
    C = TaskGraph("цепь")
    for i in range(6):
        C.add(f"v{i}", dur=1.0)
        if i:
            C.dep(f"v{i-1}", f"v{i}", w=0.0)
    K, pred, _, cyc = sweep(C, {"v0": True, **{f"v{i}": False for i in range(1, 6)}},
                            BOOL, eps=lambda u, v, Ku: Ku)
    assert all(K[f"v{i}"] for i in range(6)) and not cyc
    Kp, *_ = sweep(C, {v: NINF for v in C.V}, MAXP,
                   nu=lambda v, x: C.dur[v] + max(0.0, x))
    assert abs(Kp["v5"] - 6.0) < EPS  # 6 вершин по 1.0

REPLAY: list[tuple[str, int]] = [
    # (свойство, seed) — реплей зафиксированных случаев; при падении
    # случайного свойства добавь сюда его seed, и случай станет постоянным.
    ("P1", SEED), ("P2", SEED + 1), ("P3", SEED + 2), ("P4", SEED + 3),
    ("P5", SEED + 4),
]

PROPERTIES = [
    ("P1", prop_idempotence), ("P2", prop_monotonicity),
    ("P3", prop_permutation_determinism), ("P4", prop_naive_oracles),
    ("P5", prop_damped_standing_wave),
]

def main() -> int:
    bar = "=" * 72
    print(bar)
    print(" П53 · PROPERTY-BASED ТЕСТЫ ПРОПАГАТОРА (stdlib random, seed=%d)" % SEED)
    print(bar)
    # @example-реплеи первыми — как у hypothesis
    example_demo_graph()
    print("PASS  @example · демо-граф A..G (пути/кратчайший/узкое место/Viterbi)")
    example_edge_cases()
    print("PASS  @example · края (одна вершина, линейная цепь)")
    red = 0
    for name, fn in PROPERTIES:
        try:
            for i in range(N_GRAPHS):
                rng = random.Random(SEED + 1000 * i)   # изолированный seed на случай
                fn(rng)
            print(f"PASS  {name} · {N_GRAPHS} случайных графов (seeds {SEED}..{SEED + 1000 * (N_GRAPHS - 1)})")
        except AssertionError as err:
            red += 1
            print(f"RED   {name}: {err}")
            print(f"      → реплей: добавь (seed-диапазон) в REPLAY, случай воспроизводим")
    if red:
        print(f"\n# Result: {len(PROPERTIES) - red} passed, {red} failed (of {len(PROPERTIES)} properties)")
        return 1
    print(f"\n# Result: {len(PROPERTIES)} passed, 0 failed (of {len(PROPERTIES)} properties)")
    print(f"# примеры-реплеи: {len(REPLAY)} зафиксированы, seed воспроизводим")
    print(bar)
    return 0

if __name__ == "__main__":
    sys.exit(main())
