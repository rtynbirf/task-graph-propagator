# -*- coding: utf-8 -*-
# ============================================================================
#  TASK GRAPH ENGINE · ОБЪЕКТ=ПРЕДМЕТ «ПРОПАГАТОР»
#  ОДИН АЛГОРИТМ — ПРОПАГАТОРНАЯ РАЗВЁРТКА (Propagator Sweep)
#
#  ЯДРО:   K(v) = nu_v( phi(v)  (+)  SUM_{u->v} eps_uv( K(u) ) )   — в S
#
#  ВЫВОД (деривация):
#   1) Пропагатор = функция Грина графа задач: амплитуда дойти до v.
#      Сумма по ВСЕМ путям u⇝v, произведение переносов по пути
#      (фейнмановская сумма по путям):
#          K(v) = SUM_{p: u⇝v} PROD_{e in p} eps_e
#   2) Группировка путей по последнему ребру => K = phi + A·K,
#      т.е. K = (1 - A)^{-1}·phi — резольвента (уравнение Дайсона),
#      ряд Неймана по оператору смежности A.
#   3) На DAG ряд конечен => ОДНА проходка в порядке Кана = точное
#      решение. Порядок Кана — та же развёртка в булевом полукольце.
#   4) Каждая функция структуры = то же ядро в своём полукольце S:
#        B(∨,∧)   достижимость · транзитивное замыкание · K(v,v)=циклы
#        (max,+)  ранние финишы · расписание · критический путь
#        (max,+) на G^T  высоты H(v) · поздние сроки · резерв R=T-K+ -K-+d
#        (min,+)  кратчайший маршрут
#        (+,×)    статсумма (число путей)
#        (max,×)  надёжность (Viterbi)
#        (max,min) пропускная способность (узкое место)
#      Обратная проходка по G^T = опережающий пропагатор K-;
#      R(v) = T − K+(v) − K−(v) + d(v); R=0 <=> вершина на световом
#      конусе дедлайна (критический путь).
#   5) Циклы: ряд Неймана не обрывается => регулятор beta (затухание,
#      «масса»): K_{n+1} = phi + beta·A·K_n -> стоячая амплитуда.
#      Булев возврат K(v,v)>0 = детектор циклов (остаток Кана).
#   6) Свидетель arg (ветвь-победитель) = измерение: сворачивает сумму
#      по путям в одну мировую линию (восстановление траектории).
# ============================================================================
from collections import defaultdict

INF, NINF = float('inf'), float('-inf')

# ── Полукольца: «пространства», в которых живёт волна ──────────────────────
class Semiring:
    def __init__(self, zero, one, plus, times, name):
        self.zero, self.one = zero, one
        self.plus, self.times = plus, times
        self.name = name

BOOL  = Semiring(False, True, lambda a, b: a or b, lambda a, b: a and b,
                 "Б(∨,∧): достижимость / циклы")
MINP  = Semiring(INF, 0.0, min, lambda a, b: a + b,
                 "(min,+): кратчайший маршрут")
MAXP  = Semiring(NINF, 0.0, max, lambda a, b: a + b,
                 "(max,+): расписание / критический путь")
COUNT = Semiring(0.0, 1.0, lambda a, b: a + b, lambda a, b: a * b,
                 "(+,×): статсумма путей")
VITER = Semiring(0.0, 1.0, max, lambda a, b: a * b,
                 "(max,×): надёжность / Viterbi")
BOTT  = Semiring(NINF, INF, max, min,
                 "(max,min): узкое место")

# ── Носитель: граф задач ───────────────────────────────────────────────────
class TaskGraph:
    """V — задачи (d(v) — работа, p(v) — надёжность вершины);
       E — зависимости u→v: w — лаг, cost — цена, cap — пропускная способность."""
    def __init__(self, title=""):
        self.title = title
        self.adj, self.cost, self.cap = defaultdict(dict), defaultdict(dict), defaultdict(dict)
        self.dur, self.prob = {}, {}

    def add(self, v, dur=0.0, prob=1.0):
        self.dur[v], self.prob[v] = dur, prob
        return self

    def dep(self, u, v, w=0.0, cost=None, cap=None):
        self.adj[u][v] = w
        self.cost[u][v] = cost if cost is not None else w
        self.cap[u][v] = cap if cap is not None else INF
        return self

    @property
    def V(self):
        return list(self.dur)

    def preds(self, v):
        return [u for u in self.V if v in self.adj[u]]

    def flip(self):
        """Обратный граф G^T (носитель опережающего пропагатора K-)."""
        H = TaskGraph(self.title + "^T")
        for v in self.V:
            H.add(v, self.dur[v], self.prob[v])
        for u in self.adj:
            for v, w in self.adj[u].items():
                H.dep(v, u, w, self.cost[u][v], self.cap[u][v])
        return H

# ── Порядок волны: Кан = та же развёртка в Б(∨,∧) ──────────────────────────
def kahn(G):
    indeg = {v: len(G.preds(v)) for v in G.V}
    order, stack = [], [v for v in G.V if indeg[v] == 0]
    while stack:
        u = stack.pop()
        order.append(u)
        for v in G.adj[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                stack.append(v)
    return order, sorted(v for v in G.V if indeg[v] > 0)   # остаток = циклы

# ============================================================================
#  ★ ОДИН АЛГОРИТМ ★  Пропагаторная развёртка
#     K(v) = nu_v( phi(v) (+) SUM_{u->v} eps_uv(K(u)) )
#  phi : V→S  — граничное возбуждение (дельта-источники);
#  eps : (u,v,Ku)→S — функтор ребра (по умолчанию перенос ⊗ w);
#  nu  : (v,x)→S    — функтор вершины (по умолчанию тождественный).
#  На DAG одна проходка = точное решение K = (1-A)^{-1} phi.
# ============================================================================
def sweep(G, phi, sr, eps=None, nu=None):
    order, cyc = kahn(G)
    K, pred = {}, {}
    for v in order:
        acc, arg = phi[v], None
        for u in G.preds(v):
            val = eps(u, v, K[u]) if eps else sr.times(K[u], G.adj[u][v])
            acc = sr.plus(acc, val)          # интерференция путей
            if acc == val:
                arg = u                      # свидетель (ветвь-победитель)
        K[v] = nu(v, acc) if nu else acc
        pred[v] = arg
    return K, pred, order, cyc

def worldline(v, pred):
    """Измерение: сумма по путм -> одна мировая линия (pred-цепь)."""
    path = []
    while v is not None:
        path.append(v)
        v = pred[v]
    return path[::-1]

def damped(G, phi, sr, beta=0.9, iters=1000, tol=1e-15):
    """Циклический режим: регуляризованный ряд Неймана
       K_{n+1} = phi (+) beta·A·K_n  ->  стоячая амплитуда."""
    K = dict(phi)
    n = 0
    for n in range(1, iters + 1):
        K2 = dict(phi)
        for u in G.V:
            for v, w in G.adj[u].items():
                K2[v] = sr.plus(K2[v], sr.times(beta * K[u], w))
        delta = max(abs(K2[v] - K[v]) for v in G.V)
        K = K2
        if delta < tol:
            break
    return K, n

# ── Демонстрация: движок применяет себя к производству этого ответа ────────
def main():
    meta = {   # вершина: (работа d, надёжность p, название)
        "A": (1.0, .99, "Постановка задачи"),
        "B": (2.0, .95, "Вывод ядра K(v)"),
        "C": (1.0, .98, "Реестр полуколец"),
        "D": (3.0, .90, "Реализация движка"),
        "E": (1.0, .99, "Прогон тестов"),
        "F": (1.0, .98, "Демо-прогон"),
        "G": (1.0, .99, "Отчёт (этот ответ)"),
    }
    G = TaskGraph("Движок·самоприменение")
    for k, (d, p, _) in meta.items():
        G.add(k, dur=d, prob=p)
    G.dep("A", "B", w=0.0, cost=2, cap=9)
    G.dep("A", "C", w=0.0, cost=1, cap=7)
    G.dep("B", "D", w=0.0, cost=1, cap=8)
    G.dep("C", "D", w=0.5, cost=1, cap=5)
    G.dep("C", "G", w=0.0, cost=5, cap=4)
    G.dep("D", "E", w=0.0, cost=1, cap=6)
    G.dep("E", "F", w=0.0, cost=1, cap=6)
    G.dep("F", "G", w=0.0, cost=1, cap=7)

    bar = "=" * 72
    print(bar)
    print(" TASK GRAPH ENGINE · ОБЪЕКТ=ПРЕДМЕТ «ПРОПАГАТОР»")
    print(" ОДИН АЛГОРИТМ:  K(v) = nu_v( phi(v) (+) SUM_{u->v} eps_uv(K(u)) )")
    print(bar)

    # [1] Порядок волны + детектор циклов
    order, cyc = kahn(G)
    print("\n[1] ПОРЯДОК РАЗВЁРТКИ (Кан):", " → ".join(order))
    print("    возвратные амплитуды K(v,v):",
          "циклов нет" if not cyc else "ЦИКЛЫ: " + " ".join(cyc))

    # [2] Достижимость: S = Б(∨,∧), эпсилон нейтрален к весам
    K, *_ = sweep(G, {v: v == "A" for v in G.V}, BOOL,
                  eps=lambda u, v, Ku: Ku)
    print("\n[2] ДОСТИЖИМОСТЬ из A (Б(∨,∧)):",
          " ".join(sorted(v for v in G.V if K[v])))

    # [3][4][5] Расписание: запаздывающий K+, опережающий K-, резерв R
    phi = {v: NINF for v in G.V}
    nu_dur = lambda v, x: G.dur[v] + max(0.0, x)
    Kp, pred, _, _ = sweep(G, phi, MAXP, nu=nu_dur)
    Km, _, _, _ = sweep(G.flip(), phi, MAXP, nu=nu_dur)
    T = max(Kp.values())
    R = {v: T - Kp[v] - Km[v] + G.dur[v] for v in G.V}
    crit = worldline(max(Kp, key=Kp.get), pred)
    print("\n[3] ЗАПАЗДЫВАЮЩИЙ ПРОПАГАТОР K+ (max,+), ранние финишы:")
    print("    " + "  ".join(f"{v}={Kp[v]:g}" for v in sorted(Kp)), f" (T={T:g})")
    print("[4] ОПЕРЕЖАЮЩИЙ ПРОПАГАТОР K- (тот же алгоритм на G^T), высоты:")
    print("    " + "  ".join(f"{v}={Km[v]:g}" for v in sorted(Km)))
    print("[5] РЕЗЕРВ R = T − K+ − K- + d  (световой конус дедлайна):")
    print("    " + "  ".join(f"{v}={R[v]:g}" for v in sorted(R)))
    print("    КРИТИЧЕСКИЙ ПУТЬ (мировая линия):", " → ".join(crit))

    # [6] Кратчайший маршрут: (min,+), цена рёбер cost
    KM, predm, _, _ = sweep(G, {v: 0.0 if v == "A" else INF for v in G.V},
                            MINP, eps=lambda u, v, Ku: Ku + G.cost[u][v])
    print("\n[6] КРАТЧАЙШИЙ МАРШРУТ (min,+):", f"{KM['G']:g} по",
          " → ".join(worldline("G", predm)))

    # [7] Статсумма путей: (+,×)
    KC, *_ = sweep(G, {v: 1.0 if v == "A" else 0.0 for v in G.V}, COUNT,
                   eps=lambda u, v, Ku: Ku)
    print("[7] СТАТСУММА — ЧИСЛО ПУТЕЙ A→G (+,×):", f"{KC['G']:g}")

    # [8] Надёжность: (max,×) — Viterbi-траектория
    KV, predv, _, _ = sweep(G, {v: 1.0 if v == "A" else 0.0 for v in G.V},
                            VITER, eps=lambda u, v, Ku: Ku,
                            nu=lambda v, x: x * G.prob[v])
    print("[8] НАДЁЖНОСТЬ (max,×):", f"{KV['G']:.4f} по",
          " → ".join(worldline("G", predv)))

    # [9] Узкое место: (max,min) — пропускная способность
    KB, predb, _, _ = sweep(G, {v: INF if v == "A" else NINF for v in G.V},
                            BOTT, eps=lambda u, v, Ku: min(Ku, G.cap[u][v]))
    print("[9] УЗКОЕ МЕСТО (max,min):", f"пропуск {KB['G']:g} по",
          " → ".join(worldline("G", predb)))

    # [10] Циклический режим: регуляризованный ряд Неймана
    C3 = TaskGraph("3-цикл")
    for k in ("E0", "Z1", "Z2", "Z3"):
        C3.add(k, dur=1.0)
    C3.dep("E0", "Z1", 1.0)
    C3.dep("Z1", "Z2", 1.0)
    C3.dep("Z2", "Z3", 1.0)
    C3.dep("Z3", "Z1", 1.0)
    _, cyc3 = kahn(C3)
    Kd, n = damped(C3, {"E0": 1.0, "Z1": 0.0, "Z2": 0.0, "Z3": 0.0},
                   COUNT, beta=0.9)
    print("\n[10] ЦИКЛЫ: остаток Кана =", " ".join(cyc3),
          "→ ряд Неймана не обрывается")
    print("     регулятор β=0.9: K(Zi) = {:.6f} (теория β/(1-β^3) = {:.6f}),".format(
        Kd["Z1"], 0.9 / (1 - 0.9 ** 3)))
    print(f"     сошлось за {n} итераций — стоячая амплитуда (усиление петли β^L)")

    # [11] Самопроверка движка
    assert all(abs(R[v]) < 1e-9 for v in crit), "критический путь ↔ резерв"
    assert KC["G"] == 3.0, "число путей A→G"
    assert abs(KV["G"] - .99 * .98 * .99) < 1e-12, "Viterbi: A→C→G"
    assert KM["G"] == 5.0 and KB["G"] == 6.0
    assert abs(Kd["Z1"] - 0.9 / (1 - 0.9 ** 3)) < 1e-9, "стоячая волна β/(1-β^L)"
    print("\n[11] САМОПРОВЕРКА: все утверждения выполнены — движок верен.")

    print("\n" + bar)
    print(" РЕЕСТР: ОДНО ЯДРО → ВСЕ ФУНКЦИИ СТРУКТУРЫ")
    print("   S=(∨,∧)      достижимость · замыкание · K(v,v)=циклы")
    print("   S=(max,+)    E(v) расписание · критический путь (K+)")
    print("   S=(max,+)G^T H(v) поздние сроки · резерв R=T−K+−K-+d (K-)")
    print("   S=(min,+)    кратчайший маршрут")
    print("   S=(+,×)      статсумма: число путей")
    print("   S=(max,×)    надёжность · Viterbi-мировая линия")
    print("   S=(max,min)  пропускная способность · узкое место")
    print("   Нейман+β     циклический режим: стоячие амплитуды")
    print(bar)

if __name__ == "__main__":
    main()
