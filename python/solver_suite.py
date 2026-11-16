import numpy as np
from numba import njit


BUDGET = 1_960_400


@njit(cache=True)
def fitness(path, dose):
    value = dose[0, path[0]]
    for k in range(path.size - 1):
        value += dose[path[k], path[k + 1]]
    return value + dose[path[-1], 0]


@njit(cache=True)
def initialize_population(m, n):
    pop = np.empty((m, n), np.int64)
    fit = np.empty(m, np.float64)
    for i in range(m):
        pop[i] = np.random.permutation(n) + 1
    return pop, fit


@njit(cache=True)
def evaluate_population(pop, dose, out):
    for i in range(pop.shape[0]):
        out[i] = fitness(pop[i], dose)


@njit(cache=True)
def ox(a, b):
    n = a.size
    i = np.random.randint(0, n - 1)
    j = np.random.randint(i + 1, n)
    child = np.zeros(n, np.int64)
    used = np.zeros(n + 1, np.bool_)
    child[i:j + 1] = a[i:j + 1]
    for k in range(i, j + 1):
        used[a[k]] = True
    pos = 0
    for k in range(n):
        value = b[k]
        if not used[value]:
            if pos == i:
                pos = j + 1
            child[pos] = value
            pos += 1
    return child


@njit(cache=True)
def swap(path):
    i = np.random.randint(path.size)
    j = np.random.randint(path.size - 1)
    if j >= i:
        j += 1
    path[i], path[j] = path[j], path[i]


@njit(cache=True)
def insert(path):
    n = path.size
    i = np.random.randint(n)
    j = np.random.randint(n)
    while i == j:
        j = np.random.randint(n)
    value = path[i]
    if i < j:
        for k in range(i, j):
            path[k] = path[k + 1]
    else:
        for k in range(i, j, -1):
            path[k] = path[k - 1]
    path[j] = value


@njit(cache=True)
def invert(path):
    i = np.random.randint(0, path.size - 1)
    j = np.random.randint(i + 1, path.size)
    path[i:j + 1] = path[i:j + 1][::-1]


@njit(cache=True)
def scramble(path):
    i = np.random.randint(0, path.size - 1)
    j = np.random.randint(i + 1, path.size)
    segment = path[i:j + 1].copy()
    path[i:j + 1] = segment[np.random.permutation(segment.size)]


@njit(cache=True)
def mutate(path, ps, pi, pv, pr):
    if np.random.random() < ps:
        swap(path)
    if np.random.random() < pi:
        insert(path)
    if np.random.random() < pv:
        invert(path)
    if np.random.random() < pr:
        scramble(path)


@njit(cache=True)
def opt2(path, dose, attempts):
    accepted = 0
    n = path.size
    for _ in range(attempts):
        i = np.random.randint(0, n - 2)
        j = np.random.randint(i + 1, n - 1)
        a = 0 if i == 0 else path[i - 1]
        b = path[i]
        c = path[j]
        e = path[j + 1]
        if dose[a, c] + dose[b, e] < dose[a, b] + dose[c, e]:
            path[i:j + 1] = path[i:j + 1][::-1]
            accepted += 1
    return accepted


@njit(cache=True)
def curve_checkpoints(curve):
    out = np.empty(101, np.float64)
    last = curve.size - 1
    for i in range(101):
        out[i] = curve[int(round(i * last / 100.0))]
    return out


@njit(cache=True)
def run_cpo_dbo(dose, seed, m=400, tscale=30000.0, elite_ratio=0.02,
                twoopt_attempts=1, dynamic_probs=1, use_group_bias=1,
                use_elite=1, mutation_mode=0, use_pbest=1,
                budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    elite = max(2, int(round(elite_ratio * m))) if use_elite else 0
    iterations = int((budget - m) // (m - elite))
    pop, fits = initialize_population(m, n)
    evaluate_population(pop, dose, fits)
    initial = pop.copy()
    pbest = pop.copy()
    pfit = fits.copy()
    idx = np.argmin(fits)
    best = fits[idx]
    best_path = pop[idx].copy()
    curve = np.empty(iterations, np.float64)
    attempts = 0
    accepted = 0
    for it in range(iterations):
        t = min(1.0, (it + 1) / tscale) if dynamic_probs else 0.0
        pd = 0.70 * (1.0 - 0.60 * t)
        pa = 0.25 * (1.0 - 0.20 * t)
        pe = min(0.60, 0.05 + 0.55 * t)
        total = pd + pa + pe
        pd /= total
        pa /= total
        pe /= total
        order = np.argsort(fits, kind="mergesort")
        newpop = pop.copy()
        newfit = fits.copy()
        for e in range(elite):
            newpop[e] = pop[order[e]]
            newfit[e] = fits[order[e]]
        for i in range(elite, m):
            bias = 0.0
            if use_group_bias:
                if i < int(round(m * 0.40)):
                    bias = 0.10
                elif i < int(round(m * 0.60)):
                    bias = -0.05
                else:
                    bias = -0.10
            p_d = min(0.95, max(0.05, pd + bias))
            p_a = min(0.95, max(0.05, pa - bias / 2.0))
            scale = p_d + p_a + pe
            p_d /= scale
            p_a /= scale
            r = np.random.random()
            if r < p_d:
                child = ox(pop[i], best_path) if np.random.random() < 0.80 else pop[i].copy()
                if np.random.random() < 0.50:
                    swap(child)
                mode = 0
            elif r < p_d + p_a:
                mate = np.random.randint(m)
                donor = pbest[mate] if use_pbest else pop[mate]
                child = ox(pop[i], donor)
                mode = 1
            else:
                child = pop[i].copy()
                if np.random.random() < 0.80:
                    scramble(child)
                if np.random.random() < 0.60:
                    invert(child)
                mode = 2
            if mutation_mode == 1:
                if np.random.random() < (0.60 if mode == 0 else (0.35 if mode == 1 else 0.25)):
                    swap(child)
            else:
                if mode == 0:
                    mutate(child, 0.60, 0.15, 0.15, 0.0)
                elif mode == 1:
                    mutate(child, 0.35, 0.20, 0.30, 0.0)
                else:
                    mutate(child, 0.25, 0.25, 0.45, 0.0)
            if mode < 2 and twoopt_attempts > 0:
                attempts += twoopt_attempts
                accepted += opt2(child, dose, twoopt_attempts)
            value = fitness(child, dose)
            if value < fits[i]:
                newpop[i] = child
                newfit[i] = value
                if value < pfit[i]:
                    pbest[i] = child
                    pfit[i] = value
        pop = newpop
        fits = newfit
        idx = np.argmin(fits)
        if fits[idx] < best:
            best = fits[idx]
            best_path = pop[idx].copy()
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    evaluations = m + iterations * (m - elite)
    return best, best_path, curve_checkpoints(curve), first_hit, evaluations, attempts, accepted, initial


@njit(cache=True)
def run_cpo(dose, seed, m=400, budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    iterations = int((budget - m) // m)
    pop, fits = initialize_population(m, n)
    evaluate_population(pop, dose, fits)
    initial = pop.copy()
    pbest = pop.copy()
    pfit = fits.copy()
    idx = np.argmin(fits)
    best = fits[idx]
    best_path = pop[idx].copy()
    curve = np.empty(iterations)
    for it in range(iterations):
        t = (it + 1) / iterations
        pd = 0.70 * (1.0 - t)
        pa = 0.30 * (1.0 - t)
        pe = 0.20 + 0.30 * t
        total = pd + pa + pe
        pd /= total
        pa /= total
        for i in range(m):
            r = np.random.random()
            if r < pd:
                child = ox(pop[i], best_path) if np.random.random() < 0.70 else pop[i].copy()
                if np.random.random() < 0.50:
                    swap(child)
                mutate(child, 0.60, 0.15, 0.15, 0.10)
                opt2(child, dose, 1)
            elif r < pd + pa:
                child = ox(pop[i], pbest[np.random.randint(m)])
                mutate(child, 0.35, 0.20, 0.25, 0.10)
                opt2(child, dose, 1)
            else:
                child = pop[i].copy()
                if np.random.random() < 0.85:
                    scramble(child)
                if np.random.random() < 0.65:
                    invert(child)
                mutate(child, 0.25, 0.25, 0.35, 0.25)
            value = fitness(child, dose)
            if value < fits[i]:
                pop[i] = child
                fits[i] = value
                if value < pfit[i]:
                    pbest[i] = child
                    pfit[i] = value
        idx = np.argmin(pfit)
        if pfit[idx] < best:
            best = pfit[idx]
            best_path = pbest[idx].copy()
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    return best, best_path, curve_checkpoints(curve), first_hit, m + iterations * m, initial


@njit(cache=True)
def run_dbo(dose, seed, m=400, budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    iterations = int((budget - m) // m)
    pop, fits = initialize_population(m, n)
    evaluate_population(pop, dose, fits)
    initial = pop.copy()
    pbest = pop.copy()
    pfit = fits.copy()
    idx = np.argmin(fits)
    best = fits[idx]
    best_path = pop[idx].copy()
    curve = np.empty(iterations)
    ball = int(round(0.50 * m))
    breed = int(round(0.25 * m))
    for it in range(iterations):
        for i in range(m):
            child = pbest[i].copy()
            if i < ball:
                if np.random.random() < 0.60:
                    swap(child)
                if np.random.random() < 0.30:
                    child = ox(child, best_path)
                mutate(child, 0.30, 0.10, 0.10, 0.05)
            elif i < ball + breed:
                if np.random.random() < 0.60:
                    child = ox(child, best_path)
                mutate(child, 0.40, 0.15, 0.20, 0.10)
            else:
                if np.random.random() < 0.70:
                    insert(child)
                mutate(child, 0.35, 0.20, 0.25, 0.20)
            value = fitness(child, dose)
            if value < pfit[i]:
                pbest[i] = child
                pfit[i] = value
        idx = np.argmin(pfit)
        if pfit[idx] < best:
            best = pfit[idx]
            best_path = pbest[idx].copy()
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    return best, best_path, curve_checkpoints(curve), first_hit, m + iterations * m, initial


@njit(cache=True)
def tournament(fits, size=3):
    best = np.random.randint(fits.size)
    for _ in range(size - 1):
        c = np.random.randint(fits.size)
        if fits[c] < fits[best]:
            best = c
    return best


@njit(cache=True)
def run_ga(dose, seed, m=400, budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    iterations = int((budget - m) // m)
    pop, fits = initialize_population(m, n)
    evaluate_population(pop, dose, fits)
    initial = pop.copy()
    best = fits[np.argmin(fits)]
    best_path = pop[np.argmin(fits)].copy()
    curve = np.empty(iterations)
    elite = int(round(0.05 * m))
    for it in range(iterations):
        selected = np.empty_like(pop)
        for i in range(m):
            selected[i] = pop[tournament(fits)]
        children = selected.copy()
        for i in range(0, m - 1, 2):
            if np.random.random() < 0.85:
                children[i] = ox(selected[i], selected[i + 1])
                children[i + 1] = ox(selected[i + 1], selected[i])
        for i in range(m):
            if np.random.random() < 0.15:
                swap(children[i])
                if np.random.random() < 0.10:
                    invert(children[i])
        childfit = np.empty(m)
        evaluate_population(children, dose, childfit)
        order_old = np.argsort(fits, kind="mergesort")
        merged_pop = np.empty((m + elite, n), np.int64)
        merged_fit = np.empty(m + elite)
        merged_pop[:m] = children
        merged_fit[:m] = childfit
        for e in range(elite):
            merged_pop[m + e] = pop[order_old[e]]
            merged_fit[m + e] = fits[order_old[e]]
        keep = np.argsort(merged_fit, kind="mergesort")[:m]
        pop = merged_pop[keep].copy()
        fits = merged_fit[keep].copy()
        if fits[0] < best:
            best = fits[0]
            best_path = pop[0].copy()
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    return best, best_path, curve_checkpoints(curve), first_hit, m + iterations * m, initial


@njit(cache=True)
def inject_segment(path, guide):
    n = path.size
    i = np.random.randint(0, n - 1)
    j = np.random.randint(i + 1, n)
    return ox(guide, path) if (j - i) > n // 3 else ox(path, guide)


@njit(cache=True)
def run_ipso(dose, seed, m=400, budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    iterations = int((budget - m) // m)
    chaos_iters = min(200, iterations)
    chaos = np.random.random((m, n))
    pop = np.empty((m, n), np.int64)
    fits = np.empty(m)
    for i in range(m):
        pop[i] = np.argsort(chaos[i]) + 1
        fits[i] = fitness(pop[i], dose)
    initial = pop.copy()
    pbest = pop.copy()
    pfit = fits.copy()
    idx = np.argmin(fits)
    best = fits[idx]
    best_path = pop[idx].copy()
    curve = np.empty(iterations)
    for it in range(chaos_iters):
        chaos = 4.0 * chaos * (1.0 - chaos)
        for i in range(m):
            child = np.argsort(chaos[i]) + 1
            value = fitness(child, dose)
            if value < pfit[i]:
                pbest[i] = child
                pfit[i] = value
        idx = np.argmin(pfit)
        if pfit[idx] < best:
            best = pfit[idx]
            best_path = pbest[idx].copy()
        curve[it] = best
    for it in range(chaos_iters, iterations):
        for i in range(m):
            child = ox(pop[i], pbest[i])
            child = ox(child, best_path)
            if np.random.random() < 0.10:
                child = inject_segment(child, best_path)
            mutate(child, 0.60, 0.15, 0.15, 0.10)
            value = fitness(child, dose)
            pop[i] = child
            fits[i] = value
            if value < pfit[i]:
                pbest[i] = child
                pfit[i] = value
        idx = np.argmin(pfit)
        if pfit[idx] < best:
            best = pfit[idx]
            best_path = pbest[idx].copy()
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    return best, best_path, curve_checkpoints(curve), first_hit, m + iterations * m, initial


@njit(cache=True)
def run_ga_pso(dose, seed, m=400, budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    iterations = int((budget - m) // m)
    pop, fits = initialize_population(m, n)
    evaluate_population(pop, dose, fits)
    initial = pop.copy()
    pbest = pop.copy()
    pfit = fits.copy()
    idx = np.argmin(fits)
    best = fits[idx]
    best_path = pop[idx].copy()
    curve = np.empty(iterations)
    for it in range(iterations):
        children = np.empty_like(pop)
        for i in range(m):
            if i % 2 == 0:
                a = tournament(fits)
                b = tournament(fits)
                child = ox(pop[a], pop[b])
            else:
                child = ox(pop[i], pbest[i])
                child = ox(child, best_path)
            if np.random.random() < 0.20:
                swap(child)
            if np.random.random() < 0.10:
                invert(child)
            children[i] = child
        childfit = np.empty(m)
        evaluate_population(children, dose, childfit)
        for i in range(m):
            if childfit[i] < fits[i]:
                pop[i] = children[i]
                fits[i] = childfit[i]
            if childfit[i] < pfit[i]:
                pbest[i] = children[i]
                pfit[i] = childfit[i]
        idx = np.argmin(pfit)
        if pfit[idx] < best:
            best = pfit[idx]
            best_path = pbest[idx].copy()
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    return best, best_path, curve_checkpoints(curve), first_hit, m + iterations * m, initial


@njit(cache=True)
def aco_construct(pheromone, dose):
    n = dose.shape[0] - 1
    path = np.empty(n, np.int64)
    used = np.zeros(n + 1, np.bool_)
    current = 0
    for pos in range(n):
        candidates = np.empty(5, np.int64)
        count = 0
        tries = 0
        while count < 5 and tries < 100:
            node = np.random.randint(1, n + 1)
            tries += 1
            if not used[node]:
                duplicate = False
                for q in range(count):
                    if candidates[q] == node:
                        duplicate = True
                if not duplicate:
                    candidates[count] = node
                    count += 1
        if count == 0:
            for node in range(1, n + 1):
                if not used[node]:
                    candidates[0] = node
                    count = 1
                    break
        weights = np.empty(count)
        total = 0.0
        for q in range(count):
            node = candidates[q]
            w = pheromone[current, node] * (1.0 / (dose[current, node] + 1e-9))
            weights[q] = w
            total += w
        r = np.random.random() * total
        chosen = candidates[count - 1]
        acc = 0.0
        for q in range(count):
            acc += weights[q]
            if r <= acc:
                chosen = candidates[q]
                break
        path[pos] = chosen
        used[chosen] = True
        current = chosen
    return path


@njit(cache=True)
def run_aco_ga(dose, seed, m=400, budget=BUDGET):
    np.random.seed(seed)
    n = dose.shape[0] - 1
    iterations = int((budget - m) // m)
    pop, fits = initialize_population(m, n)
    evaluate_population(pop, dose, fits)
    initial = pop.copy()
    idx = np.argmin(fits)
    best = fits[idx]
    best_path = pop[idx].copy()
    curve = np.empty(iterations)
    pheromone = np.ones((n + 1, n + 1), np.float64)
    for it in range(iterations):
        children = np.empty_like(pop)
        ant_count = int(round(0.10 * m)) if it % 20 == 0 else 0
        for i in range(ant_count):
            children[i] = aco_construct(pheromone, dose)
        for i in range(ant_count, m):
            a = tournament(fits)
            b = tournament(fits)
            child = ox(pop[a], pop[b])
            if np.random.random() < 0.15:
                swap(child)
            if np.random.random() < 0.08:
                invert(child)
            children[i] = child
        childfit = np.empty(m)
        evaluate_population(children, dose, childfit)
        elite_count = max(1, int(round(0.05 * m)))
        merged_pop = np.empty((m + elite_count, n), np.int64)
        merged_fit = np.empty(m + elite_count)
        merged_pop[:m] = children
        merged_fit[:m] = childfit
        old_order = np.argsort(fits, kind="mergesort")
        for e in range(elite_count):
            merged_pop[m + e] = pop[old_order[e]]
            merged_fit[m + e] = fits[old_order[e]]
        keep = np.argsort(merged_fit, kind="mergesort")[:m]
        pop = merged_pop[keep].copy()
        fits = merged_fit[keep].copy()
        if fits[0] < best:
            best = fits[0]
            best_path = pop[0].copy()
        if it % 20 == 0:
            pheromone *= 0.90
            top = min(5, m)
            for e in range(top):
                p = pop[e]
                amount = 1.0 / (fits[e] + 1e-9)
                prev = 0
                for k in range(n):
                    node = p[k]
                    pheromone[prev, node] += amount
                    pheromone[node, prev] += amount
                    prev = node
                pheromone[prev, 0] += amount
                pheromone[0, prev] += amount
        curve[it] = best
    first_hit = int(np.where(curve <= best + 1e-12)[0][0] + 1)
    return best, best_path, curve_checkpoints(curve), first_hit, m + iterations * m, initial
