from typing import List, Tuple
import random
import math
import time
import geometria_toolkit as gt
import sys
sys.setrecursionlimit(30000)

Point = Tuple[float, float]
Circle = Tuple[Point, float]

def punto_en_circulo(p: Point, C: Circle, eps: float = gt.EPS) -> bool:
    (cx, cy), r = C
    dx = p[0] - cx
    dy = p[1] - cy
    return dx*dx + dy*dy <= (r + eps) * (r + eps)

def circulo_trivial(R: List[Point]) -> Circle:
    if len(R) == 0:
        return (0.0, 0.0), 0.0

    if len(R) == 1:
        return R[0], 0.0

    if len(R) == 2:
        (x1, y1), (x2, y2) = R
        centro = ((x1 + x2) / 2.0, (y1 + y2) / 2.0)
        radio = math.hypot(x2 - x1, y2 - y1) / 2.0
        return centro, radio

    else:
        p1, p2, p3 = R

        if gt.orient2d(p1, p2, p3) == 0:
            c12 = circulo_trivial([p1, p2])
            c13 = circulo_trivial([p1, p3])
            c23 = circulo_trivial([p2, p3])
            return max([c12, c13, c23], key=lambda c: c[1])

        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3

        d = 2.0 * (x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))

        if abs(d) <= gt.EPS:
            c12 = circulo_trivial([p1, p2])
            c13 = circulo_trivial([p1, p3])
            c23 = circulo_trivial([p2, p3])
            return max([c12, c13, c23], key=lambda c: c[1])

        s1 = x1*x1 + y1*y1
        s2 = x2*x2 + y2*y2
        s3 = x3*x3 + y3*y3

        ux = (s1*(y2 - y3) + s2*(y3 - y1) + s3*(y1 - y2)) / d
        uy = (s1*(x3 - x2) + s2*(x1 - x3) + s3*(x2 - x1)) / d

        centro = (ux, uy)
        radio = math.hypot(ux - x1, uy - y1)
        return centro, radio

def welzl(points: List[Point]) -> Circle:
    P = points[:]          
    random.shuffle(P)     

    def _welzl(n: int, R: List[Point]) -> Circle:
        if n == 0 or len(R) == 3:
            return circulo_trivial(R)

        p = P[n - 1]
        C = _welzl(n - 1, R)

        if punto_en_circulo(p, C):
            return C

        R.append(p)
        C2 = _welzl(n - 1, R)
        R.pop()
        return C2

    return _welzl(len(P), [])

# Enfoque determinista
def sec_determinista(points: List[Point]) -> Circle:
    if not points:
        return (0.0, 0.0), 0.0
    if len(points) == 1:
        return points[0], 0.0
    
    H = gt.convex_hull(points)
    h = len(H)

    mejor: Circle = ((0.0, 0.0), float("inf"))

    def contiene_todos(C: Circle) -> bool:
        return all(punto_en_circulo(p, C) for p in points)
    
    for i in range(h):
        for j in range(i + 1, h):
            C = circulo_trivial([H[i], H[j]])
            if C[1] < mejor[1] and contiene_todos(C):
                mejor = C

    for i in range(h):
        for j in range(i + 1, h):
            for k in range(j + 1, h):
                C = circulo_trivial([H[i], H[j], H[k]])
                if C[1] < mejor[1] and contiene_todos(C):
                    mejor = C

    return mejor

def gen_points_uniform(n, seed=0):
    rng = random.Random(seed)
    return [(rng.random(), rng.random()) for _ in range(n)]

def time_func(func, points, repeats=5):
    # mediana de repeats
    ts = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        func(points)
        ts.append(time.perf_counter() - t0)
    ts.sort()
    return ts[len(ts)//2]

def main():
    Ns = [10, 100, 1000, 10000]
    repeats = 7

    print(f"{'n':>7} | {'Welzl (s)':>10} | {'Det (s)':>10} | {'Det/Welzl':>10}")
    print("-"*48)

    for n in Ns:
        pts = gen_points_uniform(n, seed=12345 + n)

        tw = time_func(welzl, pts, repeats=repeats)

        rep_det = 3 if n >= 10000 else repeats
        td = time_func(sec_determinista, pts, repeats=rep_det)

        ratio = td / tw if tw > 0 else float("inf")
        print(f"{n:7d} | {tw:10.6f} | {td:10.6f} | {ratio:10.2f}x")

if __name__ == "__main__":
    main()
