from __future__ import annotations

from fractions import Fraction
from typing import Tuple, Iterable, List, TypeVar, Optional, Set
from dataclasses import dataclass
import bisect
from heapq import heappush, heappop

EPS = 1e-9

Point = Tuple[float, float]
SegmentTuple = Tuple[Point, Point]

try:
    import gmpy2
except Exception:
    gmpy2 = None


def orient2d_float(p: Point, q: Point, r: Point) -> float:
    """
    Determinante 2D usando floats (rápido).

    El signo indica orientación:
      det > 0 -> izquierda (CCW)
      det < 0 -> derecha (CW)
      det = 0 -> colineal
    """
    (px, py), (qx, qy), (rx, ry) = p, q, r
    return (qx - px) * (ry - py) - (qy - py) * (rx - px)


def orient2d_fraction(p: Point, q: Point, r: Point) -> Fraction:
    """
    Determinante 2D exacto usando Fraction (fallback robusto).
    """
    (px, py), (qx, qy), (rx, ry) = p, q, r
    px, py = Fraction(px), Fraction(py)
    qx, qy = Fraction(qx), Fraction(qy)
    rx, ry = Fraction(rx), Fraction(ry)
    return (qx - px) * (ry - py) - (qy - py) * (rx - px)


def orient2d_gmpy2(p: Point, q: Point, r: Point):
    """
    Determinante 2D exacto usando gmpy2.mpq (fallback robusto más rápido que Fraction).

    Lanza RuntimeError si gmpy2 no está instalado.
    """
    if gmpy2 is None:
        raise RuntimeError("gmpy2 no está instalado.")
    (px, py), (qx, qy), (rx, ry) = p, q, r
    px, py = gmpy2.mpq(px), gmpy2.mpq(py)
    qx, qy = gmpy2.mpq(qx), gmpy2.mpq(qy)
    rx, ry = gmpy2.mpq(rx), gmpy2.mpq(ry)
    return (qx - px) * (ry - py) - (qy - py) * (rx - px)


def orient2d(p: Point, q: Point, r: Point, eps: float = 1e-12) -> int:
    """
    Predicado de orientación robusto (adaptativo con fallback).
    - Estima una escala típica del problema (magnitud de los productos en det).
    - Si |det| > eps * scale, el signo del float es confiable.
    - Si no, usa cálculo exacto (gmpy2 si existe, caso contrario Fraction).
    """
    det = orient2d_float(p, q, r)

    (x1, y1), (x2, y2), (x3, y3) = p, q, r
    scale = max(
        abs(x2 - x1) * abs(y3 - y1),
        abs(y2 - y1) * abs(x3 - x1),
        1.0
    )

    # Importante: usar ">" (si está cerca de 0, mejor NO confiar en float).
    if abs(det) > eps * scale:
        return 1 if det > 0 else -1

    # Fallback exacto
    if gmpy2 is not None:
        detx = orient2d_gmpy2(p, q, r)
        return 1 if detx > 0 else (-1 if detx < 0 else 0)

    detx = orient2d_fraction(p, q, r)
    return 1 if detx > 0 else (-1 if detx < 0 else 0)


def orientacion_str(signo: int) -> str:
    """
    Convierte el signo de orient2d() a un mensaje legible.
    """
    if signo > 0:
        return "Va hacia la izquierda"
    if signo < 0:
        return "Va hacia la derecha"
    return "Son colineales"


def cross(o: Point, a: Point, b: Point) -> float:
    """
    Producto cruzado 2D (orientación) de OA x OB.

    Retorna:
        > 0  si o->a->b es giro a la izquierda (CCW)
        < 0  si o->a->b es giro a la derecha (CW)
        = 0  si son colineales
    """
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_hull(points: Iterable[Point], *, include_collinear: bool = False) -> List[Point]:
    """
    Calcula la envolvente convexa (convex hull) de un conjunto de puntos en 2D
    usando Monotone Chain (O(n log n)).

    Returns:
        Lista de puntos del hull en orden antihorario (CCW), sin repetir el primer punto al final.
        Si hay 0 o 1 punto distinto, retorna esos puntos.

    """
    P = sorted(set(points))  # orden lexicográfico y elimina duplicados
    n = len(P)
    if n <= 1:
        return P

    # Para excluir colineales: quitamos cuando orient <= 0 (giro horario o colineal).
    # Para incluir colineales: quitamos cuando orient < 0 (solo giro horario).
    def should_pop(o: Point, a: Point, b: Point) -> bool:
        s = orient2d(o, a, b)   # -1,0,+1 robusto
        return s < 0 if include_collinear else s <= 0

    lower: List[Point] = []
    for p in P:
        while len(lower) >= 2 and should_pop(lower[-2], lower[-1], p):
            lower.pop()
        lower.append(p)

    upper: List[Point] = []
    for p in reversed(P):
        while len(upper) >= 2 and should_pop(upper[-2], upper[-1], p):
            upper.pop()
        upper.append(p)

    # lower y upper incluyen extremos duplicados (primer y último), los quitamos al concatenar
    hull = lower[:-1] + upper[:-1]

    # Caso degenerado: todos colineales y include_collinear=False => hull queda con 2 puntos (extremos).
    # Con include_collinear=True => devuelve todos los puntos del borde.
    return hull


def segment_intersection_point(a: Point, b: Point, c: Point, d: Point) -> Optional[Point]:
    """Intersección única de segmentos (ignora solapamientos colineales largos)."""
    def on_segment(p, q, r):
        return (min(p[0], q[0]) - EPS <= r[0] <= max(p[0], q[0]) + EPS and
                min(p[1], q[1]) - EPS <= r[1] <= max(p[1], q[1]) + EPS)

    def segments_intersect(p1, q1, p2, q2):
        o1 = orient2d(p1, q1, p2)
        o2 = orient2d(p1, q1, q2)
        o3 = orient2d(p2, q2, p1)
        o4 = orient2d(p2, q2, q1)

        # Caso general: signos opuestos estrictos
        if (o1 * o2 < 0) and (o3 * o4 < 0):
            return True
        if o1 == 0 and on_segment(p1, q1, p2): return True
        if o2 == 0 and on_segment(p1, q1, q2): return True
        if o3 == 0 and on_segment(p2, q2, p1): return True
        if o4 == 0 and on_segment(p2, q2, q1): return True
        return False

    if not segments_intersect(a, b, c, d):
        return None

    ax, ay = a; bx, by = b; cx, cy = c; dx, dy = d
    r = (bx-ax, by-ay)
    s = (dx-cx, dy-cy)
    rxs = r[0]*s[1] - r[1]*s[0]
    if abs(rxs) <= EPS:
        return None

    q_p = (cx-ax, cy-ay)
    t = (q_p[0]*s[1] - q_p[1]*s[0]) / rxs
    x = ax + t*r[0]
    y = ay + t*r[1]
    return (x, y)


@dataclass(frozen=True)
class Seg:
    i: int
    a: Point
    b: Point

    def upper_lower(self) -> Tuple[Point, Point]:
        # "Upper": mayor y; si empate, menor x
        if (self.a[1] > self.b[1]) or (abs(self.a[1]-self.b[1]) <= EPS and self.a[0] < self.b[0]):
            return self.a, self.b
        return self.b, self.a


def x_at_y(seg: Seg, y: float) -> float:
    """x donde el segmento corta la horizontal y (si y está en su rango)."""
    (x1, y1), (x2, y2) = seg.a, seg.b
    if abs(y2 - y1) <= EPS:
        return min(x1, x2)
    t = (y - y1) / (y2 - y1)
    return x1 + t*(x2 - x1)


def event_order(p: Point) -> Tuple[float, float]:
    return (-p[1], p[0])


def key_point(p: Point, ndigits: int = 9) -> Tuple[float, float]:
    return (round(p[0], ndigits), round(p[1], ndigits))


def segment_intersection(segments: Iterable[SegmentTuple], *, trace: bool = False) -> List[Tuple[Point, int, int]]:
    """
    Plane sweep (Bentley–Ottmann didáctico mejorado).
    Devuelve lista de (punto, i, j) con intersecciones entre segmentos i y j.
    """
    segs = [Seg(i, a, b) for i, (a, b) in enumerate(segments)]
    n = len(segs)

    # Event queue: (order_key, point, kind, payload)
    # kind: "UP", "LOW", "X"
    heap = []
    in_heap: Set[Tuple[Tuple[float, float], str, object]] = set()

    def push_event(p: Point, kind: str, payload):
        pr = key_point(p)
        k = (pr, kind, payload)
        if k in in_heap:
            return
        in_heap.add(k)
        heappush(heap, (event_order(pr), pr, kind, payload))

    # Insert endpoints
    for s in segs:
        up, low = s.upper_lower()
        push_event(up, "UP", s.i)
        push_event(low, "LOW", s.i)

    status: List[int] = []
    y_sweep: float = float("inf")
    delta = 10 * EPS

    # Para evitar duplicar intersecciones
    reported: Set[Tuple[Tuple[float, float], int, int]] = set()
    out: List[Tuple[Point, int, int]] = []

    def status_key(seg_id: int, y: float) -> float:
        return x_at_y(segs[seg_id], y)

    def locate(seg_id: int, y: float) -> int:
        """Posición donde insertar seg_id en status según x_at_y en y."""
        x = status_key(seg_id, y)
        keys = [status_key(sid, y) for sid in status]
        return bisect.bisect_left(keys, x)

    def neighbors(pos: int) -> Tuple[Optional[int], Optional[int]]:
        left = status[pos-1] if pos-1 >= 0 else None
        right = status[pos+1] if pos+1 < len(status) else None
        return left, right

    def schedule_intersection(i: Optional[int], j: Optional[int], p_event: Point):
        if i is None or j is None or i == j:
            return
        a, b = segs[i].a, segs[i].b
        c, d = segs[j].a, segs[j].b
        p = segment_intersection_point(a, b, c, d)
        if p is None:
            return

        px, py = p
        ex, ey = p_event

        # Solo programar eventos “por debajo” del evento actual,
        # o misma y pero a la derecha (convención estándar)
        if (py < ey - EPS) or (abs(py - ey) <= EPS and px > ex + EPS):
            ii, jj = (i, j) if i < j else (j, i)
            push_event(p, "X", (ii, jj))

    while heap:
        _, p, kind, payload = heappop(heap)

        # Agrupar todos los eventos con el mismo punto p
        batch = [(kind, payload)]
        while heap and heap[0][1] == p:
            _, _, k2, pl2 = heappop(heap)
            batch.append((k2, pl2))

        x_p, y_p = p
        y_sweep = y_p

        if trace:
            print(f"\nPUNTO EVENTO {p} -> batch={batch}")

        # Separar eventos por tipo
        ups = [pl for (k, pl) in batch if k == "UP"]
        lows = [pl for (k, pl) in batch if k == "LOW"]
        xs = [pl for (k, pl) in batch if k == "X"]

        # 1) Reportar intersecciones del batch
        for (i, j) in xs:
            kk = (key_point(p), i, j)
            if kk not in reported:
                reported.add(kk)
                out.append((p, i, j))
                if trace:
                    print(f"  -> report X entre s{i} y s{j} en {p}")

        # 2) Remover LOW (segmentos que terminan aquí)
        #    (remover antes de insertar ayuda a estabilidad)
        for sid in lows:
            if sid in status:
                pos = status.index(sid)
                left, right = (status[pos-1] if pos-1 >= 0 else None,
                               status[pos+1] if pos+1 < len(status) else None)
                status.pop(pos)
                # al remover, left y right se vuelven vecinos
                schedule_intersection(left, right, p)

        # 3) Insertar UP (segmentos que empiezan aquí)
        for sid in ups:
            y_probe = y_sweep - delta
            pos = locate(sid, y_probe)
            status.insert(pos, sid)
            left, right = neighbors(pos)
            schedule_intersection(left, sid, p)
            schedule_intersection(sid, right, p)

        # 4) Procesar swaps por intersección (X)
        #    En Bentley–Ottmann real, el orden cambia en el punto de intersección.
        #    Aquí hacemos una corrección local: reordenar los involucrados si ambos están en status.
        for (i, j) in xs:
            if i in status and j in status:
                pi = status.index(i)
                pj = status.index(j)
                if pi > pj:
                    pi, pj = pj, pi
                    i, j = j, i
                status[pi], status[pj] = status[pj], status[pi]

                # Revisar nuevos vecinos cerca de i y j
                for sid in (i, j):
                    pos = status.index(sid)
                    left = status[pos-1] if pos-1 >= 0 else None
                    right = status[pos+1] if pos+1 < len(status) else None
                    schedule_intersection(left, sid, p)
                    schedule_intersection(sid, right, p)

    return out
