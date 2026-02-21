from __future__ import annotations

import time, json, math, urllib.parse, urllib.request
from typing import List, Tuple, Optional

import matplotlib.pyplot as plt
import geopandas as gpd
try:
    import contextily as ctx
    USE_MAP = True
except:
    USE_MAP = False
from shapely.geometry import Point as ShpPoint 

import geometria_toolkit as gt
import SEC

Point2D = Tuple[float, float]
Circle = Tuple[Point2D, float]

PARROQUIAS_URBANAS_32 = [
    "Belisario Quevedo, Quito, Ecuador",
    "Carcelén, Quito, Ecuador",
    "Centro Histórico, Quito, Ecuador",
    "Cochapamba, Quito, Ecuador",
    "Comité del Pueblo, Quito, Ecuador",
    "Cotocollao, Quito, Ecuador",
    "Chilibulo, Quito, Ecuador",
    "Chillogallo, Quito, Ecuador",
    "Chimbacalle, Quito, Ecuador",
    "El Condado, Quito, Ecuador",
    "Guamaní, Quito, Ecuador",
    "Iñaquito, Quito, Ecuador",
    "Itchimbía, Quito, Ecuador",
    "Jipijapa, Quito, Ecuador",
    "Kennedy, Quito, Ecuador",
    "La Argelia, Quito, Ecuador",
    "La Concepción, Quito, Ecuador",
    "La Ecuatoriana, Quito, Ecuador",
    "La Ferroviaria, Quito, Ecuador",
    "La Libertad, Quito, Ecuador",
    "La Magdalena, Quito, Ecuador",
    "La Mena, Quito, Ecuador",
    "Mariscal Sucre, Quito, Ecuador",
    "Ponceano, Quito, Ecuador",
    "Puengasí, Quito, Ecuador",
    "Quitumbe, Quito, Ecuador",
    "Rumipamba, Quito, Ecuador",
    "San Bartolo, Quito, Ecuador",
    "San Isidro del Inca, Quito, Ecuador",
    "San Juan, Quito, Ecuador",
    "Solanda, Quito, Ecuador",
    "Turubamba, Quito, Ecuador",
]

PARROQUIAS_RURALES_33 = [
    "Alangasí, Quito, Ecuador",
    "Amaguaña, Quito, Ecuador",
    "Atahualpa, Quito, Ecuador",
    "Calacalí, Quito, Ecuador",
    "Calderón, Quito, Ecuador",
    "Chavezpamba, Quito, Ecuador",
    "Checa, Quito, Ecuador",
    "Conocoto, Quito, Ecuador",
    "Cumbayá, Quito, Ecuador",
    "El Quinche, Quito, Ecuador",
    "Gualea, Quito, Ecuador",
    "Guangopolo, Quito, Ecuador",
    "Guayllabamba, Quito, Ecuador",
    "La Merced, Quito, Ecuador",
    "Llano Chico, Quito, Ecuador",
    "Lloa, Quito, Ecuador",
    "Mindo, Quito, Ecuador",
    "Nanegal, Quito, Ecuador",
    "Nanegalito, Quito, Ecuador",
    "Nayón, Quito, Ecuador",
    "Nono, Quito, Ecuador",
    "Pacto, Quito, Ecuador",
    "Perucho, Quito, Ecuador",
    "Pifo, Quito, Ecuador",
    "Pintag, Quito, Ecuador",
    "Pomasqui, Quito, Ecuador",
    "Puéllaro, Quito, Ecuador",
    "Puembo, Quito, Ecuador",
    "San Antonio de Pichincha, Quito, Ecuador",
    "San José de Minas, Quito, Ecuador",
    "Tababela, Quito, Ecuador",
    "Tumbaco, Quito, Ecuador",
    "Yaruquí, Quito, Ecuador",
]

PARROQUIAS = PARROQUIAS_URBANAS_32 + [p for p in PARROQUIAS_RURALES_33 if "Pedro Vicente" not in p]

# ========= 2) Geocoding (Nominatim) =========
NOMINATIM_URL = "https://nominatim.openstreetmap.org/search"
USER_AGENT = "SEC-Quito-Universidad/1.0 (fteran@estud.usfq.edu.ec)"
SLEEP_SEC = 1.1

def geocode_one(query: str) -> Optional[Point2D]:
    params = {"q": query, "format": "json", "limit": "1"}
    url = NOMINATIM_URL + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=30) as resp:
        data = json.loads(resp.read().decode("utf-8"))
    if not data:
        return None
    lat = float(data[0]["lat"])
    lon = float(data[0]["lon"])
    return (lat, lon)

def get_points_latlon(names: List[str]) -> Tuple[List[str], List[Point2D]]:
    ok_names, pts = [], []
    for nm in names:
        p = geocode_one(nm)
        if p is not None:
            ok_names.append(nm)
            pts.append(p)
        time.sleep(SLEEP_SEC)
    return ok_names, pts

# ========= 4) Área del convex hull (en km^2) =========
def polygon_area(poly: List[Point2D]) -> float:
    if len(poly) < 3:
        return 0.0
    area2 = 0.0
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        area2 += x1 * y2 - y1 * x2
    return abs(area2) / 2.0

def main():
    # 1) Geocode
    names_ok, pts_latlon = get_points_latlon(PARROQUIAS)
    if len(pts_latlon) < 50:
        raise RuntimeError(f"Solo obtuve {len(pts_latlon)} puntos. Agrega más nombres o revisa geocoding.")

    # 2) GeoDataFrame (lon,lat) -> proyectar a metros (EPSG:3857)
    gdf = gpd.GeoDataFrame(
        {"name": names_ok},
        geometry=[ShpPoint(lon, lat) for (lat, lon) in pts_latlon], 
        crs="EPSG:4326"
    ).to_crs("EPSG:3857")

    # 3) Puntos en metros y km (para tu SEC)
    points_m = [(geom.x, geom.y) for geom in gdf.geometry]
    points_km = [(x / 1000.0, y / 1000.0) for (x, y) in points_m]

    # 4) SEC + Hull (en km)
    C_km = SEC.welzl(points_km)
    hull_km = gt.convex_hull(points_km)

    (cx_km, cy_km), r_km = C_km
    cx_m, cy_m = cx_km * 1000.0, cy_km * 1000.0

    centro_m = ShpPoint(cx_m, cy_m)
    centro_gdf = gpd.GeoSeries([centro_m], crs="EPSG:3857").to_crs("EPSG:4326")
    lon_c, lat_c = centro_gdf.iloc[0].x, centro_gdf.iloc[0].y

    print(f"Centro óptimo (lat, lon): ({lat_c:.6f}, {lon_c:.6f})")
    print(f"Radio mínimo: {r_km:.3f} km")
    centro_gdf = gpd.GeoSeries([centro_m], crs="EPSG:3857").to_crs("EPSG:4326")
    lon_c, lat_c = centro_gdf.iloc[0].x, centro_gdf.iloc[0].y

    print(f"Centro óptimo (lat, lon): ({lat_c:.6f}, {lon_c:.6f})")

    # Áreas (km^2)
    (_, _), r_km = C_km
    area_hull = polygon_area(hull_km)
    area_circle = math.pi * (r_km ** 2)
    desperdicio_pct = 0.0 if area_hull == 0 else (area_circle - area_hull) / area_hull * 100.0

    print(f"Puntos usados: {len(points_km)}")
    print(f"Radio mínimo: {r_km:.3f} km")
    print(f"Área hull: {area_hull:.3f} km^2")
    print(f"Área círculo: {area_circle:.3f} km^2")
    print(f"Desperdicio vs hull: {desperdicio_pct:.2f}%")

    # 5) Convertir hull y círculo a metros para dibujar encima del mapa
    (cx_km, cy_km), r_km = C_km
    cx, cy, r = cx_km * 1000.0, cy_km * 1000.0, r_km * 1000.0
    hull_m = [(x * 1000.0, y * 1000.0) for (x, y) in hull_km]

    # 6) Plot sobre mapa real
    fig, ax = plt.subplots(figsize=(8, 8))

    xs = [x for x, y in points_m]
    ys = [y for x, y in points_m]
    pad = 20000  # 3 km
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)

    ax.set_aspect('equal', adjustable='box')

    if USE_MAP:
        ctx.add_basemap(ax, source=ctx.providers.OpenStreetMap.Mapnik)

    ax.scatter(xs, ys, s=20)

    if len(hull_m) >= 2:
        hx = [p[0] for p in hull_m] + [hull_m[0][0]]
        hy = [p[1] for p in hull_m] + [hull_m[0][1]]
        ax.plot(hx, hy, linewidth=2)

    t = [2 * math.pi * i / 360 for i in range(361)]
    cxp = [cx + r * math.cos(a) for a in t]
    cyp = [cy + r * math.sin(a) for a in t]
    ax.plot(cxp, cyp, linewidth=2)

    ax.set_title("Quito: puntos + convex hull + SEC (sobre mapa real)")
    plt.show()

if __name__ == "__main__":
    main()
