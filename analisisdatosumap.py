pip install -U loompy

pip install anndata

pip install anndata numpy matplotlib scipy requests

pip install scanpy

import anndata
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from scipy.spatial import ConvexHull
import requests
import io

# Función para descargar los datos
def descargar_datos(url):
    response = requests.get(url)
    if response.status_code == 200:
        return io.BytesIO(response.content)
    else:
        raise Exception(f"Error al descargar los datos. Código de estado: {response.status_code}")

# Función para calcular el hull convexo
def calcular_hull_convexo(puntos):
    hull = ConvexHull(puntos)
    return puntos[hull.vertices]

# URL para descargar el conjunto de datos de las coordenadas
data_url = "https://datasets.cellxgene.cziscience.com/1c2d14d8-32d4-41be-b38d-ba975ad10efa.h5ad"

print("Descargando y cargando los datos...")
archivo_datos = descargar_datos(data_url)
adata = anndata.read_h5ad(archivo_datos)

# Verificar si 'X_umap' existe en adata.obsm, si no, calcularlo
if 'X_umap' not in adata.obsm.keys():
    # Calcular el grafo de vecinos
    sc.pp.neighbors(adata)
    # Calcular las coordenadas UMAP y almacenarlas en adata.obsm
    sc.tl.umap(adata)  # Esto almacenará las coordenadas UMAP en adata.obsm['X_umap']

# Extraer coordenadas UMAP e ID de clúster
print("Extrayendo coordenadas UMAP e IDs de clúster...")
coordenadas_umap = adata.obsm['X_umap']
ids_cluster = adata.obs['cluster_id']

# Obtener colores únicos para cada clúster
clusters_unicos = np.unique(ids_cluster)
colores = plt.cm.rainbow(np.linspace(0, 1, len(clusters_unicos)))

# Crear la figura
plt.figure(figsize=(12, 8))

# Graficar puntos y calcular hull convexo para cada clúster
print("Graficando datos y calculando hulls convexos...")
for i, cluster in enumerate(clusters_unicos):
    puntos_cluster = coordenadas_umap[ids_cluster == cluster]

    # Graficar puntos
    plt.scatter(puntos_cluster[:, 0], puntos_cluster[:, 1], c=[colores[i]], label=f'Clúster {cluster}', alpha=0.6)

    # Calcular y graficar hull convexo
    if len(puntos_cluster) > 2:  # Necesitamos al menos 3 puntos para un hull convexo
        puntos_hull = calcular_hull_convexo(puntos_cluster)
        plt.fill(puntos_hull[:, 0], puntos_hull[:, 1], color=colores[i], alpha=0.2)

plt.title('Proyección UMAP de Tipos de Células Cerebrales con Hulls Convexos')
plt.xlabel('UMAP1')
plt.ylabel('UMAP2')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Guardar la figura
print("Guardando la figura...")
plt.savefig('proyeccion_umap_celulas_cerebrales_con_hulls.png', dpi=300, bbox_inches='tight')

print("Análisis completo. Figura guardada como 'proyeccion_umap_celulas_cerebrales_con_hulls.png'")