

############################################
# 0) Paquetes
############################################
# install.packages(c("readxl","FNN","igraph","cluster"))
library(readxl)
library(FNN)
library(igraph)
library(cluster)

#df<-read_excel("C:/Users/COJEDIAZ/Downloads/RESULTADOS/KL_matrix2.xlsx", ) # interior

df<-read_excel("C:/Users/COJEDIAZ/Downloads/RESULTADOS/KL_matrix2.xlsx",sheet = "Hoja2") # superior



# Asumo 1ra columna = nombres; filas y columnas en mismo orden
monedas   <- df[[1]]
KL_matrix <- as.matrix(df[,-1])
rownames(KL_matrix) <- monedas
colnames(KL_matrix) <- monedas

# Si tu KL no es perfectamente simétrica, puedes promediar (opcional):
# KL_matrix <- (KL_matrix + t(KL_matrix)) / 2

KL_dist <- as.dist(KL_matrix)

# Proyección MDS para visualizar / K-means / FNN
set.seed(123)
coords <- cmdscale(KL_dist, k = 2)
if (is.null(rownames(coords))) rownames(coords) <- monedas

############################################
# 2) ANÁLISIS PARA ELEGIR K
#    2.1 Codo (K-means sobre coords)
#    2.2 Silhouette (hclust sobre KL_dist)
############################################
N    <- nrow(coords)
kmax <- max(3, min(10, N-1))     # límite superior razonable

## --- 2.1 Codo con K-means (coords) ---
set.seed(123)
wss <- numeric(kmax)
for (k in 1:kmax) {
  km <- kmeans(coords, centers = k, nstart = 50)
  wss[k] <- km$tot.withinss
}

plot(1:kmax, wss, type="b", pch=19,
     xlab="Número de clusters K",
     ylab="Total within-cluster sum of squares (WSS)",
     main="Gráfico de codo - K-means sobre MDS")

## --- 2.2 Silhouette con hclust (KL_dist) ---
hc_all <- hclust(KL_dist, method = "average")
ks <- 2:kmax
avg_sil <- numeric(length(ks))
for (i in seq_along(ks)) {
  gr <- cutree(hc_all, k = ks[i])
  sil <- silhouette(gr, dist = KL_dist)
  avg_sil[i] <- mean(sil[, "sil_width"])
}

plot(ks, avg_sil, type="b", pch=19,
     xlab="Número de clusters K",
     ylab="Silhouette promedio",
     main="Evaluación con silhouette")

# Elegimos K por silhouette (puedes cambiar el criterio si lo prefieres)
K <- ks[which.max(avg_sil)]
message(sprintf("K elegido por silhouette = %d (silhouette promedio = %.3f)",
                K, max(avg_sil, na.rm=TRUE)))



############################################
# 3) GRÁFICOS FINALES USANDO EL K ELEGIDO
############################################

## --- 3.1 Dendrograma y corte en K ---
hc <- hc_all  # ya calculado arriba
plot(hc, main = sprintf("Dendrograma - KL divergence (K=%d)", K),
     xlab = "", sub = "")
rect.hclust(hc, k = K, border = 2:6)

grupos_hc <- cutree(hc, k = K)
print(grupos_hc)   # asignación de cluster por moneda

## --- 3.2 MDS coloreado por hclust (K elegido) ---
cols <- as.factor(grupos_hc[rownames(coords)])
plot(coords,
     col = cols, pch = 19,
     main = sprintf("MDS (cmdscale) coloreado por hclust K=%d", K),
     xlab = "Dim 1", ylab = "Dim 2")
text(coords, labels = rownames(coords), pos = 3, cex = 0.8)
legend("topright",
       legend = sort(unique(grupos_hc)),
       col = 1:length(unique(cols)), pch = 19, title = "Cluster")

## --- 3.3 K-means sobre coords (mismo K) ---
set.seed(123)
km <- kmeans(coords, centers = K, nstart = 50)
print(km$cluster)

plot(coords,
     col = km$cluster, pch = 19,
     main = sprintf("K-means en MDS (K=%d)", K),
     xlab = "Dim 1", ylab = "Dim 2")
text(coords, labels = rownames(coords), pos = 3, cex = 0.8)
points(km$centers, pch = 4, cex = 2, lwd = 2)
legend("topright", legend = 1:K, col = 1:K, pch = 19, title = "Cluster")

## --- 3.4 Grafo k-NN (k=3) coloreado por hclust (K elegido) ---
kNN <- 3
nn <- FNN::get.knn(coords, k = kNN)

edges <- do.call(rbind, lapply(seq_along(monedas), function(i){
  data.frame(
    from = monedas[i],
    to   = monedas[nn$nn.index[i, ]],
    d    = nn$nn.dist[i, ],
    stringsAsFactors = FALSE
  )
}))
edges$d <- as.numeric(edges$d)
edges$w <- 1/(1 + edges$d)   # similitud

g <- igraph::graph_from_data_frame(edges[, c("from","to","w")], directed = FALSE)
g <- igraph::simplify(g, edge.attr.comb = list(w = "mean"))

# colorear por clusters del dendrograma
V(g)$color <- as.numeric(as.factor(grupos_hc[V(g)$name]))
E(g)$width <- 2 * E(g)$w

lay <- coords[V(g)$name, , drop = FALSE]
plot(g,
     layout = as.matrix(lay),
     vertex.label = V(g)$name,
     vertex.size  = 10,
     edge.width   = E(g)$width,
     main = sprintf("Grafo k-NN (k=%d) coloreado por hclust K=%d", kNN, K))
