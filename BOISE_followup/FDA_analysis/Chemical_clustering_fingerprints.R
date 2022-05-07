dist = read.csv('Similarity_mat_FDA.csv', header = T)
dist = as.matrix(dist)
dist = 1 - dist
cl = diana(x = dist, diss = T, stop.at.k = 20, keep.diss = T)
# tmp = NbClust(diss = cl$diss, distance = NULL, min.nc = 2, max.nc = 25, method = 'average', index = 'cindex')
hcd = as.dendrogram(cl)
tmp = cut(hcd, h = 0.95)
plot(tmp$upper)

hac = hclust(cl$diss, method = 'average')
hcd = as.dendrogram(hac)
tmp = cut(hcd, h = 0.95)
plot(tmp$upper)

## Silhouette score
silhouette_score <- function(k){
  tmp_cl = cutree(cl, k)
  ss <- silhouette(x = tmp_cl, dmatrix = dist)
  mean(ss[, 3])
}
k <- 2:25
avg_sil <- sapply(k, silhouette_score)
plot(k, type='b', avg_sil, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE)
tmp_cl = cutree(cl, k=19)
table(tmp_cl)
res = data.frame(cid = colnames(dist), clust = tmp_cl)
write.csv(res, file = 'chemical_clustering_res.csv', row.names = F)

tSNE = Rtsne(cl$diss, is_distance = T, dims = 2, perplecity = 50, max_iter=10000)

fps = read.csv('fp_FDA.csv')
fps = t(as.matrix(fps))
tSNE = Rtsne(fps, check_duplicates = F,dims = 2, perplecity = 25, max_iter=10000)
tsne_plot = data.frame(x= tSNE$Y[,1], y = tSNE$Y[,2], col = as.factor(tmp_cl))
p = ggplot(tsne_plot)+
  geom_point(aes(x=x, y=y, color = col))
p     
        