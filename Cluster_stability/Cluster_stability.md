Cluster_stability
================
M, Marcao
2024-06-17

``` r
# Load
adata_obs <- fread("/mnt/scratch1/maycon/genesis_embryo_lego/round_4/objects/cluster0_subclustering_multiple_res.csv", header = TRUE) %>% as.data.frame()
head(adata_obs)[, 14:17]
```

    ##   cl0_leiden_02 cl0_leiden_05 cl0_leiden_1.0 cl0_leiden_1.5
    ## 1             0             1              0              0
    ## 2             0             4              7             11
    ## 3             1             3              6              5
    ## 4             1             5              4              4
    ## 5             0             2              2              1
    ## 6             1             3             10             10

``` r
# Subset
cluster_df_sub <- adata_obs[, c("cl0_leiden_02",
                                "cl0_leiden_05",
                                "cl0_leiden_1.0",
                                "cl0_leiden_1.5")]
# Change column names
names(cluster_df_sub)
```

    ## [1] "cl0_leiden_02"  "cl0_leiden_05"  "cl0_leiden_1.0" "cl0_leiden_1.5"

``` r
names(cluster_df_sub) <- c("leiden_row_1", "leiden_row_2", "leiden_row_3", "leiden_row_4") # same order as in names(cluster_df_sub)
```

``` r
# Run cluster stability 
library(clustree)
clustree_out <- clustree(cluster_df_sub[, ], prefix = "leiden_row_", node_colour = "sc3_stability")

clustree_out
```

![Alt Text](ToolBox/Cluster_stability/Run clustree-1.png)<!-- -->

``` r
# Save cluster purity output 
# saveRDS(clustree_out, , file = "~your/path/clustree_output.rds")
```
