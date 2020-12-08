# Update Seurat object form v2 to v3

> Hoang, 2020, science .Rdata was provided for Seurat v2.  
> Because of size (~7GB) impossible to update locally  
> Attempting to do so in Biowulf

Run in Terminal:

```bash
ssh bioqlulf.nih.gov
sinteractive --mem=32g --cpus-per-task=16
module load R
R
```

Once in R:
```R
setwd("/data/angueyraaristjm/2020_zfDevHoang/")
library(Seurat)
load("Zebrafish_development_pbmc_seurat.RData")
upPBMC = UpdateSeuratObject(pbmc)
pbmc = upPBMC
saveRDS(pbmc,"./zfDev_pbmc_v3.rds")
```