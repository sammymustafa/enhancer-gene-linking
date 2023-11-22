from setup import *

!wget -c https://www.dropbox.com/sh/d11ndep99xc6aka/AAByrCihV_V7Izes0wmdmiJLa?dl=1 -O linking.zip
!unzip -o linking.zip -d linking
os.unlink("linking.zip")

rna = anndata.read("linking/rna.h5ad")
h3k27ac = anndata.read("linking/H3K27ac_reduced.h5ad")
dnase = anndata.read("linking/DNase-seq_reduced.h5ad")