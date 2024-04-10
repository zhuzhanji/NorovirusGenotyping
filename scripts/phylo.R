
if(system.file(package='phangorn') == ""){
  install.packages("phangorn")
}
library(phangorn)

if(system.file(package='seqinr') == ""){
  install.packages("seqinr")
}
library(seqinr)


## set command line arguments ----
args <- commandArgs(trailingOnly = TRUE)

#stop the script if no command line argument
if(length(args)==0){
  stop("Requires command line argument.")
}


noro <- read.alignment(file =  args[1], format = "phylip", type = "DNA")
noro2 <- as.phyDat(noro, type="DNA")


fitGTR <- pml_bb(noro2, model="HKY")

# bs value should be set to 100 according to [Kroneman 11]
bs <- bootstrap.pml(fitGTR, bs=5, optNni=TRUE,
                    control = pml.control(trace = 0))

# assigning standard bootstrap values to our tree; this is the default method
# supporting values are written to tree_stdbs as labels 
tree_stdbs <- plotBS(fitGTR$tree, bs, type = "n")

# save tree to file
write.tree(tree_stdbs, file=args[2])


