
data = read.csv("data/dseq.csv")

cat(apply(data[,2:ncol(data)],2,function(x){sum(x)/length(x)}),"\n")
