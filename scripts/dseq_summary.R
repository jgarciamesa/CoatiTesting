data = read.csv("data/dseq.csv")
summary = apply(data[,-1],2,sum)/nrow(data)
write.table(file="data/dseq_summary", x = summary, quote = FALSE, sep = "\t\t", col.names = FALSE)
