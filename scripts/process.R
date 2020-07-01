library(seqinr)
library(stringr)

omega = function(x) {
    o = kaks(x)
    ka = o$ka[1]
    ks = o$ks[1]
    if(ka < 0 || ks < 0) {
        return(0)
    }
    if(ka == 0 && ks == 0) {
        return(0)
    }
    if(ks == 0) {
        return(10)
    }
    ka/ks
}

ka = function(x) {
    if(class(x) != "alignment") {
        return(NA)
    }
    o = kaks(x)
    ka = o$ka[1]
    ka
}

ks = function(x) {
    if(class(x) != "alignment") {
        return(NA)
    }
    o = kaks(x)
    ks = o$ks[1]
    ks
}

tstat = function(x) {
    x = x$cleaned
    if(class(x) != "alignment") {
        return(NA)
    }
    o = kaks(x)
    tt = (o$ka[1]-o$ks[1])/sqrt(o$vka[1]+o$vks[1])
    tt
}

has_gaps = function(x) {
    x = x$raw
    if(class(x) != "alignment") {
        return(c(TRUE,TRUE))
    }
    x$seq %>% str_detect(pattern="[^acgtn]")
}

gap_lens = function(x) {
    x = x$raw
    if(class(x) != "alignment") {
        return(NA)
    }
    x$seq %>% str_extract_all(pattern="[^acgtn]+")
}

fasta = list.files("sim2",pattern="*.fasta")
fasta = list.files("fasta",pattern="*.rawcds.fasta")

process_aln = function(d,fasta) {
    ret = list()

    for(f in fasta) {
        x = try(read.alignment(str_c(d, "/", f),"fasta"),silent=TRUE)
        if(class(x) == "try-error") {
            break
        } else {
            seq = x$seq %>% str_split(pattern="")
            o = (seq[[which(x$nam == "homo_sapiens")]] %in% c("a","c","g","t","n"))
            seq = lapply(seq, "[", o)
            x$seq = lapply(seq, str_c, collapse="")
            ret[[f]] = x
        }
    }
    ret   
}

coati = process_aln("coati",fasta)
clustal = process_aln("clustaloaa",fasta)
prank = process_aln("prankcodon",fasta)
macse = process_aln("macse",fasta)
mafft = process_aln("mafft",fasta)
sim = process_aln("sim",fasta)


fasta = list.files("sim2",pattern="*.fasta")

process_fasta = function(d,fasta) {
    ret = list()

    for(f in fasta) {
        x = try(read.fasta(str_c(d, "/", f)),silent=TRUE)
        if(class(x) == "try-error") {
            ret[[f]] = list()
        } else {
            ret[[f]] = x
        }
    }
    ret   
}

coati = process_fasta("coati",fasta)
prankdna = process_fasta("prankdna",fasta)
prankcodon = process_fasta("prankcodon",fasta)
macse = process_fasta("macse",fasta)
mafft = process_fasta("mafft",fasta)


coati = process("coati")
prankdna = process("prankdna")
macse = process("macse")

coati_g = sapply(coati,has_gaps)
prankdna_g = sapply(prankdna,has_gaps)
prankcodon_g = sapply(prankcodon,has_gaps)
macse_g = sapply(macse,has_gaps)

g = apply(coati_g | macse_g | prankdna_g | prankcodon_g,2,any)

coati_ka = sapply(coati,ka)
prank_ka = sapply(prank,ka)
macse_ka = sapply(macse,ka)
mafft_ka = sapply(mafft,ka)
clustal_ka = sapply(clustal,ka)

coati_ks = sapply(coati,ks)
prankdna_ks = sapply(prankdna,ks)
prankcodon_ks = sapply(prankcodon,ks)
macse_ks = sapply(macse,ks)

coati_t = sapply(coati,tstat)
prankdna_t = sapply(prankdna,tstat)
prankcodon_t = sapply(prankcodon,tstat)
macse_t = sapply(macse,tstat)

coati_k = 2*coati_ka/(coati_ka+coati_ks)-1
prankdna_k = 2*prankdna_ka/(prankdna_ka+prankdna_ks)-1
prankcodon_k = 2*prankcodon_ka/(prankcodon_ka+prankcodon_ks)-1
macse_k = 2*macse_ka/(macse_ka+macse_ks)-1

plot(coati_ka[g],prankdna_ka[g],ylim=c(0,0.05),xlim=c(0,0.05))
points(coati_ka[g],macse_ka[g],col="blue")
points(coati_ka[g],prankcodon_ka[g],col="red")

plot(coati_ks,prankdna_ks,ylim=c(0,0.1),xlim=c(0,0.1))
points(coati_ks,macse_ks,col="blue")
points(coati_ks,prankcodon_ks,col="red")

mka = cbind(prankdna_ka,prankcodon_ka,macse_ka)-coati_ka
matplot(coati_ka[g],mka[g,],pch=19,xlim=c(0,0.03))
matplot(seq_len(sum(g)),mka[g,],pch=19)

m = cbind(coati_k,prankdna_k,macse_k,prankcodon_k)
m[is.nan(m)] = 0

coati_gaps = lapply(coati,gap_lens)
prankdna_gaps = lapply(prankdna,gap_lens)
macse_gaps = lapply(macse,gap_lens)
prankcodon_gaps = lapply(prankcodon,gap_lens)

gaps_ins = c(
    nchar(unlist(lapply(coati_gaps,"[[",1))),
    nchar(unlist(lapply(prankdna_gaps,"[[",1))),
    nchar(unlist(lapply(macse_gaps,"[[",1)))
)

gaps_del = c(
    nchar(unlist(lapply(coati_gaps,"[[",2))),
    nchar(unlist(lapply(prankdna_gaps,"[[",2))),
    nchar(unlist(lapply(macse_gaps,"[[",2)))
)

gaps_ins_count = c(
    sapply(lapply(coati_gaps,"[[",1),length),
    sapply(lapply(prankdna_gaps,"[[",1),length),
    sapply(lapply(macse_gaps,"[[",1),length)
)

gaps_del_count = c(
    sapply(lapply(coati_gaps,"[[",2),length),
    sapply(lapply(prankdna_gaps,"[[",2),length),
    sapply(lapply(macse_gaps,"[[",2),length)
)

## RENAME COATI SEQS
# for(f in fasta) {
#     f = str_c("coati/", f)
#     a = read.fasta(f,forceDNAtolower = FALSE)
#     n = names(a)
#     n[n == "input"] = "homo_sapiens"
#     n[n == "output"] = "pan_troglodytes"
#     write.fasta(a, n, f)
# }

fasta = list.files("sim",pattern="*.fasta")

process_fasta = function(d,fasta) {
    ret = list()

    for(f in fasta) {
        x = try(read.fasta(str_c(d, "/", f)),silent=TRUE)
        if(class(x) == "try-error") {
            ret[[f]] = list()
        } else {
            ret[[f]] = x
        }
    }
    ret   
}

coati = process_fasta("coati",fasta)
clustal = process_fasta("clustaloaa",fasta)
prank = process_fasta("prankcodon",fasta)
macse = process_fasta("macse",fasta)
mafft = process_fasta("mafft",fasta)
sim = process_fasta("sim",fasta)

# gaps count as 0
aln_dist = function(A,B) {
    if(length(A) != length(B)) {
        return(1)
    }

    h = function(hs,pt) {
        hs_b = tolower(hs) %in% letters
        pt_b = tolower(pt) %in% letters
        hs_o = rep(0,length(hs))
        pt_o = rep(0,length(pt))
        hs_o[hs_b] = seq(1,sum(hs_b))
        pt_o[pt_b] = seq(1,sum(pt_b))
        list("hs"=pt_o[hs_b],"pt"=hs_o[pt_b])
    }
    ha = h(A$homo_sapiens, A$pan_troglodytes)
    hb = h(B$homo_sapiens, B$pan_troglodytes)
    mean(c(ha$hs != hb$hs, ha$pt != hb$pt))
}

# gaps count using position
aln_dist = function(A,B) {
    if(length(A) != length(B)) {
        return(1)
    }

    h = function(hs,pt) {
        hs_b = tolower(hs) %in% letters
        pt_b = tolower(pt) %in% letters
        hs_o = rep(0,length(hs))
        pt_o = rep(0,length(pt))
        hs_o[hs_b] = seq(1,sum(hs_b))
        pt_o[pt_b] = seq(1,sum(pt_b))
        r = rle(!hs_b)
        s = cumsum(r$lengths)-r$lengths+1
        e = cumsum(r$lengths)
        for( o in which(r$values)) {
            hs_o[seq.int(s[o],e[o])] = ifelse(s[o]==1,0,-hs_o[s[o]-1])
        }

        r = rle(!pt_b)
        s = cumsum(r$lengths)-r$lengths+1
        e = cumsum(r$lengths)
        for( o in which(r$values)) {
            pt_o[seq.int(s[o],e[o])] = ifelse(s[o]==1,0,-pt_o[s[o]-1])
        }

        list("hs"=pt_o[hs_b],"pt"=hs_o[pt_b])
    }
    ha = h(A$homo_sapiens, A$pan_troglodytes)
    hb = h(B$homo_sapiens, B$pan_troglodytes)
    mean(c(ha$hs != hb$hs, ha$pt != hb$pt))
}


alignment_dist = function(A,B) {
    n = names(A)
    ret = sapply(n, function(x) aln_dist(A[[x]], B[[x]]))
    names(ret) = n
    ret
}


coati_dist = alignment_dist(sim,coati)
macse_dist = alignment_dist(sim,macse)
mafft_dist = alignment_dist(sim,mafft)
prank_dist = alignment_dist(sim,prank)
clustal_dist = alignment_dist(sim,clustal)

per_acgt = sapply(sim, function(x) {
    x = tolower(x$pan_troglodytes)
    sum(x %in% c('a','c','g','t'))/sum(x %in% letters)
})

mat_dist = cbind(coati_dist,macse_dist,mafft_dist,clustal_dist,prank_dist)

mat_qual = pmin(-10*log10(mat_dist),40)

apply(mat_dist,2,mean)

best = apply(mat_dist,1,min)

apply(mat_dist == best,2,sum)
