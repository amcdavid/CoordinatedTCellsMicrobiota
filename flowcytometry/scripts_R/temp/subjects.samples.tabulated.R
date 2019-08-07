subjects <- list(tphe = data.frame(id = unlist(lapply(rownames(ist.md$tphe), function(i) (strsplit(i, "_")[[1]][1]))),
                                  visit = unlist(lapply(rownames(ist.md$tphe), function(i) (strsplit(i, "_")[[1]][2])))),
                 ics = data.frame(id = unlist(lapply(rownames(ist.md$ics), function(i) (strsplit(i, "_")[[1]][1]))),
                                 visit = unlist(lapply(rownames(ist.md$ics), function(i) (strsplit(i, "_")[[1]][2])))),
                 tphe_ics = data.frame(id = unlist(lapply(intersect(rownames(ist.md$tphe), rownames(ist.md$ics)), function(i) (strsplit(i, "_")[[1]][1]))),
                                      visit = unlist(lapply(intersect(rownames(ist.md$tphe), rownames(ist.md$ics)), function(i) (strsplit(i, "_")[[1]][2])))))

table(table(subjects$tphe$id))
table(table(subjects$ics$id))
table(table(subjects$tphe_ics$id))

table(subjects$tphe$visit)
table(subjects$ics$visit)
table(subjects$tphe_ics$visit)

sapply(c(1, 7, 19), function(i) length(which(subjects$tphe$visit == i)))
sapply(c(1, 7, 19), function(i) length(which(subjects$ics$visit == i)))
sapply(c(1, 7, 19), function(i) length(which(subjects$tphe_ics$visit == i)))

i <- subjects$tphe_ics

length(grep("RPRC053", i[["id"]][i[["visit"]] == 1], invert = TRUE))
length(grep("RPRC053", i[["id"]][i[["visit"]] == 1]))
length(grep("RPRC053", i[["id"]][i[["visit"]] == 7], invert = TRUE))
length(grep("RPRC053", i[["id"]][i[["visit"]] == 7]))
length(grep("RPRC053", i[["id"]][i[["visit"]] == 19], invert = TRUE))
length(grep("RPRC053", i[["id"]][i[["visit"]] == 19]))
