
setwd("/Volumes/apocrita/autoScratch/2016-10-11_Bombus_population_genomics/results/2017-01-17_CNV_presence")

## Read data in:
data.del <- read.table("lowRD.site_list.low_comp.regions_removed.length.counts.txt")

## Subset only count data:
data.del.subset <- data.del[, 4:26]

## Divide each column by column sum:
plot(colSums(data.del.subset))

autoplot(prcomp(as.data.frame(t(data.del.subset)/colSums(data.del.subset))))

data.del.subset.transformed.df <- as.data.frame(t(data.del.subset))

test <- table(round(((colSums(data.del.subset.transformed.df >= 3))/23), digits=1))

plot(test/nrow(data.del))

## Add a status identifier for each CNV:
count <- 0
test <- vector()

while (count < nrow(data.del)){
        count <- count + 1
        test <- c(test, paste("DEL", count, sep=""))
}

## Add the output to input:
data.del$identifier <- test

## Write output to file:
write.table(data.del, "lowRD.site_list.low_comp.regions_removed.plus_names.txt", col.names = F, row.names = F, sep="\t")
