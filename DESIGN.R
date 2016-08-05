scanonevar(cross,
           mean.formula,
           var.formula)

validate.input

wrangle.input

do.scan




validate.input(pheno,
               geno,
               mean.formula,
               var.formula)

pheno and geno must imply the same number of individuals

all chrs to scan must be in geno

all vars from mean.formula must be in pheno or geno or be mean.QTL.add or mean.QTL.dom

all vars from var.formula must be in pheno or geno or be var.QTL.add or var.QTL.dom





wrangle.input(cross,
              mean.formula,
              var.formula)

genoprob.tbl <- make.genoprob.tbl(cross)

modeling.df <- data.frame(whichever QTL effects are present to modeling.df as NA)

a <- find vars from mean.formula that are in pheno
b <- find vars from var.formula that are in pheno
add unique(a, b) to modeling.df from pheno

c <- find vars from mean.formula that are in geno
d <- find vars from var.formula that are in geno
add unique(c, d) to modeling.df from geno

return modeling.df


make.genoprob.table(cross)

names(genoprob.tbl) == c(marker.name, chr, pos, iid, marker)
return(genoprob.tbl)