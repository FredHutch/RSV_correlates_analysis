# @ by Youyi Fong
# start R in the correlates_report folder or make sure working directory is here
rm(list=ls())    
save.results.to="input"    
library(kyotil)
library(RSVcorr)
#dat.wide$trt=ifelse(dat.wide$trt=="Placebo", 0, 1)
trt.labels=c("Placebo","Vaccine")
dat.wide.s=as.data.frame(subset(dat.wide, sampled=="Y"))
markers=c(t(outer(assays,"."%.%times,paste0))) # first 8 are RSV markers


# strata sizes
tab=NULL
for (a in 0:1) {
    tmp=with(subset(dat.wide, trt==a), table(sitegrp, vacc2birthLESSTHAN30, y1))    
    tab=cbind(tab, cbind(tmp[,2:1,1], tmp[,2:1,2]))
}
sum(tab[,c(3,4,7,8)])# number of cases
colnames(tab)=rep(c("<30d",">=30d"),4); tab
#            Plac                  Vacc
#      cntl       case       cntl       case
#   <30d >=30d <30d >=30d <30d >=30d <30d >=30d
#NN   39   371    3    19   62   761    1    26
#PP   12    79    1    11   30   155    2     7
#QQ   15    82    2     6   39   170    5     5
#RR   13    71    1     7   23   154    1     3
#SS   11    40    0     1   18    89    0     2


# strata sizes for endpoint 2 
tab=NULL
for (a in 0:1) {
    tmp=with(subset(dat.wide, trt==a), table(sitegrp, vacc2birthLESSTHAN30, y2))    
    tab=cbind(tab, cbind(tmp[,2:1,1], tmp[,2:1,2]))
}
sum(tab[,c(3,4,7,8)])# number of cases
colnames(tab)=rep(c("<30d",">=30d"),4); tab
#            Plac                  Vacc
#      cntl       case       cntl       case
#   <30d >=30d <30d >=30d <30d >=30d <30d >=30d
#NN   40   378    2    12   63   778    0     9
#PP   13    85    0     5   31   162    1     0
#QQ   16    86    1     2   43   173    1     2
#RR   14    73    0     5   24   156    0     1
#SS   11    41    0     0   18    91    0     0


# VISC samples distribution by strata
# remove 10 subjects in the vaccine arm that miss either d0 or d14/cord
miss.d0.rsv=      apply(is.na(dat.wide.s[dat.wide.s$trt==1, c("RSVB.d0"),drop=F]),              1, all)
miss.d14.cord.rsv=apply(is.na(dat.wide.s[dat.wide.s$trt==1, c("RSVB.d14","RSVB.cord"),drop=F]), 1, all)
dat.wide.a=rbind(subset(dat.wide.s, trt==1)[!(miss.d0.rsv | miss.d14.cord.rsv),], subset(dat.wide.s, trt==0))
tab.a=NULL
for (a in 0:1) {
    tmp=with(subset(dat.wide.a, trt==a), table(sitegrp, vacc2birthLESSTHAN30, y1))    
    tab.a=cbind(tab.a, cbind(tmp[,2:1,1], tmp[,2:1,2]))
}
colnames(tab.a)=rep(c("<30d",">=30d"),4); tab.a



# boxplots for RSVA and RSVB
ylim=range(subset(dat.wide, select=c(RSVA.d0, RSVA.d14, RSVA.cord)), na.rm=T)
myfigure(mfcol=c(2,2), oma=c(0,2,0,0))
for (trt. in c(0,1)) {
    myboxplot (subset(dat.wide, trt==trt. & y1==0, select=c(RSVA.d0, RSVA.d14, RSVA.cord)), main=trt.labels[trt.+1], log="y", ylim=ylim)
    myboxplot (subset(dat.wide, trt==trt. & y1==0, select=c(RSVB.d0, RSVB.d14, RSVB.cord)), main=trt.labels[trt.+1], log="y", ylim=ylim)
}
mtext(assays[1:2], side = 2, line = 0, outer = T, at = c(3,1)/4)
mydev.off(file="input/boxplot_rsv_by_time_trt")

# boxplots for EIA and PCA
ylim.1=range(subset(dat.wide, select=c(EIA.d0, EIA.d14, EIA.cord)), na.rm=T)
ylim.2=range(subset(dat.wide, select=c(PCA.d0, PCA.d14, PCA.cord)), na.rm=T)
myfigure(mfcol=c(2,2), oma=c(0,2,0,0))
for (trt. in c(0,1)) {
    myboxplot (subset(dat.wide, trt==trt. & y1==0, select=c(EIA.d0, EIA.d14, EIA.cord)), main=trt.labels[trt.+1], log="y", ylim=ylim.1)
    myboxplot (subset(dat.wide, trt==trt. & y1==0, select=c(PCA.d0, PCA.d14, PCA.cord)), main=trt.labels[trt.+1], log="y", ylim=ylim.2)
}
mtext(assays[3:4], side = 2, line = 0, outer = T, at = c(3,1)/4)
mydev.off(file="input/boxplot_eia_pca_by_time_trt")


#boxplots for fold changes

for (a in assays) dat.wide[[paste0(a,".log10cordoverd0")]]=log10(dat.wide[[paste0(a,".cord")]]/dat.wide[[paste0(a,".d0")]])

myfigure(mfcol=c(2,2), oma=c(0,2,0,0))
for (trt. in c(0,1)) {
    ylim=range(subset(dat.wide, select=c(RSVA.log10d14overd0, RSVA.log10cordoverd0)), na.rm=T)
    myboxplot (subset(dat.wide, trt==trt. & y1==0, select=c(RSVA.log10d14overd0, RSVA.log10cordoverd0)), main=trt.labels[trt.+1], log="", ylim=ylim, cex.axis=.8)
    myboxplot (subset(dat.wide, trt==trt. & y1==0, select=c(RSVB.log10d14overd0, RSVB.log10cordoverd0)), main=trt.labels[trt.+1], log="", ylim=ylim, cex.axis=.8)
}
mtext(assays[1:2], side = 2, line = 0, outer = T, at = c(3,1)/4)
mydev.off(file="input/boxplot_overd0_rsv")

myfigure(mfcol=c(2,2), oma=c(0,2,0,0))
for (trt. in c(0,1)) {
    ylim=range(subset(dat.wide, select=c(EIA.log10d14overd0, PCA.log10cordoverd0)), na.rm=T)
    myboxplot (subset(dat.wide, trt==trt. & y1==0, select=c(EIA.log10d14overd0, EIA.log10cordoverd0)), main=trt.labels[trt.+1], log="", ylim=ylim, cex.axis=.8)
    myboxplot (subset(dat.wide, trt==trt. & y1==0, select=c(PCA.log10d14overd0, PCA.log10cordoverd0)), main=trt.labels[trt.+1], log="", ylim=ylim, cex.axis=.8)
}
mtext(assays[3:4], side = 2, line = -1, outer = T, at = c(3,1)/4)
mydev.off(file="input/boxplot_overd0_eia_pca")



# scatterplots between time points
for (t in assays) {
    ylim=range(dat.wide[,c(t%.%c(".d0", ".d14", ".cord"))], na.rm=T)
    myfigure(mfrow=c(2,3), oma=c(0,2,0,0))
    for (trt. in c(0,1)) {
        corplot (as.formula(paste0(t,".d14~",t,".d0")), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=ylim, xlim=ylim, method="spearman")
        corplot (as.formula(paste0(t,".cord~",t,".d0")), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=ylim, xlim=ylim, method="spearman")
        corplot (as.formula(paste0(t,".cord~",t,".d14")), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=ylim, xlim=ylim, method="spearman")
    }
    mydev.off(file=paste0("input/scatterplot_",t))
}


# scatterplots between D14 and fold change
myfigure(mfrow=c(2,2), oma=c(0,2,0,0))
for (t in assays[1:2]) {
    ylim=range(dat.wide[,c(t%.%c(".d14", ".d14overd0"))], na.rm=T)
    for (trt. in c(0,1)) {
        corplot (as.formula(paste0(t,".d14overd0~",t,".d14")), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=ylim, xlim=ylim, method="spearman")
    }
}
mydev.off(file=paste0("input/scatterplot_rsvAB_d14_foldchange"))


# comparing how well D14 is correlated with birth level, vs. how well fold-rise is correlated with birth level.  
# If fold-rise were more strongly correlated with birth level, it could help explain why fold-rise was a stronger correlate.
tab=sapply (assays, function(a) {
    c(
        cor(subset(dat.wide, trt==1)[[paste0(a,".d14")]], subset(dat.wide, trt==1)[[paste0(a,".cord")]], method="spearman", use="pairwise"),
        cor(subset(dat.wide, trt==1)[[paste0(a,".d14overd0")]], subset(dat.wide, trt==1)[[paste0(a,".cord")]], method="spearman", use="pairwise")
    )
})
rownames(tab)=c("D14-Cord","D14overD0-Cord")
tab


myfigure(mfrow=c(2,2), oma=c(0,2,0,0))
for (t in assays[2+1:2]) {
    ylim=range(dat.wide[,c(t%.%c(".d14", ".d14overd0"))], na.rm=T)
    for (trt. in c(0,1)) {
        corplot (as.formula(paste0(t,".d14overd0~",t,".d14")), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=ylim, xlim=ylim, method="spearman")
    }
}
mydev.off(file=paste0("input/scatterplot_EIApca_d14_foldchange"))

# restrict to four fold change or above
myfigure(mfrow=c(2,2), oma=c(0,2,0,0))
for (t in assays[2+1:2]) {
    ylim=range(dat.wide[,c(t%.%c(".d14", ".d14overd0"))], na.rm=T)
    for (trt. in c(0,1)) {
        corplot (as.formula(paste0(t,".d14overd0~",t,".d14")), dat.wide[dat.wide$trt==trt. & dat.wide$y1==0 & dat.wide[[paste0(t,".d14overd0")]]>4,], main=trt.labels[trt.+1], log="xy", ylim=ylim, xlim=ylim, method="spearman")
    }
}
mydev.off(file=paste0("input/scatterplot_EIApca_d14_foldchange_4"))


# scatterplots between RSV A and B
ylim=range(subset(dat.wide, select=c(RSVA.d0, RSVA.d14, RSVA.cord)), na.rm=T)
myfigure(mfrow=c(2,3), oma=c(0,2,0,0))
for (trt. in c(0,1)) {
for (t in times[1:3]) {
    corplot (as.formula(paste0("RSVA.",t,"~RSVB.",t)), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=ylim, xlim=ylim, method="spearman")
}
}
mydev.off(file="input/scatterplot_rsva_rsvb")


# scatterplots between EIA and PCA
myfigure(mfrow=c(2,3), oma=c(0,2,0,0))
for (trt. in c(0,1)) {
for (t in times[1:3]) {
    corplot (as.formula(paste0("EIA.",t,"~PCA.",t)), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=NULL, xlim=NULL, method="spearman", add.diagonal.line=F)
}
}
mydev.off(file="input/scatterplot_eia_pca")


# scatterplots between RSVA and EIA
myfigure(mfrow=c(2,3), oma=c(0,2,0,0))
for (trt. in c(0,1)) {
for (t in times[1:3]) {
    corplot (as.formula(paste0("RSVA.",t,"~EIA.",t)), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=NULL, xlim=NULL, method="spearman", add.diagonal.line=F)
}
}
mydev.off(file="input/scatterplot_rsva_eia")


# scatterplots between RSVA and PCA
myfigure(mfrow=c(2,3), oma=c(0,2,0,0))
for (trt. in c(0,1)) {
for (t in times[1:3]) {
    corplot (as.formula(paste0("RSVA.",t,"~PCA.",t)), subset(dat.wide, trt==trt. & y1==0), main=trt.labels[trt.+1], log="xy", ylim=NULL, xlim=NULL, method="spearman", add.diagonal.line=F)
}
}
mydev.off(file="input/scatterplot_rsva_pca")


# make a table for correlation between assays
tab=
sapply(assays, function(a) {
sapply(assays, function(b) {
    cors=sapply (c(0,1), function(trt.) {
        sapply (times[1:3], function (t) {
            #myprint(a,b,trt.,t)
            cor (subset(dat.wide, trt==trt. & y1==0)[,paste0(a,".",t)], subset(dat.wide, trt==trt. & y1==0)[,paste0(b,".",t)], method="spearman", use="pairwise.c")
        })
    })
    mean(cors)
})
})
mytex(tab, file="tables/spearcor_assays", digit=2, input.folder=save.results.to)


# make a table for correlation between d0 and d14
res=
sapply(assays, function(a) {
sapply (c(0,1), function(trt.) {
    dat.tmp=subset(dat.wide, trt==trt. & y1==0)
    c(
        d14=cor(dat.tmp[[paste0(a,".d0")]], dat.tmp[[paste0(a,".d14")]], method="spearman", use="pairwise.c"),
        foldchange=cor(dat.tmp[[paste0(a,".d0")]], dat.tmp[[paste0(a,".log10d14overd0")]], method="spearman", use="pairwise.c")
    )
})
})
tab=res[c(1,3),]
rownames(tab)=trt.labels
mytex(tab, file="tables/spearcor_time_d14", digit=2, input.folder=save.results.to)
tab=res[c(2,4),]
rownames(tab)=trt.labels
mytex(tab, file="tables/spearcor_time_foldchange", digit=2, input.folder=save.results.to)




###################
# missing data 

# RSVA.d14 and RSB.d14 are either both present or both absent
all(with(dat.wide, !xor(is.na(RSVA.d14), is.na(RSVB.d14))))

tab.1=with(subset(dat.wide, sampled=='Y'), table(is.na(RSVA.d14), trt, y1))
tab.2=with(subset(dat.wide, sampled=='N'), table(is.na(RSVA.d14), trt)) # no cases
tab=cbind(tab.1[,,2], tab.1[,,1], tab.2)[2:1,]
colnames(tab)=rep(trt.labels,3)
tab=tab[2:1,]
rownames(tab)=c("Present","Missing")
mytex(tab, file="tables/samples_distr", col.headers = "\\hline\n  &  \\multicolumn{2}{c}{VISC Cases}  &  \\multicolumn{2}{c}{VISC Controls} &  \\multicolumn{2}{c}{Non-VISC Controls}  \\\\  \n", input.folder=save.results.to)

# only 2 out of two thousand have different status between RSVA.d0 and RSVA.d14
with(dat.wide, table(is.na(RSVA.d0), is.na(RSVB.d0)))

# 
with(dat.wide, table(is.na(RSVA.d0), is.na(RSVA.d14), sampled))
with(dat.wide, table(is.na(RSVA.d0), is.na(RSVA.cord), sampled))




# strata sizes
tab=NULL
for (a in 0:1) {
    tmp=with(subset(dat.wide, trt==a), table(sitegrp, vacc2birthLESSTHAN30, y1))    
    tab=cbind(tab, cbind(tmp[,2:1,1], tmp[,2:1,2]))
}
sum(tab[,c(3,4,7,8)])# number of cases
colnames(tab)=rep(c("<30d",">=30d"),4); tab
mytex(tab, file="tables/strata_sizes", input.folder=save.results.to, align=c("c","c","c","c","c|","c","c","c","c"),
    col.headers = "\\hline\n  & \\multicolumn{4}{c|}{Placebo} & \\multicolumn{4}{c}{Vaccine} \\\\\n  
                              & \\multicolumn{2}{c}{Control} & \\multicolumn{2}{c|}{Case}  & \\multicolumn{2}{c}{Control} & \\multicolumn{2}{c}{Case} \\\\\n")

# VISC samples distribution by strata
tab=NULL
for (a in 0:1) {
    tmp=with(subset(dat.wide, trt==a & sampled=="Y"), table(sitegrp, vacc2birthLESSTHAN30, y1))    
    tab=cbind(tab, cbind(tmp[,2:1,1], tmp[,2:1,2]))
}
colnames(tab)=rep(c("<30d",">=30d"),4); tab
mytex(tab, file="tables/strata_sizes_VISC", input.folder=save.results.to, align=c("c","c","c","c","c|","c","c","c","c"),
    col.headers = "\\hline\n  & \\multicolumn{4}{c|}{Placebo} & \\multicolumn{4}{c}{Vaccine} \\\\\n  
                              & \\multicolumn{2}{c}{Control} & \\multicolumn{2}{c|}{Case}  & \\multicolumn{2}{c}{Control} & \\multicolumn{2}{c}{Case} \\\\\n")

# VISC samples distribution by strata
tab=NULL
for (a in 0:1) {
    tmp=with(subset(dat.wide, trt==a & sampled=="N" & !is.na(RSVA.d14)), table(sitegrp, vacc2birthLESSTHAN30))    
    tab=cbind(tab, tmp[,2:1])
}
colnames(tab)=rep(c("<30d",">=30d"),2); tab
mytex(tab, file="tables/strata_sizes_nonVISC_hasRSVAd14", input.folder=save.results.to, align=c("c","c","c|","c","c"),
    col.headers = "\\hline\n  & \\multicolumn{2}{c|}{Placebo} & \\multicolumn{2}{c}{Vaccine} \\\\\n  
                              & \\multicolumn{2}{c|}{Control} & \\multicolumn{2}{c}{Control} \\\\\n")

tab=NULL
for (a in 0:1) {
    tmp=with(subset(dat.wide, trt==a & sampled=="N" & !is.na(RSVA.d0)), table(sitegrp, vacc2birthLESSTHAN30))    
    tab=cbind(tab, tmp[,2:1])
}
colnames(tab)=rep(c("<30d",">=30d"),2); tab
mytex(tab, file="tables/strata_sizes_nonVISC_hasRSVAd0", input.folder=save.results.to, align=c("c","c","c|","c","c"),
    col.headers = "\\hline\n  & \\multicolumn{2}{c|}{Placebo} & \\multicolumn{2}{c}{Vaccine} \\\\\n  
                              & \\multicolumn{2}{c|}{Control} & \\multicolumn{2}{c}{Control} \\\\\n")

tab=NULL
for (a in 0:1) {
    tmp=with(subset(dat.wide, trt==a & sampled=="N" & !is.na(RSVA.cord)), table(sitegrp, vacc2birthLESSTHAN30))    
    tab=cbind(tab, tmp[,2:1])
}
colnames(tab)=rep(c("<30d",">=30d"),2); tab
mytex(tab, file="tables/strata_sizes_nonVISC_hasRSVAcord", input.folder=save.results.to, align=c("c","c","c|","c","c"),
    col.headers = "\\hline\n  & \\multicolumn{2}{c|}{Placebo} & \\multicolumn{2}{c}{Vaccine} \\\\\n  
                              & \\multicolumn{2}{c|}{Control} & \\multicolumn{2}{c}{Control} \\\\\n")



#####################################################
# immune biomarkers missingness pattern


# All rows have at one of 16 markers
which(apply(is.na(dat.wide.s[,markers]), 1, all)) 

with (dat.wide.s, table(is.na(RSVA.d0), is.na(RSVB.d0)))
with (dat.wide.s, table(is.na(RSVA.d14), is.na(RSVB.d14)))
with (dat.wide.s, table(is.na(RSVA.cord), is.na(RSVB.cord)))


# number of missingness for each marker
tab=matrix(colSums(is.na(dat.wide.s[,markers])), ncol=4, byrow=T, dimnames=list(assays, times)); tab
mytex(tab, file="tables/missing_cnt", digit=0, input.folder=save.results.to)

mis.row=rowSums(is.na(dat.wide.s[,markers]))
table(mis.row)

# these have no cord blood markers
dat.wide.s[mis.row==4,markers]    
# these have no day14 markers
dat.wide.s[mis.row==8,markers]    
# one has only d14 markers, one has only d0 markers, three have only cord blood markers
dat.wide.s[mis.row==12,markers]    


# missinging
tmp.assays=c("RSVB")#,"RSVA")
tmp=t(outer(tmp.assays,"."%.%c("d0"),paste0)); miss.d0.rsv=apply(is.na(dat.wide.s[, tmp,drop=F]), 1, all)
tmp=t(outer(tmp.assays,"."%.%c("d14"),paste0)); miss.d14.rsv=apply(is.na(dat.wide.s[, tmp,drop=F]), 1, all)
tmp=t(outer(tmp.assays,"."%.%c("cord"),paste0)); miss.cord.rsv=apply(is.na(dat.wide.s[, tmp,drop=F]), 1, all)
tmp=t(outer(tmp.assays,"."%.%c("d14","cord"),paste0)); miss.d14.cord.rsv=apply(is.na(dat.wide.s[, tmp,drop=F]), 1, all)
tmp=t(outer(tmp.assays,"."%.%c("d14","cord","d0"),paste0)); miss.all.rsv=apply(is.na(dat.wide.s[, tmp,drop=F]), 1, all)

sum(miss.d0.rsv)
sum(miss.d14.rsv)
sum(miss.cord.rsv)
sum(miss.d14.cord.rsv)
sum(miss.all.rsv)

tab=table(miss.d14.rsv, miss.cord.rsv)
names(dimnames(tab))=c("miss d14", "miss cord")
mytex(tab, file="tables/miss_d14_cor", digit=0, input.folder=save.results.to)


# restrict to one arm
tmp.assays=c("RSVB")#,"RSVA")
tab=NULL
for (t in 1:0) {
    tmp=t(outer(tmp.assays,"."%.%c("d0"),paste0));         miss.d0.rsv=      apply(is.na(dat.wide.s[dat.wide.s$trt==t, tmp,drop=F]), 1, all)
    tmp=t(outer(tmp.assays,"."%.%c("d14","cord"),paste0)); miss.d14.cord.rsv=apply(is.na(dat.wide.s[dat.wide.s$trt==t, tmp,drop=F]), 1, all)
    res=table(miss.d0.rsv, miss.d14.cord.rsv)
    tab=cbind(tab, res)
}
tab
names(dimnames(tab))=c("miss d0", "miss d14/cord")
mytex(tab, file="tables/miss_d0_d14cord", input.folder=save.results.to, digit=0,     col.headers = "\\hline\n  miss d0&  \\multicolumn{3}{c}{miss d14/cord}  \\\\  \n &  \\multicolumn{2}{c}{Vaccine}  &  \\multicolumn{1}{c}{Placebo}  \\\\  \n")

# restrict to controls. the result is that 2 that miss both d0 and d14/cord are cases, all others are controls. 
tmp.assays=c("RSVB")#,"RSVA")
tab=NULL
for (t in 1:0) {
    myprint(t)
    tmp=t(outer(tmp.assays,"."%.%c("d0"),paste0));         miss.d0.rsv=      apply(is.na(dat.wide.s[dat.wide.s$trt==t & dat.wide.s$y1==0, tmp,drop=F]), 1, all)
    tmp=t(outer(tmp.assays,"."%.%c("d14","cord"),paste0)); miss.d14.cord.rsv=apply(is.na(dat.wide.s[dat.wide.s$trt==t & dat.wide.s$y1==0, tmp,drop=F]), 1, all)
    res=table(miss.d0.rsv, miss.d14.cord.rsv); print(res)
    tab=cbind(tab, res)
}
tab

# 
with(dat.wide.s, table(y1, y2))
with(dat.wide.s, table(y2, y3))




#################

with(dat.wide, table(is.na(EIA.d0), is.na(PCA.d0)))
with(dat.wide, table(is.na(EIA.d14), is.na(PCA.d14)))
with(dat.wide, table(is.na(EIA.cord), is.na(PCA.cord)))

with(dat.wide, table(is.na(EIA.d14), is.na(EIA.cord), y1))




tmp.assays=c("EIA")
tab=NULL
for (t in 1:0) {
    tmp=t(outer(tmp.assays,"."%.%c("d0"),paste0));         miss.d0.rsv=      apply(is.na(dat.wide[dat.wide$trt==t, tmp,drop=F]), 1, all)
    tmp=t(outer(tmp.assays,"."%.%c("d14","cord"),paste0)); miss.d14.cord.rsv=apply(is.na(dat.wide[dat.wide$trt==t, tmp,drop=F]), 1, all)
    res=table(miss.d0.rsv, miss.d14.cord.rsv)
    #print(res)
    tab=cbind(tab, res)
}
tab

names(dimnames(tab))=c("miss d0", "miss d14/cord")
mytex(tab, file="tables/miss_d0_d14cord_EIA", input.folder=save.results.to, digit=0,     col.headers = "\\hline\n  miss d0&  \\multicolumn{4}{c}{miss d14/cord}  \\\\  \n &  \\multicolumn{2}{c}{Vaccine}  &  \\multicolumn{2}{c}{Placebo}  \\\\  \n")


tmp.assays=c("PCA")
tab=NULL
for (t in 1:0) {
    tmp=t(outer(tmp.assays,"."%.%c("d0"),paste0));         miss.d0.rsv=      apply(is.na(dat.wide[dat.wide$trt==t, tmp,drop=F]), 1, all)
    tmp=t(outer(tmp.assays,"."%.%c("d14","cord"),paste0)); miss.d14.cord.rsv=apply(is.na(dat.wide[dat.wide$trt==t, tmp,drop=F]), 1, all)
    res=table(miss.d0.rsv, miss.d14.cord.rsv)
    #print(res)
    tab=cbind(tab, res)
}
tab

names(dimnames(tab))=c("miss d0", "miss d14/cord")
mytex(tab, file="tables/miss_d0_d14cord_PCA", input.folder=save.results.to, digit=0,     col.headers = "\\hline\n  miss d0&  \\multicolumn{4}{c}{miss d14/cord}  \\\\  \n &  \\multicolumn{2}{c}{Vaccine}  &  \\multicolumn{2}{c}{Placebo}  \\\\  \n")


subset(dat.wide, trt==0 & EIA.log10d14>3.7)
# there are five subjects, all are Y in all flags including ppefmfl ppefifl ppimifl ppimmfl itteffl ittimfl, except 1 subject,
# who has NA baseline EIA or PCA or RSVA/B values and is not in ppimmfl
