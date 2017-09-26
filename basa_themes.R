## theme for BASA trait paper plots
##function that applies theme to genotpyic variation plots
geno_theme<-function(){
  list(
    geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0),
    geom_point(aes(fill=SEX), shape=21),
    scale_shape_nell(),
    theme_nell(),
    scale_fill_manual(values=c('black','white')),
    theme(axis.text.x=element_blank(), axis.title.x=element_blank(), legend.position='none', axis.text.y = element_text(size = 10), 
          axis.title.y = element_text(size = 12), axis.ticks.x=element_blank(), plot.margin=margin(0.25,0,0.25,0,'cm'))
  )
}

##function for sex plots
sex_theme<-function(){
  list(
    geom_errorbar(aes(ymin=lsm-lsm_se, ymax=lsm+lsm_se), width=0), geom_point(aes(fill=SEX), size=1.85, shape=21),
    scale_fill_manual(values=c('black','white')),
    scale_shape_nell(),
    theme_nell(),
    theme(axis.text=element_blank(),axis.title=element_blank(), legend.position='none', plot.margin=margin(0.25,0,0.25,0,'cm'), axis.ticks=element_blank())
  )
}

##make genotypic plot, sex plot, plot together
geno_sex_plot<-function(tr8, ymin, ymax, ylabel){
  plot_grid(
    ggplot(data=lsm.final%>%filter(trait==paste(tr8)), aes(x=interaction(reorder(FAM,lsmean),SEX), y=lsmean))+labs(y=paste(tr8))+ylim(c(ymin,ymax))+geno_theme()+labs(y=paste(ylabel)),
    ggplot(data=lsm.sex%>%filter(trait==paste(tr8)), aes(x=SEX, y=lsm))+sex_theme()+ylim(c(ymin,ymax)), align='hv', nrow=1, ncol=2, rel_widths=c(1,.4))
}

geno_plot<-function(tr8, ymin, ymax, ylabel){
  ggplot(data=lsm.final%>%filter(trait==paste(tr8)), aes(x=interaction(reorder(FAM,lsmean),SEX), y=lsmean))+ylim(c(ymin,ymax))+geno_theme()+labs(y=ylabel)
}

sex_plot<-function(tr8, ymin, ymax, ylabel){
  ggplot(data=lsm.sex%>%filter(trait==paste(tr8)), aes(x=SEX, y=lsm))+sex_theme()+ylim(c(ymin,ymax))
}

#######################################################c
##functions for arthropod community data


##Supply functional group (arths, Pred, Her, Omni, Detrit, or Non), order of aggregation (taxa, Order, Family,Suborder, Superfamily,Guild??), distance, whether or not to log transform abundances
##returns list containing 
##at the individual plant level: full dataframe, community matrix, log transofrmed matrix, distances, 

arth_comm_ind<-function(trophic = 'arths', group = 'Order', dist = 'bray', log = FALSE, rel.abun = TRUE){
  arthropods$agg.level<-arthropods[,group]
  ##community dataframe by individual plant
  if (trophic == 'arths'){
    comm.df <- arthropods%>%dcast(ID~agg.level, value.var='total_arths', fun.aggregate=sum)%>%left_join(traits%>%dplyr::select(ID, FAM, SEX, WATER, BLOCK, ROW,COL,RGR, flowers_all, Per_H2O, MT_sum, SLA, leaf.toughness,CN)%>%transform(ID=as.character(ID)), by='ID')
  } else {
    comm.df <- arthropods%>%filter(Feeding == trophic)%>%dcast(ID~agg.level, value.var='total_arths', fun.aggregate=sum)%>%
    left_join(traits%>%dplyr::select(ID, FAM, SEX, WATER, BLOCK, ROW, COL, RGR, flowers_all, Per_H2O, MT_sum, SLA, leaf.toughness,CN)%>%transform(ID=as.character(ID)), by='ID')%>%dplyr::select(ID, FAM, SEX, WATER, BLOCK, ROW,COL, everything())
  }
  ##community matrix
  if (log == TRUE) {
    comm.mat<-as.data.frame(comm.df)%>%dplyr::select(-ID, -FAM, -SEX, -WATER, -BLOCK, -ROW, -COL, -RGR, -flowers_all, -Per_H2O, -MT_sum, -SLA, -leaf.toughness, -CN)%>%mutate_each(funs(log(.+2)))
  ##relative densities
  } else if (rel.abun == TRUE){
    comm.mat<-as.data.frame(comm.df)%>%dplyr::select(-ID, -FAM, -SEX, -WATER, -BLOCK, -ROW, -COL, -RGR, -flowers_all, -Per_H2O, -MT_sum, -SLA, -leaf.toughness, -CN)
    comm.mat<-comm.mat%>%mutate_all(funs(./rowSums(comm.mat,na.rm=TRUE)))
  } else {
    comm.mat<-as.data.frame(comm.df)%>%dplyr::select(-ID, -FAM, -SEX, -WATER, -BLOCK, -ROW, -COL, -RGR, -flowers_all, -Per_H2O, -MT_sum, -SLA, -leaf.toughness, -CN)
  }

  ##ditance matrix
  comm.dist<-vegdist(comm.mat, method = dist, na.rm=TRUE)
  
  ##returned list of dfs
  end.df<-list(comm.df = comm.df, comm.mat =comm.mat, comm.dist = comm.dist)
  ##
}

#at the genotype level: lsmean community df, relative density matrix, dissimilariity matrix

arth_comm_geno<-function(trophic = 'arths', group = 'Order', dist = 'bray', log = FALSE, rel.abun = TRUE){
  ind.comm<-arth_comm_ind(trophic=trophic, group=group, dist=dist, log=log)
  input.df<-ind.comm$comm.df
  
  ##calculate lsmeans for each arthopod grouping
  ord.list<-colnames(ind.comm$comm.mat)
  
  arth.ord.fam<-NULL
  for (i in 1:length(ord.list)){
    ordname<-ord.list[i]
    input.df$order<-input.df[,ordname]
    input.df$order<-ifelse(is.na(input.df$order), 0, input.df$order)
    ord.cc<-input.df
    trait.mod<-lmer(order~FAM+(1|ROW)+(1|WATER), data=ord.cc, na.action='na.omit')
    ord.lsm<-lsmeans::lsmeans(trait.mod, ~FAM)
    ord.df<-input.df%>%dplyr::select(SEX,FAM)%>%unique()%>%
      left_join(summary(ord.lsm), by='FAM')%>%mutate(taxa = paste(ord.list[i]))
    arth.ord.fam<-rbind(arth.ord.fam, ord.df)
  }
  ##cast df
  arth.ord.fam$lsmean<-ifelse(arth.ord.fam$lsmean < 0.05, 0, arth.ord.fam$lsmean)
  arth.cast<-arth.ord.fam%>%dcast(SEX+FAM~taxa, value.var='lsmean')
  
  ##by sex
  lsm.sex<-arth.ord.fam%>%group_by(taxa, SEX)%>%summarize(lsm = mean(lsmean), lsm_se=se(lsmean))
  
  ##community matrix
  if (log == TRUE) {
    comm.mat<-as.data.frame(arth.cast)%>%dplyr::select(-FAM, -SEX)%>%mutate_each(funs(log(.+1)))
  ##relative densities
  } else if (rel.abun == TRUE){
    comm.mat<-as.data.frame(arth.cast)%>%dplyr::select(-FAM, -SEX)
    comm.mat<-comm.mat%>%mutate_all(funs(./rowSums(comm.mat,na.rm=TRUE)))
  } else {
    comm.mat<-as.data.frame(comm.df)%>%dplyr::select(-FAM, -SEX)
  }
  
  ##dissimilarity
  ord.dist.rel<-vegdist(comm.mat, method=dist)
  
  
  output.df<-list(lsms=arth.ord.fam, comm.df = arth.cast, sex.df = lsm.sex, comm.mat = comm.mat, comm.dist = ord.dist.rel)
  return(output.df)
  
}


####################
##functions to produce dbrda plots

plot_cap_df<-function(dbrda, unc = TRUE, species = TRUE, data = arth.df$comm.df, scaly=3, fam.center=TRUE){
  #plot centroids for each genotype, unconstrained ordination
  plot.df.ind<-data.frame(SEX= data$SEX, ID=rownames(dbrda$CA$u),
                          FAM=data$FAM, 
                          MDS1 = scores(dbrda, scaling=scaly, correlation=TRUE)$sites[,1], 
                          MDS2 = scores(dbrda, scaling=scaly, correlation=TRUE)$sites[,2])
  
  ##axis labels with var explained
  eig<-eigenvals(dbrda)
  xlabel<-paste0('MDS1 (', round(100*cumsum(eig / sum(eig))[[1]],2), '% of total variation)')
  ylabel<-paste0('MDS2 (', round(100*(eig[[2]] / sum(eig)),2), '% of total variation)')
  
  ##use envfit for fitted order arrow significance and fam centroids
  pco.env<-envfit(dbrda, data%>%dplyr::select(everything()), na.rm=TRUE)
  arrow.df<-data.frame(order= rownames(pco.env$vectors$arrows), MDS1 = pco.env$vectors$arrows[,1], MDS2 = pco.env$vectors$arrows[,2], pval = pco.env$vectors$pvals, rsq = pco.env$vectors$r)%>%filter(pval <= 0.05)

  fam.c<-data.frame(FAM = rownames(pco.env$factors$centroids), MDS1 = pco.env$factors$centroids[,1], MDS2 = pco.env$factors$centroids[,2])
  fam.center<-fam.c%>%filter(grepl('FAM', FAM))%>%
    mutate(FAM =gsub('FAM','', FAM))%>%
    left_join(data%>%dplyr::select(FAM,SEX)%>%unique()%>%transform(FAM=as.character()))
  fam.center$SEX<-ifelse(fam.center$FAM > 20, 'F', 'M')
  
  sex.center<-fam.c%>%filter(FAM == 'SEXM' | FAM =='SEXF')
  
  ##this si scaling 2 (species axes correlations)
  plot.df.sp<-data.frame(order=rownames(scores(dbrda, scaling=scaly, correlation=TRUE)$species), 
                         MDS1 = scores(dbrda, scaling=scaly, correlation=TRUE)$species[,1], 
                         MDS2 = scores(dbrda, scaling=scaly, correlation=TRUE)$species[,2])%>%filter(order %in% arrow.df$order)
  
  df<-list(ind.df = plot.df.ind, fam.df = fam.center, sp.df=plot.df.sp, sp.env = arrow.df, ylab=ylabel, xlab=xlabel,pco.fit=pco.env, sex.center = sex.center)
  return(df)
}



###  unconstrained plot
##site is eitherplot_cap_df$ind.df or fam.df
unc_plot<-function(df, fam.df, ind.df, species, env, arrow.rel=.25){
  ##unconstrained ordination with jsut sites
  ##just points
  ind.plot<-ggplot(data=ind.df, aes(x=MDS1, y=MDS2))+
    geom_jitter(aes(fill=SEX),shape=21,size=2, alpha=.9)+
    theme_nell()+
    scale_fill_manual(values=c('black','white'))+
    theme(axis.line=element_blank())+
    geom_hline(aes(yintercept=0), lty='dashed')+
    geom_vline(aes(xintercept=0), lty='dashed')+
    labs(x=df$xlab, y=df$ylab)
  
  site.plot<-ggplot(fam.df, aes(x=MDS1, y=MDS2))+
    geom_jitter(aes(fill=SEX),shape=21,size=2, alpha=.9)+
    theme_nell()+scale_fill_manual(values=c('black','white'))+
    theme(axis.line=element_blank())+
    geom_hline(aes(yintercept=0), lty='dashed')+
    geom_vline(aes(xintercept=0), lty='dashed')+
    labs(x=df$xlab, y=df$ylab)

  ##with arrows
  site.sp.plot<-ggplot(fam.df, aes(x=MDS1, y=MDS2))+
    geom_jitter(aes(fill=SEX), shape=21,size=2, alpha=.9)+
    theme_nell()+
    scale_fill_manual(values=c('black','white'))+
    theme(axis.line=element_blank())+
    geom_hline(aes(yintercept=0), lty='dashed')+
    geom_vline(aes(xintercept=0), lty='dashed')+
    labs(x=df$xlab, y=df$ylab)+
    geom_segment(data = species, aes(x = 0, xend = MDS1*as.numeric(arrow.rel), y = 0, yend = MDS2*as.numeric(arrow.rel)),
                                           arrow = arrow(length = unit(0.2, "cm")), colour = "black", size=1)+
    geom_text(data=species, aes(x=MDS1*as.numeric(arrow.rel), y=MDS2*as.numeric(arrow.rel), label=order))
    
  site.env.plot<-ggplot(fam.df, aes(x=MDS1, y=MDS2))+
    geom_jitter(aes(fill=SEX), shape=21,size=2, alpha=.9)+
    theme_nell()+
    scale_fill_manual(values=c('black','white'))+
    theme(axis.line=element_blank())+
    geom_hline(aes(yintercept=0), lty='dashed')+
    geom_vline(aes(xintercept=0), lty='dashed')+
    labs(x=df$xlab, y=df$ylab)+
    geom_segment(data = env, aes(x = 0, xend = MDS1*as.numeric(rsq), y = 0, yend = MDS2*as.numeric(rsq)),
                 arrow = arrow(length = unit(0.2, "cm")), colour = "black", size=1)+
    geom_text(data=env, aes(x=MDS1*as.numeric(rsq), y=MDS2*as.numeric(rsq), label=order))


  output<-list(ind = ind.plot, sites = site.plot, site.sp = site.sp.plot, site.env =site.env.plot)
  
}




