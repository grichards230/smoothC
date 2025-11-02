source('Find_C_SmoothC.R')
library(sigmoid)
#' compare_models()
#' @param null.vars
#' @param added.vars
#' @param time 
#' @param status
#' @param data
#' @param weighted Boolean, default TRUE
#' @param dept.censor Boolean, default FALSE
#' @return Generates side-by-side plots comparing concordance gains by 
#' traditional C-index and proposed new Smooth C-index, as well as:
#'  * C.df a data frame of c-indices
compare_models = function(null.vars, added.vars, time, status, data,
                          weighted = TRUE, dept.censor = FALSE) {
  
    library(ggplot2)
    library(ggpubr)
    
    ## Calculate and store concordance values
    null.form = paste0(
        "Surv(",
        time,
        ", ",
        status,
        ") ~ ",
        paste(null.vars, collapse = " + ")
    )
    cvals.null = find_C_smoothC(
        coxph(as.formula(null.form), data = data),
        weighted = weighted,
        dept.censor = dept.censor
    )
    full.form = paste0(c(null.form, added.vars), collapse = " + ")
    cvals.full = find_C_smoothC(
        coxph(as.formula(full.form), data = data),
        weighted = weighted,
        dept.censor = dept.censor
    )
    
    c.idx1 = cvals.null$C
    c.idx2 = cvals.full$C
    nu1 = cvals.null$smooth.C
    nu2 = cvals.full$smooth.C
    
    C.df=data.frame(C=c(c.idx1,c.idx2),exp.nu=c(exp(nu1),exp(nu2)))
    
    df=data.frame(x=rep(c(-3,0,0,3),2),y=c(1-c.idx1,1-c.idx1,c.idx1,c.idx1,
                                           1-c.idx2,1-c.idx2,c.idx2,c.idx2),
                  type=rep(as.character(1:2),each=4),
                  x0=rep(as.character(c(1,1,2,2)),2))
    
    g=ggplot(df)
    g=g+geom_line(aes(x=x,y=y,linetype=type,colour=type,group=interaction(x0,type)),lwd=1.5)
    g=g+scale_colour_manual(labels=c("Initial Model","Full Model"),values=c("dodgerblue3","navy"))
    g=g+scale_linetype_manual(labels=c("Initial Model","Full Model"),values=c("dotted","solid"))
    g=g+theme_classic()
    g=g+theme(legend.position="top",legend.title=element_blank(),legend.key.width=unit(3,"cm"),text=element_text(size=24))
    g=g+xlab(expression(paste(italic(R[i]-R[j])," (SD)")))
    g=g+ylab(expression(italic(P(T[i] < T[j]))))
    g=g+scale_x_continuous(breaks=seq(-3,3,1))
    g=g+ylim(0,1)
    g=g+ggtitle("(a)")
  
    x.val = seq(-3, 3, by = 0.01)
    
    ## Generate y-values using sigmoid function
    y.val1 = sigmoid(nu1 * x.val)
    y.val2 = sigmoid(nu2 * x.val)
    
    df=data.frame(x=rep(x.val,2),y=c(y.val1,y.val2),nu=rep(as.character(1:2),each=length(x.val)))
    g2=ggplot(df)
    g2=g2+geom_line(aes(x=x,y=y,linetype=nu,colour=nu),lwd=1.5)
    g2=g2+scale_colour_manual(labels=c("Initial Model","Full Model"),values=c("dodgerblue3","navy"))
    g2=g2+scale_linetype_manual(labels=c("Initial Model","Full Model"),values=c("dotted","solid"))
    g2=g2+theme_classic()
    g2=g2+theme(legend.position="top",legend.key.width=unit(3,"cm"),text=element_text(size=24),legend.text=element_text(size=24),legend.title=element_blank())
    g2=g2+xlab(expression(paste(italic(R[i]-R[j])," (SD)")))
    g2=g2+ylab("")
    g2=g2+scale_x_continuous(breaks=seq(min(x.val),max(x.val),1))
    g2=g2+ylim(0,1)
    g2=g2+ggtitle("(b)")
    
    out=ggarrange(g,g2,nrow=1,common.legend=T)
    
    ## Finally, return data frame of values
    return(list(C.df,out))
}

