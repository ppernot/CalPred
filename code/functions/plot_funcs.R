# Load packages and define parameters

if(!require(aplpack)) # Convex hulls
  install.packages("aplpack",dependencies=TRUE)
library(aplpack)

# Define transparent colors
col2tr = function(col,alpha) 
  rgb(t(col2rgb(col)),alpha=alpha,maxColorValue=255)
orchid_tr=col2tr("orchid",20)
blue_tr  =col2tr("blue"  ,20)
yellow_tr=col2tr("yellow"  ,120)
orange_tr=col2tr("orange"  ,120)

cexPlot = 4    # Magnification factor of text and symbols
reso    = 1200 # Image width in pixels

# Printing functions ####
prettyUnc = function(y, uy, numDig=1) {
  # Print result + uncertainty in parenthesis format 
  
  if (uy <= 0) return (y)
  
  # Get scales
  n0=1+floor(log10(abs(y)))
  n1=floor(log10(uy))
  
  # Format uncertainty
  switch(sign(n1) + 2, # Map (-1,0,1) to (1,2,3)
         {
           fmt = paste0("%",n0-n1+numDig-1,".",-n1+numDig-1,"f")
           short_uy = signif(uy/10^(n1-numDig+1),numDig)
         },
         {
           fmt = paste0("%",n0-n1+numDig-1,".",-n1+numDig-1,"f")
           short_uy = signif(uy/10^n1,numDig)
         },
         {    
           fmt = paste0("%",n0,".0f")
           short_uy = signif(uy/10^(n1-numDig+1),numDig)
         }
  )
  short_y  = sprintf(fmt, y)
  
  return( paste0(short_y,'(',short_uy,')') )
}

# Plot functions ####
plot_hull <- function(xlim=NULL, ylim=NULL, lwd= 5,
                      meth_type=NULL, true_point=NULL) {
  # Plot 95% contours of LJ parameters samples
  
  par(mfrow=c(1,1),mar=c(3,3,1.6,.2),mgp=c(2,.75,0),
      pty='s',tcl=-0.5, cex=cexPlot)
  plot(data$x,data$x,type='n',
       xlim=xlim,ylim=ylim,
       xlab = expression(paste(epsilon[LJ]," / K")),
       ylab = expression(paste(sigma[LJ]," / ",ring(A)))
  )
  if(!is.null(true_point))
    points(true_point[1],true_point[2],pch=8,cex=1,lwd=2,col=2)

  if(meth_type == 'hier') {
    cols= rep(1:5,2)
    pch = c(rep(16,5),rep(17,5))
    lty = c(rep(1,5),rep(2,5))

    legend('topright',ncol=4,cex=0.75,
           bty='o', bg=NULL, box.col=0,
           legend=c(parse(text=paste0('italic(D)^(',1:max(data$tag),')')),
                    'Hyper.'),
           col=c(cols[1:max(data$tag)],'orchid'),
           lty=c(lty[1:max(data$tag)],1),
           pch=c(rep(-1,max(data$tag)),-1),
           lwd=c(rep(lwd,max(data$tag)),rep(lwd,1)),
           seg.len=1.2)

    # Group params
    mu = extract(fit,"mu")[[1]]
    for(iD in 1:max(data$tag)) {
      eps = mu[,1,iD]
      sig = mu[,2,iD]
      P=compute.bagplot(eps, sig, factor =1.96)
      poly = P$hull.loop
      polygon(poly,col=NA,border=cols[iD],lwd=lwd,lty=lty[iD])
    }
  }

  if(meth_type == 'hier' | meth_type == 'abc') {
    mup = extract(fit,"mup")[[1]]
    eps = mup[,1]
    sig = mup[,2]
  } else {
    eps = extract(fit,"eps")[[1]]
    sig = extract(fit,"sig")[[1]]
  }

  P=compute.bagplot(eps, sig, factor =1.96)
  poly = P$hull.loop
  polygon(poly,col=NA,border='orchid',lwd=lwd,lty=1)
  points(mean(eps),mean(sig),pch=3,cex=1,lwd=4,col='orchid')
  grid(); box()

}
plot_predict <- function(center=FALSE, ylim=NULL, 
                         meth_type=NULL, legend=FALSE,
                         probs=c(0.025,0.975)) {
  # Plot 95% prediction intervals in data or residuals space
  resid  = extract(fit,"resid" )[[1]]
  y_conf = extract(fit,"y_conf")[[1]]
  y_pred = extract(fit,"y_pred")[[1]]
  
  q_conf = q_pred = matrix(NA,ncol=2,nrow=length(data$x_pred))
  u_conf = rep(NA,length(data$x_pred))
  for (i in 1:length(data$x_pred)) {
    u_conf[i]  = sd(y_conf[,i])
    q_conf[i,] = quantile(y_conf[,i],probs=probs)
    q_pred[i,] = quantile(y_pred[,i],probs=probs)
  }
  
  if (meth_type == 'gp') {
    y_m = extract(fit,"y_m")[[1]]
    y_gp = extract(fit,"y_gp")[[1]]
    q_m = q_gp = matrix(NA,ncol=2,nrow=length(data$x_pred))
    m_gp = m_m = c()
    for (i in 1:length(data$x_pred)) {
      m_m[i]   = mean(y_m[,i] )
      m_gp[i]  = mean(y_gp[,i])
      q_m[i,]  = quantile(y_m[,i] ,probs=probs)
      q_gp[i,] = quantile(y_gp[,i],probs=probs)
    }
  } 
  
  par(mfrow=c(1,1),mar=c(3,3,1.6,.2),mgp=c(2,.75,0), 
      pty='s',tcl=-0.5, cex=cexPlot)
  matplot(data$x_pred, data$x_pred, type='n',ylim=ylim, 
          xlab='T / K', 
          ylab=expression(paste('Residuals / ',mu,'Pa.s'))
  )
  grid(); box()
  # Define different symbols/colors for 10 sets
  cols= rep(1:5,2)
  pch = c(rep(16,5),rep(17,5))
  if(center) {
    # Center on best fit
    lp  = extract(fit,"lp__")[[1]]
    map = which.max(lp)
    # if(meth_type == 'hier' | meth_type == 'abc' | meth_type == 'gp') {
    #   eps = mean(extract(fit,"eps")[[1]])
    #   sig = mean(extract(fit,"sig")[[1]])
    # } else {
      eps = extract(fit,"eps")[[1]][map]
      sig = extract(fit,"sig")[[1]][map]
    # }
    c_new = phys_mod(data$x_pred,eps,sig,inadequacy)
    c_exp = phys_mod(data$x,eps,sig,inadequacy)
    res = resid[map,]
  
    # Centered CI
    # if (meth_type == 'gp') {
    #   polygon(cbind(data$x_pred,rev(data$x_pred)),
    #           cbind(q_pred[,1]-c_new,rev(q_pred[,2]-c_new)),
    #           col='gray90', border=NA)
    #   polygon(cbind(data$x_pred,rev(data$x_pred)),
    #           cbind(q_conf[,1]-c_new,rev(q_conf[,2]-c_new)),
    #           col='gray70', border=NA)
      # polygon(cbind(data$x_pred,rev(data$x_pred)),
      #         cbind(q_pred[,1]-m_gp-m_m,rev(q_pred[,2]-m_gp-m_m)),
      #         col='gray90', border=NA)
      # polygon(cbind(data$x_pred,rev(data$x_pred)),
      #         cbind(q_gp[,1]-m_gp,rev(q_gp[,2]-m_gp)),
      #         col='gray70', border=NA)
      # Decomposition of GP into contributions
      # polygon(cbind(data$x_pred,rev(data$x_pred)),
      #         cbind(q_m[,1]-c_new,rev(q_m[,2]-c_new)),
      #         col=col2tr('blue',60), border=NA)
      # polygon(cbind(data$x_pred,rev(data$x_pred)),
      #         cbind(q_gp[,1],rev(q_gp[,2])),
      #         col=col2tr('orange',120), border=NA)
    # } else {
      polygon(cbind(data$x_pred,rev(data$x_pred)),
              cbind(q_pred[,1]-c_new,rev(q_pred[,2]-c_new)),
              col='gray95', border=NA)
      polygon(cbind(data$x_pred,rev(data$x_pred)),
              cbind(q_conf[,1]-c_new,rev(q_conf[,2]-c_new)),
              col='gray70', border=NA)
    # }
    
    # Centered data
    points(data$x,res,pch=pch[data$tag], cex=0.75, 
           col=cols[data$tag])
    segments(data$x,res-2*data$uy,data$x,res+2*data$uy,
             lwd=2,col=cols[data$tag])
    
  } else {
    polygon(cbind(data$x_pred,rev(data$x_pred)),
            cbind(q_pred[,1],rev(q_pred[,2])),
            col='gray95', border=NA)
    polygon(cbind(data$x_pred,rev(data$x_pred)),
            cbind(q_conf[,1],rev(q_conf[,2])),
            col='gray70', border=NA)
    if (meth_type == 'gp') {
      polygon(cbind(data$x_pred,rev(data$x_pred)),
              cbind(q_m[,1],rev(q_m[,2])),
              col=col2tr('blue',60), border=NA)
      polygon(cbind(data$x_pred,rev(data$x_pred)),
              cbind(q_gp[,1],rev(q_gp[,2])),
              col=col2tr('orange',120), border=NA)
    }
    points(data$x,data$y,type='p',
           pch=pch[data$tag],col=cols[data$tag],cex=0.5)
    segments(data$x,data$y-2*data$uy,data$x,data$y+2*data$uy,
             lwd=2,col=cols[data$tag])
  }
  if(legend) {
    if (meth_type == 'gp') {
      legend('topleft',ncol=6,cex=0.75, 
             bty='o', bg=NULL, box.col=0,
             legend=c(parse(text=paste0('italic(D)^(',1:max(data$tag),')')),
                      'M(.)','K(.)'),
             col=c(cols[1:max(data$tag)],
                   col2tr('blue',60),
                   col2tr('orange',120)),
             pch=c(pch[1:max(data$tag)],15,15)
      )
    } else {
      legend('topleft',ncol=5,cex=0.75, 
             bty='o', bg=NULL, box.col=0,
             legend=parse(text=paste0('italic(D)^(',1:max(data$tag),')')),
             col=cols[1:max(data$tag)],
             pch=pch[1:max(data$tag)]
      )
    }
  }
}
plot_residuals <- function(ylim=NULL, sample=TRUE, 
                           legend=FALSE) {
  # Plot residuals
  par(mfrow=c(1,1),mar=c(3,3,1.6,.2),mgp=c(2,.75,0),
      pty='s',tcl=-0.5, cex=cexPlot)
  plot(data$x,data$x,type='n', ylim=ylim,
       xlab='T / K',
       ylab=expression(paste('Residuals / ',mu,'Pa.s'))
  )
  resid = extract(fit,"resid")[[1]]
  if(sample) {
    matplot(x,t(resid),type='p',pch=19,cex=1,col=col2tr('brown',5),
            add=TRUE)
  }
  abline(h=0)
  lp  = extract(fit,"lp__")[[1]]
  map = which.max(lp)
  cols= rep(1:5,2)
  pch = c(rep(16,5),rep(17,5))
  res = t(resid)[,map]
  
  points(data$x,res,pch=pch[data$tag], 
         cex=1, col=cols[data$tag])
  segments(data$x,res-2*data$uy,data$x,res+2*data$uy,
           lwd=2,col=cols[data$tag])
  
  if(legend)
    legend('topleft',ncol=5,cex=0.6,
           bty='o', bg=NULL, box.col=0,
           legend=parse(text=paste0('italic(D)^(',1:max(data$tag),')')),
           col=cols[1:max(data$tag)],
           pch=pch[1:max(data$tag)]
    )
  grid(); box()
}
plot_all  <- function(legend=FALSE) {
  # Generate all plots for current method
  dir_out = paste0('simulation/',case,'/figures')
  
  Title=model_tag
  
  # Type of method
  meth_type = 
    switch( substr(model_tag,1,2),
            Hi      = 'hier',
            AB      = 'abc',
            Ma      = 'abc',
            GP      = 'gp',
            'other'
    )
  
  if(substr(case,1,2)=='SD')
    true_point = data$true_point
  else
    true_point = NULL
  
  tiff(file=paste0(dir_out,'/sample_',model_tag,'.tiff'),
      width=reso, height=reso, compression="lzw")
  plot_hull(xlim=eps_lim, ylim=sig_lim, 
            meth_type=meth_type, true_point=true_point)
  title(main=Title)
  dev.off()
  
  tiff(file=paste0(dir_out,'/residuals_',model_tag,'.tiff'),
      width=reso, height=reso, compression="lzw")
  plot_residuals(ylim=res_lim,sample=FALSE,
                 legend=legend)
  title(main=Title)
  dev.off()
  
  tiff(file=paste0(dir_out,'/predict_',model_tag,'.tiff'),
      width=reso, height=reso, compression="lzw")
  plot_predict(center=FALSE, ylim=c(0,max(data$y)), 
               meth_type=meth_type, legend=legend)
  title(main=Title)
  dev.off()
  
  tiff(file=paste0(dir_out,'/predict_centered_',model_tag,'.tiff'),
      width=reso, height=reso, compression="lzw")
  plot_predict(center=TRUE, ylim=res_lim, 
               meth_type=meth_type, legend=legend)
  title(main=Title)
  dev.off()
  
}
plot_shifts <- function(model_tags, true_shifts=NULL, 
                        cexPlot=1,ylim=c(-1,1)) {
  # Print the optimized shifts for a series of methods
  
  dir_out = paste0('simulation/',case,'/figures')
  tiff(file=paste0(dir_out,'/compare_shifts.tiff'),
      width=3/2*reso,height=reso, compression="lzw")
  par(mfrow=c(1,1),mar=c(3,3,1.6,.2),mgp=c(2,.75,0),
      tcl=-0.5, cex=cexPlot)
  imod = 0
  for (model_tag in model_tags) {
    imod = imod + 1
    fit_tag = paste0("fit",model_tag)
    load(paste0('simulation/',case,'/',fit_tag,'.rda'))
    shifts = extract(fit,"shiftc")[[1]]
    nShift = ncol(shifts)
    qShift = matrix(NA,ncol=5,nrow=nShift)
    for (i in 1:nShift)
      qShift[i,] = quantile(shifts[,i] ,
                            probs=c(0.025,0.25,0.5,0.75,0.975))

    x = 1:nShift + 0.1*(imod-1) - 0.1*length(model_tags)/2
    if(imod == 1) {
      plot(x,qShift[,3],type='n',
           xlim=c(0.5,nShift+0.5),ylim=ylim,
           xlab='Series Index',            
           ylab=expression(paste('Shift parameters / ',mu,'Pa.s')),
           xaxt='n')
      axis(1,at=1:nShift)
      if(!is.null(true_shifts))
        segments(1:nShift-0.4,true_shifts,
                 1:nShift+0.4,true_shifts,
                 lwd=4,col=1,lty=2)
    }
    segments(x,qShift[,1],x,qShift[,5],lwd=2,col=imod)
    segments(x,qShift[,2],x,qShift[,4],lwd=5,col=imod)
    points(x,qShift[,3],pch=23,col=imod,cex=0.5)
  }
  # abline(v=1:nShift)
  grid()
  legend('topleft',legend=model_tags,lty=1,col=1:imod,ncol=3, lwd=3,
         bty='o',bg=NULL,box.col='white')
  box()
  dev.off()
}

plotVirial <- function(model_tag, cexPlot=1) {
  # Plot predicted Virial coefficients
  tiff(file=paste0('Figs/predictVirial_',model_tag,'.tiff'),
      width=reso,height=reso, compression="lzw")
  par(mfrow=c(1,1),mar=c(3,3,1.6,.2),mgp=c(2,.75,0),
      tcl=-0.5, cex=cexPlot)
  fit_tag = paste0("fit",model_tag)
  fit = get(fit_tag)
  mup = extract(fit,"mup")[[1]]
  nMC=500
  eps = mup[1:nMC,1]
  sig = mup[1:nMC,2]
  plot(virData[,1],virData[,2],pch=19,col=4,
       xlab='T / K',ylab=expression('Second Virial / L.mol-1'),
       main=model_tag)
  segments(virData[,1],virData[,2]-2*virData[,3],
           virData[,1],virData[,2]+2*virData[,3])
  for (i in 1:nMC) {
    points(virData[,1],virial2(virData[,1],eps[i],sig[i]),
           pch=19,col=orchid_tr)
  }

  legend('topleft',legend=c('Data','Prediction'),
         lty=0,pch=19,col=c('blue','orchid'),
         bty='o',bg=NULL,box.col='white')
  box(); grid()
  dev.off()
}