##########################################################
# INPUTS
#  ProxyDat := An N vector of the proxy variable
#  GoodDat  := An N vector of the actual variable (observations that were not labelled should be set to NA)
#  NOTE: ProxyDat and GoodDat should have the same column names 

#  InputWeights := These are the W_i values in the corresponding paper

#  alpha := Statistical significance level. The default is 0.1
#  q := the quantile to estimate


#  TuningScheme := character scheme to use
#  "OptDiagonal" Estimates optimal tuning matrix among diagonal tuning matrices
#  "None" sets Omega to the identity
#  "Zero" sets tuning matrix to zero (used for debugging purposes or to study classical estimator)

##########################################################
library(ggdist)
library(topolow)
weighted_PTD_quantile <-function(ProxyDat,GoodDat,InputWeights,q,alpha=0.1, TuningScheme="OptDiagonal"){
  

  
  #Helpful variable definitions
  N <- length(ProxyDat)

  IdxComplete <- complete.cases(GoodDat)
  #Double check that this corresponds to set where W is not zero
  WisNotZero <- InputWeights !=0
  if(mean(IdxComplete==WisNotZero)<1){print("Warning: complete cases and samples with W \neq 0 don't match")}
  
  #Extract Good and Proxy Data on the Complete cases 
  GoodDatCalib <- GoodDat[IdxComplete]
  ProxyDatCalib <- ProxyDat[IdxComplete]
  WCalib <- InputWeights[IdxComplete]
  
  ##################### Calculating point initial point estimates ####################
  
    #betaHat
    thetaHatCalib <- weighted_quantile(x = GoodDatCalib,probs = q,weights = WCalib)

    #GammaHat
    gammaHatCalib <- weighted_quantile(x = ProxyDatCalib,probs = q,weights = WCalib)

    #gamma hat on all samples (naive approach)
    gammaHatAll <- weighted_quantile(x = ProxyDat,probs = q,weights = rep(1,N))

    ##################### Calculating asymptotic covariance matrix components ####################
    # If you study the derivation of asymptotic linear expansions for quantiles you will see that
    # \dot{l}_{\theta}(Y)=I{Y \leq \theta} - q
    # H_{\theta} = f(\theta) #f is the pdf of Y
    
    
    
    dotLossGoodCalib <- I(GoodDatCalib <= thetaHatCalib)-q
    dotLossProxyCalib <- I(ProxyDatCalib <= gammaHatAll)-q
    dotLossProxyAll <- I(ProxyDat <= gammaHatAll)-q
    

    #Sigma Terms
    Sig11 <- sum(WCalib* dotLossGoodCalib* WCalib*dotLossGoodCalib)/N
    Sig12 <- sum(WCalib* dotLossGoodCalib* WCalib*dotLossProxyCalib)/N
    Sig22 <- sum(WCalib* dotLossProxyCalib* WCalib*dotLossProxyCalib)/N

    Sig13 <-  sum(WCalib* dotLossGoodCalib* dotLossProxyCalib)/N
    Sig33 <- sum(dotLossProxyAll*dotLossProxyAll)/N
    
    #Hessians (NEED TO ESTIMATE DENSITY AT THE QUANTILES)
    HessTheta <- weighted_kde(x = GoodDatCalib,weights = WCalib,n = 1,from = thetaHatCalib,to = thetaHatCalib)$y
    HessGamma <- weighted_kde(x = ProxyDat,weights = rep(1,N),n = 1,from = gammaHatAll,to = gammaHatAll)$y
    #HessGamma <- weighted_kde(x = ProxyDatCalib,weights =WCalib,n = 1,from = gammaHatCalib,to = gammaHatCalib)$y
    
    

    ##################### Calculating desired tuning matrices ####################
      if(TuningScheme=="None"){OmegaTuned <- 1}
      if(TuningScheme=="Zero"){OmegaTuned <- 0}
      if(TuningScheme %in% c("OptDiagonal","Diagonal") ){
          OmegaTuned <- (1/HessTheta) * (Sig12-Sig13) * (1/(Sig22-Sig33)) * HessGamma
      }
     
      
      
    ##################### Calculating estimator and CIs ####################
    
      #PTD point estimate
      thetaPTDest <- thetaHatCalib+ OmegaTuned * (gammaHatAll-gammaHatCalib)
      
      
      #Calculating asymptotic variance using CLT
      SigMatComb <- rbind(cbind(Sig11,Sig12,Sig13),
                          cbind(t(Sig12),Sig22,Sig33),
                          cbind(t(Sig13),Sig33,Sig33))
      OmegaInvHessGamma <- OmegaTuned/HessGamma
      AMat <- cbind(-1/HessTheta, OmegaInvHessGamma,-OmegaInvHessGamma)
      
      AsympCovPTD <- AMat %*% SigMatComb %*% t(AMat)
      
      #Calculating CIs
      SEPTDAsympEst <- sqrt(diag(AsympCovPTD/N))
      CI_lbs <- thetaPTDest-qnorm(1-alpha/2)*SEPTDAsympEst
      CI_ubs <- thetaPTDest+qnorm(1-alpha/2)*SEPTDAsympEst
      
      CIsOut <- c(thetaPTDest,CI_lbs,CI_ubs)
      names(CIsOut) <- c("Estimate","CI_lb","CI_ub") 
      
      # For complete samples calculating estimates of H_{\theta}^{-1} \dot{l}_{\theta}(Y)-H_{\gamma}^{-1} \dot{l}_{\gamma}(\tilde{Y})
      DiffInfluenceCalib <- dotLossGoodCalib/HessTheta -dotLossProxyCalib/HessGamma
      datForLabelingRule <- list(DiffInfluenceCalib=DiffInfluenceCalib,IdxComplete=IdxComplete)
      
      #Calculating CIs for classical baseline
      
      AsympCovClassical <-  Sig11/(HessTheta^2) 
      SEPTDAsympEstClassical <- sqrt(AsympCovClassical/N)
      CI_lbsClass <- thetaHatCalib-qnorm(1-alpha/2)*SEPTDAsympEstClassical
      CI_ubsClass <- thetaHatCalib+qnorm(1-alpha/2)*SEPTDAsympEstClassical
      
      CIsClassical <- c(thetaHatCalib,CI_lbsClass,CI_ubsClass)
      names(CIsClassical) <- c("Estimate","CI_lb","CI_ub") 
      

  return(list(OmegaTuned=OmegaTuned,CIsPTD=CIsOut,CIsClassical=CIsClassical,datForLabelingRule=datForLabelingRule))
}