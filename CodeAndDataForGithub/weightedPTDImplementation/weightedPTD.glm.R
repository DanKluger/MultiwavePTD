##########################################################
# INPUTS
#  ProxyDat := An N x (p+1) data frame where all variables that are not widely available are replaced by their widely available proxies
#  GoodDat  := An N x (p+1) data frame with the actual data. Rows corresponding to data points that were not labelled should be set to NA
#  NOTE: ProxyDat and GoodDat should have the same column names 

#  InputWeights := These are the W_i values in the corresponding paper

#  OutcomeVarName   := The name of the outcome variable in your data frame. This function will estimate the regression of that variable on the remaining variables.

#  RegType   := Use "linear" for linear regression or "logistic" for logistic regression

#  alpha := Statistical significance level. The default is 0.1

#  TuningScheme := character scheme to use
#  "OptDiagonal" Estimates optimal tuning matrix among diagonal tuning matrices
#  "None" sets Omega to the identity
#  "Zero" sets tuning matrix to zero (used for debugging purposes or to study classical estimator)

##########################################################

weighted_PTD_glm <-function(ProxyDat,GoodDat,InputWeights,OutcomeVarName="Y",RegType="linear",alpha=0.1, TuningScheme="OptDiagonal"){
  
  #Setting up glm implementation within function
  RegFormulaLocal <- as.formula(paste0(OutcomeVarName,"~ ."))
  if(RegType %in% c("linear","Linear","Gaussian","gaussian","Normal","normal","OLS","ols")){
    psidot <- function(s){return(s)}
    psidotdot <- function(s){return(rep(1,length(s)))}
    familyLocal <- gaussian(link = "identity")
  } else if(RegType %in% c("logistic","Logistic","logit")){
    psidot <- function(s){return(1/(1+exp(-s)))}
    psidotdot <- function(s){return(exp(s)/((1+exp(s))^2))}
    familyLocal <- binomial(link = "logit")
  }
  
  #Helpful variable definitions
  N <- nrow(ProxyDat)
  d <- ncol(GoodDat)
  
  IdxComplete <- complete.cases(GoodDat)
  #Double check that this corresponds to set where W is not zero
  WisNotZero <- InputWeights !=0
  if(mean(IdxComplete==WisNotZero)<1){print("Warning: complete cases and samples with W \neq 0 don't match")}
  
  #Extract Good and Proxy Data on the Complete cases 
  GoodDatCalib <- GoodDat[IdxComplete,]
  ProxyDatCalib <- ProxyDat[IdxComplete,]
  WCalib <- InputWeights[IdxComplete]
  
  ##################### Calculating point initial point estimates ####################
  
    #betaHat
    thetaHatGlmFit <- glm(formula = RegFormulaLocal,family = familyLocal,data = GoodDatCalib,weights = WCalib)
    thetaHatCalib <- thetaHatGlmFit$coefficients
    CoefNamesThetaHat <- names(thetaHatCalib)

    #GammaHat
    gammaHatCalib <- glm(formula = RegFormulaLocal,family = familyLocal,data = ProxyDatCalib,weights = WCalib)$coefficients

    #gamma hat on all samples (naive approach)
    gammaHatAll <- glm(formula = RegFormulaLocal,family = familyLocal,data = ProxyDat)$coefficients

    ##################### Calculating asymptotic covariance matrix components ####################
    
    GoodDatCovsCalib <- cbind(rep(1,nrow(GoodDatCalib)),as.matrix(GoodDatCalib)[,which(names(GoodDatCalib) != OutcomeVarName )] )
    ProxyDatCovsCalib <- cbind(rep(1,nrow(ProxyDatCalib)),as.matrix(ProxyDatCalib)[,which(names(ProxyDatCalib) != OutcomeVarName )] )
    ProxyDatCovs <- cbind(rep(1,nrow(ProxyDat)),as.matrix(ProxyDat)[,which(names(ProxyDat) != OutcomeVarName )] )
    
    
    dotLossGoodCalib <- matrix(rep((psidot( GoodDatCovsCalib %*% thetaHatCalib)-GoodDatCalib[[OutcomeVarName]])
                                  ,d),ncol = d,byrow = F)*GoodDatCovsCalib
    dotLossProxyCalib <- matrix(rep((psidot( ProxyDatCovsCalib %*% gammaHatAll)-ProxyDatCalib[[OutcomeVarName]])
                                   ,d),ncol = d,byrow = F)*ProxyDatCovsCalib
    dotLossProxyAll <- matrix(rep((psidot( ProxyDatCovs %*% gammaHatAll)-ProxyDat[[OutcomeVarName]])
                                    ,d),ncol = d,byrow = F)*ProxyDatCovs
    
    WeightsCalibMat <- matrix(rep(WCalib,d),ncol=d,byrow=F)
    
    #Sigma Terms
    Sig11 <- (t(WeightsCalibMat*dotLossGoodCalib) %*% (WeightsCalibMat*dotLossGoodCalib))/N
    Sig12 <- (t(WeightsCalibMat*dotLossGoodCalib) %*% (WeightsCalibMat*dotLossProxyCalib))/N
    Sig22 <- (t(WeightsCalibMat*dotLossProxyCalib) %*% (WeightsCalibMat*dotLossProxyCalib))/N
    
    Sig13 <-  (t(WeightsCalibMat*dotLossGoodCalib) %*% (dotLossProxyCalib))/N
    Sig33 <- (t(dotLossProxyAll) %*% (dotLossProxyAll))/N
    
    #Hessians
    #H_{\theta}=E[W psidotdot(\theta^\tran X) XX^\tran]
    #H_{\gamma}= E[psidotdot(\gamma^\tran \tilde{X}) \tilde{X} \tilde{X}^\tran]
    
    weightHessTheta <- WCalib * psidotdot(GoodDatCovsCalib %*% thetaHatCalib)
    weightHessGamma <- psidotdot(ProxyDatCovs %*% gammaHatAll)
    
    weightHessThetaMat <- matrix(rep(weightHessTheta,d),ncol=d,byrow=F)
    weightHessGammaMat <- matrix(rep(weightHessGamma,d),ncol=d,byrow=F)
    
    HessTheta <- (t(weightHessThetaMat * GoodDatCovsCalib) %*% GoodDatCovsCalib)/N
    HessGamma <- (t(weightHessGammaMat * ProxyDatCovs) %*% ProxyDatCovs)/N
    
    

    ##################### Calculating desired tuning matrices ####################
      if(TuningScheme=="None"){OmegaTuned <- diag(d)}
      if(TuningScheme=="Zero"){OmegaTuned <- matrix(0,nrow= d,ncol= d)}
      if(TuningScheme %in% c("OptDiagonal","Diagonal") ){
        numVec <- diag(solve(HessTheta) %*% (Sig12-Sig13) %*% solve(HessGamma))
        denVec <- diag(solve(HessGamma) %*% (Sig22-Sig33) %*% solve(HessGamma))
        OmegaTuned <- diag(numVec/denVec)
      }
    
      
      
    ##################### Calculating estimator and CIs ####################
    
      #PTD point estimate
      thetaPTDest <- as.vector(thetaHatCalib+ OmegaTuned %*% (gammaHatAll-gammaHatCalib))
      
      
      #Calculating asymptotic variance using CLT
      SigMatComb <- rbind(cbind(Sig11,Sig12,Sig13),
                          cbind(t(Sig12),Sig22,Sig33),
                          cbind(t(Sig13),Sig33,Sig33))
      OmegaInvHessGamma <- OmegaTuned %*% solve(HessGamma)
      InvHessTheta <- solve(HessTheta)
      AMat <- cbind(-InvHessTheta, OmegaInvHessGamma,-OmegaInvHessGamma)
      
      AsympCovPTD <- AMat %*% SigMatComb %*% t(AMat)
      
      #Calculating CIs
      SEPTDAsympEst <- sqrt(diag(AsympCovPTD/N))
      CI_lbs <- matrix(thetaPTDest-qnorm(1-alpha/2)*SEPTDAsympEst,ncol = 1)
      CI_ubs <- matrix(thetaPTDest+qnorm(1-alpha/2)*SEPTDAsympEst,ncol = 1)
      
      CIMatOut <- cbind(matrix(thetaPTDest,ncol=1),CI_lbs,CI_ubs)
      rownames(CIMatOut) <- CoefNamesThetaHat
      colnames(CIMatOut) <- c("Estimate","CI_lb","CI_ub") 
      
      
      # For complete samples calculating estimates of H_{\theta}^{-1} \dot{l}_{\theta}(Y)-H_{\gamma}^{-1} \dot{l}_{\gamma}(\tilde{Y})
      DiffInfluenceCalib <- t(solve(HessTheta) %*% t(dotLossGoodCalib) -solve(HessGamma) %*% t(dotLossProxyCalib) )
      datForLabelingRule <- list(DiffInfluenceCalib=DiffInfluenceCalib,IdxComplete=IdxComplete)

      
      #Calculating CIs for classical baseline
      
      AsympCovClassical <- InvHessTheta %*% Sig11 %*% t(InvHessTheta)
      SEPTDAsympEstClassical <- sqrt(diag(AsympCovClassical/N))
      CI_lbsClass <- matrix(thetaHatCalib-qnorm(1-alpha/2)*SEPTDAsympEstClassical,ncol = 1)
      CI_ubsClass <- matrix(thetaHatCalib+qnorm(1-alpha/2)*SEPTDAsympEstClassical,ncol = 1)
      
      CIMatClassical <- cbind(matrix(thetaHatCalib,ncol=1),CI_lbsClass,CI_ubsClass)
      rownames(CIMatClassical) <- CoefNamesThetaHat
      colnames(CIMatClassical) <- c("Estimate","CI_lb","CI_ub") 
      
      
      
      #ClassicalFitTable <- summary(thetaHatGlmFit)$coefficients
      #CIMatClassical <- cbind(ClassicalFitTable[,'Estimate'],
       #                       ClassicalFitTable[,'Estimate']-qnorm(1-alpha/2)*ClassicalFitTable[,'Std. Error'],
      #                        ClassicalFitTable[,'Estimate']+qnorm(1-alpha/2)*ClassicalFitTable[,'Std. Error'])
      #colnames(CIMatClassical) <- c("Estimate","CI_lb","CI_ub") 
      
  return(list(OmegaTuned=OmegaTuned,CIsPTD=CIMatOut,CIsClassical=CIMatClassical,datForLabelingRule=datForLabelingRule))
}