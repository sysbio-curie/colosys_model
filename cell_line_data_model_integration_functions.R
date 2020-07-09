###################################
# LOOP to calculate synergy values

# start_cell_line<-96; stop_cell_line<-length(unique(viab_df_tidy_all$cell_line))
# bb<-function_synergy_calculator(start_cell_line,stop_cell_line,viab_df_tidy_all,chek1_conc,mk2_conc, 
#                             n_interpol, fine_chek1_conc, fine_mk2_conc, mesh_fine)

function_synergy_calculator <- function(start_cell_line,stop_cell_line,viab_df_tidy_all,chek1_conc,mk2_conc, 
                                        n_interpol, fine_chek1_conc, fine_mk2_conc, mesh_fine,rev_val){

for (variable in seq(start_cell_line,stop_cell_line) ) {

  viab_df_single <- viab_df_tidy_all[viab_df_tidy_all$cell_line %in% unique(viab_df_tidy_all$cell_line)[variable],]
  if (rev_val=="rev"){
  y_chek1=rev(pull(viab_df_single[viab_df_single$mk2==0,],value)); y_mk2=rev(pull(viab_df_single[viab_df_single$chek1==0,],value))
  } else {
    y_chek1=pull(viab_df_single[viab_df_single$mk2==0,],value); y_mk2=pull(viab_df_single[viab_df_single$chek1==0,],value)
  }
  
  # plot(chek1_conc[2:length(x)],y_chek1[2:length(x)],lty=3,xlab="chek1",ylab="viability",log="x",ylim=c(0,max(c(y_mk2,y_chek1))))
  # lines(mk2_conc[2:length(mk2_conc)],y_mk2[2:length(y_mk2)],col="red")
  
  # fit sigmoidal to CHEK1 only: gamma + (1-gamma)/(1+(x/alpha)^b)
  x<-chek1_conc; y <- y_chek1
  a_init <- approx(x,y)$x[which.min(abs(approx(x,y)$y-(1+min(y))/2) )]; if (a_init==0){a_init=1e-3}
  init_vals <- list(a=a_init, b=2)
  m<-NULL; try(m <- nls(y ~ min(y) + (1-min(y))/(1+(x/a)^b), start=init_vals))
  if (!is.null(m)){
    params <- summary(m)$coefficients[,1]; alpha_chek1=params["a"]; beta_chek1=params["b"]; gamma_chek1=min(y)
  } else { 
    try(m <- gnls(y ~ min(y) + (1-min(y))/(1+(x/a)^b), start=init_vals))
    if (!is.null(m)){ alpha_chek1=summary(m)$coefficients["a"]; beta_chek1=summary(m)$coefficients["b"]; gamma_chek1=min(y) }
    }
  if (is.null(m)){ gamma_chek1<-beta_chek1<-alpha_chek1<-NA; 
  if (variable==1) {write.table("new simul",paste(folderpath,"errorlog.csv",sep=""),quote=F,row.names=F,col.names=F) } else {
  write.table(paste(variable,"chek1 cannot fit",sep = ","),paste(folderpath,"errorlog.csv",sep=""),quote=F,row.names=F,col.names=F,append=T) 
    }
  }
  
  # calculate expected x value for y=(1-gamma)*0.5, with inverse of the function
  # y_50=(1+min(y))/2; x_ex50_chek1 = alpha_chek1*( (y_50-1)/(gamma_chek1 - y_50) )^(1/beta_chek1)
  
  # fit sigmoidal to MK2 only: gamma + (1-gamma)/(1+(x/alpha)^b)
  x<-mk2_conc; y <- y_mk2; a_init <- approx(x,y)$x[which.min(abs(approx(x,y)$y-(1+min(y))/2))]; if (a_init==0){a_init=1e-3}
  init_vals <- list(a=a_init, b=2)
  m2 <-NULL; try(m2 <- nls(y ~ min(y) + (1-min(y))/(1+(x/a)^b), start=init_vals),silent = T)
  if (!is.null(m2)){
    params <- summary(m2)$coefficients[,1]; gamma_mk2=min(y); alpha_mk2=params["a"]; beta_mk2=params["b"]
  } else { 
    try(m2 <- gnls(y ~ min(y) + (1-min(y))/(1+(x/a)^b), start=init_vals))
    if (!is.null(m2)){ alpha_mk2=summary(m2)$coefficients["a"]; beta_mk2=summary(m2)$coefficients["b"]; gamma_mk2=min(y) }
  }
  if (is.null(m2)){ gamma_mk2<-beta_mk2<-alpha_mk2<-NA
  write.table(paste(variable,"mk2 cannot fit",sep = ","),paste(folderpath,"errorlog.csv",sep=""),quote=F,row.names=F,col.names=F,append=T)}
  
  # plot(x,y_mk2); lines(x,predict(m2),col="red",lty=2,lwd=3); cor(y,predict(m2))
  
  # all interpolated chek1 values
  # viab_df_single[viab_df_single$mk2!=0 & viab_df_single$chek1!=0,"value"]
  v_ab_vals=viab_df_single[viab_df_single$mk2!=0 & viab_df_single$chek1!=0,]
  if ( is.null(m)+is.null(m2) ==0 ){
    chek1_interpol_v_ab = alpha_chek1*( ( v_ab_vals$value -1)/(gamma_chek1 - v_ab_vals$value) )^(1/beta_chek1)
    # all interpolated mk2 values
    mk2_interpol_v_ab = alpha_mk2*( (v_ab_vals$value - 1)/(gamma_mk2 - v_ab_vals$value) )^(1/beta_mk2)
    # check interpolated points vs fitted sigmoidal
    # plot(chek1_interpol_v_ab, v_ab_vals$value, type='p', col = 'blue')
    # lines(chek1_interpol_v_ab[order(chek1_interpol_v_ab)],
    # min(y_chek1) + (1-min(y_chek1))/(1+(chek1_interpol_v_ab[order(chek1_interpol_v_ab)]/alpha_chek1)^beta_chek1), col="red")
    
    # Loewe synergy values (smaller than 1 means synergy)
    v_ab_vals$loewe_syn = v_ab_vals$mk2/mk2_interpol_v_ab + v_ab_vals$chek1/chek1_interpol_v_ab # <1 means synergy, >1 antag
  } else {
    v_ab_vals$loewe_syn <- NA
  }
  # where we cannot interpolate to this effect from single-agent effects, we can plug in 0 instead of NA
  # t1 = v_ab_vals$mk2/mk2_interpol_v_ab; t1[is.na(t1)]<-0; t2 = v_ab_vals$chek1/chek1_interpol_v_ab; t2[is.na(t2)]<-0
  # v_ab_vals$loewe_syn_corr <- t1+t2
  
  # matrix of viability values
  viab_matr_single <- acast(viab_df_single[,c(2,3,4)], mk2~chek1, value.var="value")
  # linear interpolation for loewe: syn = effect_cA/alpha + effect_cB/alpha
  effect_matrix=1-viab_matr_single; effect_matrix[effect_matrix<0]=0
  effect_matrix_first_row <- t(matrix(effect_matrix[1,2:ncol(effect_matrix)],
                                  length(effect_matrix[1,2:ncol(effect_matrix)]),length(2:nrow(effect_matrix)) ) )
  
  effect_matrix_first_col <- matrix(effect_matrix[2:nrow(effect_matrix),1],
                                    length(2:nrow(effect_matrix)),length(effect_matrix[1,2:ncol(effect_matrix)]))
  # 
  loewe_syn_linear_matrix <- effect_matrix_first_row/effect_matrix[2:nrow(effect_matrix),2:ncol(effect_matrix)] + 
    effect_matrix_first_col/effect_matrix[2:nrow(effect_matrix),2:ncol(effect_matrix)]
  loewe_syn_linear_df <- melt(loewe_syn_linear_matrix); colnames(loewe_syn_linear_df) <- c("mk2","chek1","loewe_syn_linear")
  if (sum(colnames(v_ab_vals) %in% "loewe_syn_linear")==0){
    v_ab_vals <- merge(v_ab_vals, loewe_syn_linear_df, by = c("mk2","chek1"),all.x=TRUE)
  }
  
  # Bliss-indep synergies (SMALLER THAN 0 means synergy)
  all_bliss_synergy_matrix <- viab_matr_single[2:nrow(viab_matr_single),2:ncol(viab_matr_single)] - 
    outer(viab_matr_single[2:nrow(viab_matr_single),1],viab_matr_single[1,2:ncol(viab_matr_single)])
  all_bliss_synergy_df <- melt(all_bliss_synergy_matrix); colnames(all_bliss_synergy_df) <- c("mk2","chek1","bliss_syn")
  ##
  # bliss synergy score at highest conc of inhibitor A can be calculated as area between the viab. curves with
  # both inhibitors and only B
  bliss_syn_area_maxmk2 <- trapz( as.numeric(colnames(viab_matr_single)),viab_matr_single[nrow(viab_matr_single),]/viab_matr_single[nrow(viab_matr_single),1] )-
    trapz(as.numeric(colnames(viab_matr_single)),viab_matr_single[1,] )
  bliss_syn_area_maxchek1 <- trapz(as.numeric(rownames(viab_matr_single)),viab_matr_single[,ncol(viab_matr_single)]/viab_matr_single[1,ncol(viab_matr_single)])-
    trapz(as.numeric(rownames(viab_matr_single)),viab_matr_single[,1] )
  all_bliss_synergy_df$bliss_area_maxmk2 <- bliss_syn_area_maxmk2; all_bliss_synergy_df$bliss_area_maxchek1 <- bliss_syn_area_maxchek1
  # merge with synergy values
  z <- merge(merge(viab_df_single, all_bliss_synergy_df, by = c("mk2","chek1"),all.x=TRUE), 
             v_ab_vals[,c("mk2","chek1","loewe_syn_linear","loewe_syn")], by=c("mk2","chek1"),all.x=TRUE)
  z$bliss_area_maxmk2[is.na(z$bliss_area_maxmk2)] <- unique(z$bliss_area_maxmk2)[!is.na(unique(z$bliss_area_maxmk2))]
  if (length(unique(z$bliss_area_maxchek1)[!is.na(unique(z$bliss_area_maxchek1))])>0){
    z$bliss_area_maxchek1[is.na(z$bliss_area_maxchek1)] <- unique(z$bliss_area_maxchek1)[!is.na(unique(z$bliss_area_maxchek1))]
  }
  
########################
# area-based synergy
  
 # interpolate matrix of drug effects
 interpol_effect_matr <- matrix( interp2(chek1_conc,mk2_conc,viab_matr_single,as.vector(mesh_fine$X),as.vector(mesh_fine$Y),method="linear"),
                                  dim(mesh_fine$X)[1], dim(mesh_fine$Y)[1] )
 colnames(interpol_effect_matr)<-fine_chek1_conc; rownames(interpol_effect_matr)<-fine_mk2_conc
  
  # EC50 curve
  chek1_conc_interpol_ec50 <- apply(abs(interpol_effect_matr-0.5),1,which.min)
  ec50_interpol_value <- unlist(lapply(1:nrow(interpol_effect_matr), function(n) {interpol_effect_matr[n,chek1_conc_interpol_ec50[n]]}))
  ec50_chek1_mk2_values <- data.frame(cbind(cbind(fine_chek1_conc[chek1_conc_interpol_ec50], fine_mk2_conc), ec50_interpol_value) )
  colnames(ec50_chek1_mk2_values)[1]<-"fine_chek1_conc"; rownames(ec50_chek1_mk2_values)<-c()
  
  # if no values close to 0.5, remove
  ec50_chek1_mk2_values <- ec50_chek1_mk2_values[!abs(ec50_chek1_mk2_values$ec50_interpol_value-0.5)>0.05,]
  # ec50_chek1_mk2_values <- ec50_chek1_mk2_values[!is.na(ec50_chek1_mk2_values[,1]),]
  ec50_chek1_mk2_values <- ec50_chek1_mk2_values[order(ec50_chek1_mk2_values$fine_mk2_conc),]
  # SMOOTH for outliers
  ec50_chek1_mk2_values$fine_mk2_conc <- smooth(ec50_chek1_mk2_values$fine_mk2_conc)
  ec50_chek1_mk2_values$fine_chek1_conc <- smooth(ec50_chek1_mk2_values$fine_chek1_conc)
  # if there are irregul.s in the curve
  x=ec50_chek1_mk2_values$fine_mk2_conc; y=ec50_chek1_mk2_values$fine_chek1_conc
  if (length(y)>1){
  if ( sum( y[2:length(y)]<0.5*y[1:(length(y)-1)] | y[2:length(y)]>2*y[1:(length(y)-1)] ) >1 ){
    min_abr=min(which(y[2:length(y)]>2*y[1:(length(y)-1)] | y[2:length(y)]<0.5*y[1:(length(y)-1)]))
    max_abr=max(which(y[2:length(y)]>2*y[1:(length(y)-1)] | y[2:length(y)]<0.5*y[1:(length(y)-1)]))
    shift_correction_min<-shift_correction_max<-5
    if (min_abr<=5) { shift_correction_min=min_abr-1 }; if (max_abr>=length(y)-5) { shift_correction_max=length(y)-max_abr }
    
    y[min_abr:max_abr] = mean(c(y[(min_abr-shift_correction_min):(min_abr-1)], y[(max_abr+1):(max_abr+5)] ))
    ec50_chek1_mk2_values$fine_mk2_conc=x; ec50_chek1_mk2_values$fine_chek1_conc=y
  }
    }

l_pred_mk2 <- function_predict_mk2(z,interpol_effect_matr,ec50_chek1_mk2_values,mk2_conc,fine_chek1_conc,alpha_mk2,gamma_mk2,beta_mk2)
l_pred_chek1 <- function_predict_chek1(z,interpol_effect_matr,ec50_chek1_mk2_values,chek1_conc,fine_mk2_conc,alpha_chek1,gamma_chek1,beta_chek1)

# function_predict_mk2 <- function(z,interpol_effect_matr,ec50_chek1_mk2_values,mk2_conc,fine_chek1_conc,alpha_mk2,gamma_mk2,beta_mk2)

# end of area-based synergy
########################

  print(variable)
  
if (variable==start_cell_line) {
  viab_df_tidy_all_synerg_pred_mk2 <- l_pred_mk2[[1]] # 
  viab_df_tidy_all_synerg_pred_chek1 <- l_pred_chek1[[1]] # viab_df_tidy_all_synerg_chek1 <- 
    
  } else {
    viab_df_tidy_all_synerg_pred_mk2 <- rbind(viab_df_tidy_all_synerg_pred_mk2, l_pred_mk2[[1]])
    viab_df_tidy_all_synerg_pred_chek1 <- rbind(viab_df_tidy_all_synerg_pred_chek1, l_pred_chek1[[1]])
  }
  
} # end of for loop

# single or multiple cell lines?
  if (start_cell_line==stop_cell_line) {
    list(l_pred_mk2, l_pred_chek1, ec50_chek1_mk2_values)
  } else {
    list(viab_df_tidy_all_synerg_pred_mk2, viab_df_tidy_all_synerg_pred_chek1)
  }

} # end of function

################################################
################################################

function_create_per_cell_line_table <- function(viab_df_tidy_all_synerg){
  x <- viab_df_tidy_all_synerg[[1]][viab_df_tidy_all_synerg[[1]]$mk2>0 & viab_df_tidy_all_synerg[[1]]$chek1>0,] %>% 
    group_by(cell_line) %>% summarise_all(funs(mean(.,na.rm=T)))
  x<-x[,!grepl("bliss_area",colnames(x))]; x <- x[,!colnames(x) %in% c("mk2","chek1","value")]
  # add "pred_mk2" tag
  colnames(x)[grepl("synergy_area",colnames(x))] <- paste(colnames(x)[grepl("synergy_area",colnames(x))],"_pred_mk2",sep = "")
  # add pred_chek1 area based synergy measures
  y <- unique(viab_df_tidy_all_synerg[[2]][,c(3,which(grepl("synergy_area",colnames(viab_df_tidy_all_synerg[[2]]))))])
  if ( !length(intersect(colnames(x), colnames(y)))>1 ){
    x <- cbind(x,y[,grepl("synergy_area",colnames(y))])
    colnames(x)[(ncol(x)-3):ncol(x)] <- paste(colnames(x)[(ncol(x)-3):ncol(x)], "_pred_chek1",sep=""); x[is.na(x)]<-NA; x[x==Inf]<-NA
    all_synergy_measures_per_cell_line<-x; rownames(all_synergy_measures_per_cell_line)<-c()
    all_synergy_measures_per_cell_line
  }
}

################################################
################################################
################################################

# calculate predicted EC50 curve at CHEK1=0 concentrations from the function viab=f(CHEK1)

function_predict_chek1 <- function(z,interpol_effect_matr,ec50_chek1_mk2_values,chek1_conc,fine_mk2_conc,
                                   alpha_chek1,gamma_chek1,beta_chek1){

# sort ec50_chek1_mk2_values by MK2 (x) values
ec50_chek1_mk2_values=ec50_chek1_mk2_values[order(ec50_chek1_mk2_values$fine_mk2_conc),]
  
# chek1 conc.s of EC50
mk2_only_interpol_concs_viab_above50 <- fine_mk2_conc[interpol_effect_matr[,1]>0.5]
# don't include values smaller than interpolated EC50 curve
mk2_only_interpol_concs_viab_above50<-mk2_only_interpol_concs_viab_above50[
    mk2_only_interpol_concs_viab_above50>=min(ec50_chek1_mk2_values$fine_mk2_conc,na.rm=T) & 
      mk2_only_interpol_concs_viab_above50<=max(ec50_chek1_mk2_values$fine_mk2_conc,na.rm=T)]

mk2_only_interpol_effects_viab_above50 <- as.numeric(interpol_effect_matr[fine_mk2_conc %in% mk2_only_interpol_concs_viab_above50,1])
# larger than 1 should not be allowed
if (length(which(mk2_only_interpol_effects_viab_above50>1))){
mk2_only_interpol_effects_viab_above50[1:max(which(mk2_only_interpol_effects_viab_above50>1))] <- 1
}
# mk2_only_interpol_effects_viab_above50>1 1:max(which(mk2_only_interpol_effects_viab_above50>1))
# what would be the required amount of MK2 to get to 0.5 viability? FIT with viab=f(mk2) sigmoidal
res_effect <- 1.5-mk2_only_interpol_effects_viab_above50
chek1_resid_sigmoid_fit = alpha_chek1*( (res_effect-1)/(gamma_chek1-res_effect) )^(1/beta_chek1)

# expected EC 50 curves
# take only part of the curve where predicted mk2 concentration is smaller than maximal in experiment
pred_ec50_mk2_restrictive <- mk2_only_interpol_concs_viab_above50[mk2_only_interpol_concs_viab_above50<=max(mk2_conc) & 
                                                                    !is.na(mk2_only_interpol_concs_viab_above50) & 
                                                                    chek1_resid_sigmoid_fit<max(chek1_conc) & 
                                                                    !is.na(chek1_resid_sigmoid_fit)]
pred_ec50_chek1_restrictive<-chek1_resid_sigmoid_fit[mk2_only_interpol_concs_viab_above50<=max(mk2_conc) & 
                                                       !is.na(mk2_only_interpol_concs_viab_above50) &
                                                       chek1_resid_sigmoid_fit<max(chek1_conc) & 
                                                       !is.na(chek1_resid_sigmoid_fit)]
# lines(pred_ec50_mk2_restrictive,pred_ec50_chek1_restrictive, col="red")

# take part of curve where inferred mk2 finite, but if larger than max(mk2), then make it equal to this
pred_ec50_chek1_finite_inner <- chek1_resid_sigmoid_fit[!is.na(chek1_resid_sigmoid_fit)]
# include area outside of square
pred_ec50_chek1_finite_outer<-pred_ec50_chek1_finite_inner
# only include area inside
pred_ec50_chek1_finite_inner[pred_ec50_chek1_finite_inner>max(chek1_conc)]<-max(chek1_conc)
pred_ec50_mk2_finite <- mk2_only_interpol_concs_viab_above50[!is.na(chek1_resid_sigmoid_fit)]
# include all area inside, also if mk2 predicted is NA  
pred_ec50_chek1_broad=chek1_resid_sigmoid_fit; 
pred_ec50_chek1_broad[is.na(pred_ec50_chek1_broad)|pred_ec50_chek1_broad>max(chek1_conc)]<-max(chek1_conc)

# have the interpolated curve with same resolution
if (sum(interpol_effect_matr<=0.5)>0 & sum(!sapply(ec50_chek1_mk2_values$fine_chek1_conc,is.na))>1 & 
    sum(!sapply(ec50_chek1_mk2_values$fine_mk2_conc,is.na))>1 ) {
  
  # area 1
  if ( length(pred_ec50_chek1_restrictive)>0 ) {
    
    truth_vector <- !(ec50_chek1_mk2_values$fine_mk2_conc>max(pred_ec50_mk2_restrictive)) &
      !(ec50_chek1_mk2_values$fine_mk2_conc<min(pred_ec50_mk2_restrictive))
    truth_vector[is.na(truth_vector)]<-FALSE
    
    z$synergy_area_restrict <- abs(trapz(pred_ec50_mk2_restrictive, pred_ec50_chek1_restrictive)) - 
      abs(trapz(ec50_chek1_mk2_values$fine_mk2_conc[truth_vector], ec50_chek1_mk2_values$fine_chek1_conc[truth_vector]) )
  } else {z$synergy_area_restrict <- NA}
  
  # if restrictive spans the whole curve, no need to calculate the others
  if (!exists("truth_vector")) {truth_vector=NA}
  if (sum(truth_vector,na.rm = T)==length(fine_mk2_conc)) {
    z$synergy_area_broad <- z$synergy_area_finite_inner <-z$synergy_area_finite_outer <- z$synergy_area_restrict
  } else {
  
  # area 2
  if (length(pred_ec50_mk2_finite)>0 )   {
    
    truth_vector <- !(ec50_chek1_mk2_values$fine_mk2_conc>max(pred_ec50_mk2_finite)) & 
      !(ec50_chek1_mk2_values$fine_mk2_conc<min(pred_ec50_mk2_finite))
    
    truth_vector[is.na(truth_vector)]<-FALSE
    
    z$synergy_area_finite_inner <- abs(trapz(pred_ec50_mk2_finite, pred_ec50_chek1_finite_inner)) - 
      abs(trapz(ec50_chek1_mk2_values$fine_mk2_conc[truth_vector], ec50_chek1_mk2_values$fine_chek1_conc[truth_vector]))
    # area 3
    z$synergy_area_finite_outer <- abs(trapz(pred_ec50_mk2_finite, pred_ec50_chek1_finite_outer)) - 
      abs(trapz(ec50_chek1_mk2_values$fine_mk2_conc[truth_vector], ec50_chek1_mk2_values$fine_chek1_conc[truth_vector] ))
    
  } else {z$synergy_area_finite_inner<-NA; z$synergy_area_finite_outer<-NA}
  
  # area 4
  if ( length(pred_ec50_chek1_broad)>0 & length(mk2_only_interpol_concs_viab_above50)>0 ) {
    truth_vector<-!(ec50_chek1_mk2_values$fine_mk2_conc>max(mk2_only_interpol_concs_viab_above50)) & 
      !(ec50_chek1_mk2_values$fine_mk2_conc<min(mk2_only_interpol_concs_viab_above50))
    truth_vector[is.na(truth_vector)]<-FALSE
    z$synergy_area_broad <- abs(trapz(mk2_only_interpol_concs_viab_above50, pred_ec50_chek1_broad)) - 
      abs(trapz(rev(ec50_chek1_mk2_values$fine_mk2_conc[truth_vector]), rev(ec50_chek1_mk2_values$fine_chek1_conc[truth_vector])))
  } else {z$synergy_area_broad <- NA}
  
    } # sum(truth_vector)==length(fine_mk2_conc)
  
} else {  z$synergy_area_restrict <- NA; z$synergy_area_finite_inner<-NA; z$synergy_area_finite_outer<-NA; z$synergy_area_broad <- NA }

list(z,pred_ec50_mk2_restrictive, pred_ec50_chek1_restrictive, pred_ec50_mk2_finite, pred_ec50_chek1_finite_inner, 
     pred_ec50_chek1_finite_outer, mk2_only_interpol_concs_viab_above50, pred_ec50_chek1_broad)

}

################################################
################################################

# calculate predicted EC50 curve at MK2=0 concentrations from the function viab=f(MK2)
# ec50_chek1_mk2_values

function_predict_mk2 <- function(z,interpol_effect_matr,ec50_chek1_mk2_values,mk2_conc,fine_chek1_conc,alpha_mk2,gamma_mk2,beta_mk2){

# sort
ec50_chek1_mk2_values=ec50_chek1_mk2_values[order(ec50_chek1_mk2_values$fine_chek1_conc),]
  
# chek1 conc.s of EC50
chek1_only_interpol_concs_viab_above50 <- fine_chek1_conc[interpol_effect_matr[1,]>0.5]
# don't include values smaller than interpolated EC50 curve
chek1_only_interpol_concs_viab_above50<-chek1_only_interpol_concs_viab_above50[
    chek1_only_interpol_concs_viab_above50>=min(ec50_chek1_mk2_values$fine_chek1_conc,na.rm=T)&
    chek1_only_interpol_concs_viab_above50<=max(ec50_chek1_mk2_values$fine_chek1_conc,na.rm=T)]

chek1_only_interpol_effects_viab_above50 <- as.numeric(interpol_effect_matr[1,fine_chek1_conc %in% chek1_only_interpol_concs_viab_above50])
# larger than 1 should not be allowed
if (length(which(chek1_only_interpol_effects_viab_above50>1))>1){
chek1_only_interpol_effects_viab_above50[1:max(which(chek1_only_interpol_effects_viab_above50>1))] <- 1
# chek1_only_interpol_effects_viab_above50>1
}
# what would be the required amount of MK2 to get to 0.5 viability? FIT with viab=f(mk2) sigmoidal
res_effect <- 1.5-chek1_only_interpol_effects_viab_above50
mk2_resid_sigmoid_fit = alpha_mk2*( (1 - res_effect)/(res_effect - gamma_mk2) )^(1/beta_mk2)

# expected EC 50 curves
# take only part of the curve where predicted mk2 concentration is smaller than maximal in experiment
pred_ec50_mk2_restrictive <- mk2_resid_sigmoid_fit[mk2_resid_sigmoid_fit<max(mk2_conc) & !is.na(mk2_resid_sigmoid_fit)]
pred_ec50_chek1_restrictive <- chek1_only_interpol_concs_viab_above50[mk2_resid_sigmoid_fit<max(mk2_conc) & !is.na(mk2_resid_sigmoid_fit)]
# lines(pred_ec50_chek1_restrictive,pred_ec50_mk2_restrictive, col="red")
upper_limit_pred_curve <- min(max(ec50_chek1_mk2_values$fine_mk2_conc),max(mk2_conc))

# take part of curve where inferred mk2 finite, but if larger than max(mk2), then make it equal to this
pred_ec50_mk2_finite_inner <- mk2_resid_sigmoid_fit[!is.na(mk2_resid_sigmoid_fit)]
# include area outside of square
pred_ec50_mk2_finite_outer<-pred_ec50_mk2_finite_inner
# only include area inside
pred_ec50_mk2_finite_inner[pred_ec50_mk2_finite_inner>upper_limit_pred_curve]<- upper_limit_pred_curve # max(mk2_conc)
pred_ec50_chek1_finite <- chek1_only_interpol_concs_viab_above50[!is.na(mk2_resid_sigmoid_fit)]
# include all area inside, also if mk2 predicted is NA  
pred_ec50_mk2_broad=mk2_resid_sigmoid_fit; 
pred_ec50_mk2_broad[is.na(pred_ec50_mk2_broad)|pred_ec50_mk2_broad>upper_limit_pred_curve]<-upper_limit_pred_curve

# have the interpolated curve with same resolution
if (sum(interpol_effect_matr<=0.5)>0 & sum(!sapply(ec50_chek1_mk2_values$fine_chek1_conc,is.na))>1 & 
    sum(!sapply(ec50_chek1_mk2_values$fine_mk2_conc,is.na))>1 ) {
  
  # area 1
  if ( length(pred_ec50_chek1_restrictive)>0 ) {
    
    truth_vector <- !(ec50_chek1_mk2_values$fine_chek1_conc>max(pred_ec50_chek1_restrictive)) &
      !(ec50_chek1_mk2_values$fine_chek1_conc<min(pred_ec50_chek1_restrictive))
    truth_vector[is.na(truth_vector)]<-FALSE
    
    z$synergy_area_restrict <- abs(trapz(pred_ec50_chek1_restrictive, pred_ec50_mk2_restrictive)) - 
      abs(trapz(ec50_chek1_mk2_values$fine_chek1_conc[truth_vector],ec50_chek1_mk2_values$fine_mk2_conc[truth_vector] ))
  } else {z$synergy_area_restrict <- NA}
  
  if (!exists("truth_vector")) {truth_vector=NA}
  # if restrictive spans the whole curve, no need to calculate the others
  if (sum(truth_vector,na.rm = T)==length(fine_chek1_conc)) {
    z$synergy_area_broad <- z$synergy_area_finite_inner <-z$synergy_area_finite_outer <- z$synergy_area_restrict
  } else {
  
  # area 2
  if (length(pred_ec50_chek1_finite)>0 ) {
    truth_vector <- !(ec50_chek1_mk2_values$fine_chek1_conc>max(pred_ec50_chek1_finite)) & 
      !(ec50_chek1_mk2_values$fine_chek1_conc<min(pred_ec50_chek1_finite))
    truth_vector[is.na(truth_vector)]<-FALSE
    z$synergy_area_finite_inner <- abs(trapz(pred_ec50_chek1_finite, pred_ec50_mk2_finite_inner)) - 
      abs(trapz(ec50_chek1_mk2_values$fine_chek1_conc[truth_vector],ec50_chek1_mk2_values$fine_mk2_conc[truth_vector]))
    
    # area 3
    z$synergy_area_finite_outer <- abs(trapz(pred_ec50_chek1_finite, pred_ec50_mk2_finite_outer)) - 
      abs(trapz(ec50_chek1_mk2_values$fine_chek1_conc[truth_vector],ec50_chek1_mk2_values$fine_mk2_conc[truth_vector] ))
    
  } else {z$synergy_area_finite_inner<-NA; z$synergy_area_finite_outer<-NA}
  
  # area 4
  if ( length(pred_ec50_mk2_broad)>0 & length(chek1_only_interpol_concs_viab_above50)>0 ) {
    
    truth_vector<-!(ec50_chek1_mk2_values$fine_chek1_conc>max(chek1_only_interpol_concs_viab_above50)) & 
      !(ec50_chek1_mk2_values$fine_chek1_conc<min(chek1_only_interpol_concs_viab_above50))
    truth_vector[is.na(truth_vector)]<-FALSE
    z$synergy_area_broad <- abs(trapz(chek1_only_interpol_concs_viab_above50, pred_ec50_mk2_broad)) - 
      abs(trapz(ec50_chek1_mk2_values$fine_chek1_conc[truth_vector],ec50_chek1_mk2_values$fine_mk2_conc[truth_vector]))
  } else {z$synergy_area_broad <- NA}
  
  }
}
  else { z$synergy_area_restrict <- NA; z$synergy_area_finite_inner<-NA; z$synergy_area_finite_outer<-NA; z$synergy_area_broad <- NA }

list(z, pred_ec50_chek1_restrictive, pred_ec50_mk2_restrictive, pred_ec50_chek1_finite, pred_ec50_mk2_finite_inner, 
     pred_ec50_mk2_finite_outer, chek1_only_interpol_concs_viab_above50, pred_ec50_mk2_broad)
}

######
# MEMO: global, local vars in functions
function_power <- function(x){
  a<-randn(1);
  z<-function_multipl(a,x)
  z }
function_multipl<-function(a,x){a*x}
# function_multipl<-function(x){a*x}
# would not work!
#####


################################################
################################################
# function plotter synergy areas

# list(l_pred_mk2, l_pred_chek1, ec50_chek1_mk2_values)
#
# predict CHEK1
# l_pred_chek1 <- list(z,pred_ec50_mk2_restrictive, pred_ec50_chek1_restrictive, pred_ec50_mk2_finite, pred_ec50_chek1_finite_inner, 
#      pred_ec50_chek1_finite_outer, mk2_only_interpol_concs_viab_above50, pred_ec50_chek1_broad)
# # predict mk2 
# l_pred_mk2 <- list(z,pred_ec50_chek1_restrictive, pred_ec50_mk2_restrictive, pred_ec50_chek1_finite, pred_ec50_mk2_finite_inner, 
# pred_ec50_mk2_finite_outer, chek1_only_interpol_concs_viab_above50, pred_ec50_mk2_broad)

function_synergy_area_plot_pred_mk2 <- function(height_width_vals, start_cell_line, xmax,ymax,
                                                l_pred, ec50_chek1_mk2_values, cex_lab_val,cex_main,cex_axis_val,margin_params){
  viab_df_tidy_indiv_cell_line_synerg<-l_pred[[1]]
  pred_ec50_chek1_restrictive=l_pred[[2]]; pred_ec50_mk2_restrictive<-l_pred[[3]]
  pred_ec50_chek1_finite<-l_pred[[4]]
  pred_ec50_mk2_finite_inner<-l_pred[[5]]; pred_ec50_mk2_finite_outer<-l_pred[[6]]; chek1_only_interpol_concs_viab_above50<-l_pred[[7]]
  pred_ec50_mk2_broad<-l_pred[[8]]
  # interpol curve same for pred mk2 and chek1 
  x<-ec50_chek1_mk2_values$fine_chek1_conc; y<-ec50_chek1_mk2_values$fine_mk2_conc
  
  syn_vals <- round(as.numeric(viab_df_tidy_indiv_cell_line_synerg[1,grepl("synergy_area",colnames(viab_df_tidy_indiv_cell_line_synerg))]),2)
  
  par(mfrow=c(2,2),mar=margin_params)
  # restrictive: until predicted curve is within area
  p %<a-% plot(x,y,xlim=c(0,xmax),xlab="chek1",ylab="mk2",ylim=c(0,ymax),cex=0.7,cex.lab=cex_lab_val,cex.axis=cex_axis_val)
  p; title(paste("restrictive,",syn_vals[1]),cex.main=cex_main,font.main=1)
  l1 %<a-% lines(seq(0,5,0.1),rep(10,51), lty="dashed"); l2 %<a-% lines(rep(5,101), seq(0,10,0.1), lty="dashed" ); l1;l2
  l3 %<a-% lines(pred_ec50_chek1_finite, pred_ec50_mk2_finite_outer,col="green",lty="dashed"); l3
  lines(pred_ec50_chek1_restrictive, pred_ec50_mk2_restrictive,xlim=c(0,5),ylim=c(0,10),col="red",lwd=3)
  # color area
  if (!is_empty(pred_ec50_chek1_restrictive) & is.finite(pred_ec50_chek1_restrictive[1]) ){
  vert_line_area<-seq(pred_ec50_mk2_restrictive[1], y[which.min(abs(x-min(pred_ec50_chek1_restrictive,na.rm=T)))], length.out = 100)
  polygon(c(x[x>=min(pred_ec50_chek1_restrictive)],rev(pred_ec50_chek1_restrictive),rep(min(pred_ec50_chek1_restrictive), length(vert_line_area))),
          c(y[x>=min(pred_ec50_chek1_restrictive)],rev(pred_ec50_mk2_restrictive), vert_line_area), col="red", border=NA)
  }
  mtext(unique(viab_df_tidy_all$cell_line)[start_cell_line],side=3,line=-3,outer=TRUE,cex=3,font.main=1)
  
  # finite inner
  p; title(paste("inner (finite values of pred. curve), ",syn_vals[2]),cex.main=cex_main,font.main=1)
  l1; l2; l3; lines(pred_ec50_chek1_finite, pred_ec50_mk2_finite_inner,col="blue",lwd=3)
  if (!is_empty(pred_ec50_mk2_finite_inner) ){ # & is.finite(pred_ec50_mk2_finite_inner[1]
  vert_line_area<-seq(pred_ec50_mk2_finite_inner[1], y[which.min(abs(x-min(pred_ec50_chek1_finite)))],length.out=100)
  polygon(c(x[x>=min(pred_ec50_chek1_finite)], rev(pred_ec50_chek1_finite), rep(min(pred_ec50_chek1_finite), length(vert_line_area)) ),
          c(y[x>=min(pred_ec50_chek1_finite)], rev(pred_ec50_mk2_finite_inner), vert_line_area), col="blue", border=NA )
  }
  
  # finite outer
  p; title(paste("outer (finite values of pred. curve),",syn_vals[3]),cex.main=cex_main,font.main=1) 
  l1; l2; lines(pred_ec50_chek1_finite, pred_ec50_mk2_finite_outer,col="green",lwd=3); l3
  if (!is_empty(pred_ec50_mk2_finite_outer) ){
    if (sum(is.infinite(pred_ec50_mk2_finite_outer))>0) { pred_ec50_mk2_finite_outer[is.infinite(pred_ec50_mk2_finite_outer)] <- ymax }
  vert_line_area<-seq(pred_ec50_mk2_finite_outer[1], y[which.min(abs(x-min(pred_ec50_chek1_finite)))],length.out=100)
  polygon(c(x[x>=min(pred_ec50_chek1_finite)], rev(pred_ec50_chek1_finite), rep(min(pred_ec50_chek1_finite), length(vert_line_area)) ),
          c(y[x>=min(pred_ec50_chek1_finite)], rev(pred_ec50_mk2_finite_outer), vert_line_area), col="green", border=NA )
  }
  
  # finite+infinite inner
  p; title(paste("broad (inner area for all values)",syn_vals[4]),cex.main=cex_main,font.main=1)
  l1;l2;lines(chek1_only_interpol_concs_viab_above50, pred_ec50_mk2_broad,col="orange",lwd=3);l3
  if (!is_empty(pred_ec50_mk2_broad) ) { # & is.finite(pred_ec50_mk2_broad[1])
  vert_line_area<-seq(pred_ec50_mk2_broad[1], y[which.min(abs(x-min(chek1_only_interpol_concs_viab_above50,na.rm = T)))],length.out=100)
  # vert_line_area<-seq(y[length(y)],pred_ec50_mk2_broad[length(pred_ec50_mk2_broad)], length.out=100)
  polygon(c( x[x>=min(chek1_only_interpol_concs_viab_above50)], rev(chek1_only_interpol_concs_viab_above50),
             rep(min(chek1_only_interpol_concs_viab_above50), length(vert_line_area))),
          c( y[x>=min(chek1_only_interpol_concs_viab_above50)],
             rev(pred_ec50_mk2_broad), vert_line_area ), col="orange", border=NA )
  }
}

##########################
# PLOT MK2-CHEK1 synergy area (CHEK1 is predicted from values at MK2 only)

# predict CHEK1
# l_pred_chek1 <- list(z,pred_ec50_mk2_restrictive, pred_ec50_chek1_restrictive, pred_ec50_mk2_finite, pred_ec50_chek1_finite_inner, 
#      pred_ec50_chek1_finite_outer, mk2_only_interpol_concs_viab_above50, pred_ec50_chek1_broad)

function_synergy_area_plot_pred_chek1 <- function(height_width_vals, start_cell_line, xmax,ymax,
                                                l_pred, ec50_chek1_mk2_values, cex_lab_val,cex_main,cex_axis_val,margin_params){
  
  viab_df_tidy_indiv_cell_line_synerg<-l_pred[[1]]
  pred_ec50_mk2_restrictive<-l_pred[[2]]; pred_ec50_chek1_restrictive<-l_pred[[3]]
  pred_ec50_mk2_finite<-l_pred[[4]]; pred_ec50_chek1_finite_inner<-l_pred[[5]]
  pred_ec50_chek1_finite_outer<-l_pred[[6]]; mk2_only_interpol_concs_viab_above50<-l_pred[[7]]
  pred_ec50_chek1_broad<-l_pred[[8]]
  # interpol curve same for pred mk2 and chek1 
  ec50_chek1_mk2_values <- ec50_chek1_mk2_values[order(ec50_chek1_mk2_values$fine_mk2_conc),]; rownames(ec50_chek1_mk2_values)<-c()
  x<-ec50_chek1_mk2_values$fine_mk2_conc; y<-ec50_chek1_mk2_values$fine_chek1_conc
  
  syn_vals <- round(as.numeric(viab_df_tidy_indiv_cell_line_synerg[1,grepl("synergy_area",colnames(viab_df_tidy_indiv_cell_line_synerg))]),2)
  
  par(mfrow=c(2,2),mar=margin_params)
  # restrictive: until predicted curve is within area
  p %<a-% plot(x,y,xlim=c(0,xmax),xlab="mk2",ylab="chek1",ylim=c(0,ymax),cex=0.7,cex.lab=cex_lab_val,cex.axis=cex_axis_val)
  p; title(paste("restrictive,",syn_vals[1]),cex.main=cex_main,font.main=1)
  l1 %<a-% lines(seq(0,10,0.1),rep(5,101), lty="dashed"); l2 %<a-% lines(rep(10,101),seq(0,5,0.05), lty="dashed" ); l1;l2
  if (!is_empty(pred_ec50_mk2_finite)){
  l3 %<a-% lines(pred_ec50_mk2_finite, pred_ec50_chek1_finite_outer,col="green",lty="dashed"); l3
  lines(pred_ec50_mk2_restrictive,pred_ec50_chek1_restrictive, xlim=c(0,5),ylim=c(0,10),col="red",lwd=3)
  } else {l3<-NULL}
  # color area
  if (!is_empty(pred_ec50_chek1_restrictive)){
    vert_line_area<-seq(y[length(y)],pred_ec50_chek1_restrictive[length(pred_ec50_chek1_restrictive)], length.out=100)
    polygon(c(x[x>=min(pred_ec50_mk2_restrictive,na.rm=T)],rep(max(pred_ec50_mk2_restrictive,na.rm=T),length(vert_line_area)),rev(pred_ec50_mk2_restrictive)),
            c(y[x>=min(pred_ec50_mk2_restrictive,na.rm=T)], vert_line_area,rev(pred_ec50_chek1_restrictive)), col="red", border=NA)
  }
  mtext(unique(viab_df_tidy_all$cell_line)[start_cell_line],side=3,line=-3,outer=TRUE,cex=3,font.main=1)
  
  # finite inner
  p; title(paste("inner (finite values of pred. curve), ",syn_vals[2]),cex.main=cex_main,font.main=1)
  l1; l2; l3;
  if (!is_empty(pred_ec50_mk2_finite)){
    lines(pred_ec50_mk2_finite,pred_ec50_chek1_finite_inner, col="blue",lwd=3)
    vert_line_area<-seq(y[length(y)],pred_ec50_chek1_finite_inner[length(pred_ec50_chek1_finite_inner)], length.out=100)
    polygon(c(x[x>=min(pred_ec50_mk2_finite,na.rm=T)],rep(max(pred_ec50_mk2_finite,na.rm=T),length(vert_line_area)),rev(pred_ec50_mk2_finite)),
            c(y[x>=min(pred_ec50_mk2_finite,na.rm=T)], vert_line_area,rev(pred_ec50_chek1_finite_inner)), col="blue", border=NA)
  }
  
  # finite outer
  p; title(paste("outer (finite values of pred. curve),",syn_vals[3]),cex.main=cex_main,font.main=1) 
  l1; l2; l3
  if (!is_empty(pred_ec50_mk2_finite)){
    lines(pred_ec50_mk2_finite, pred_ec50_chek1_finite_outer,col="green",lwd=3); 
    if (sum(is.infinite(pred_ec50_chek1_finite_outer))>0) { pred_ec50_chek1_finite_outer[is.infinite(pred_ec50_chek1_finite_outer)] <- ymax }
    
    vert_line_area<-seq(y[length(y)],pred_ec50_chek1_finite_outer[length(pred_ec50_chek1_finite_outer)], length.out=100)
    polygon(c(x[x>=min(pred_ec50_mk2_finite,na.rm=T)],rep(max(pred_ec50_mk2_finite,na.rm=T),length(vert_line_area)),rev(pred_ec50_mk2_finite)),
            c(y[x>=min(pred_ec50_mk2_finite,na.rm=T)], vert_line_area,rev(pred_ec50_chek1_finite_outer)), col="green", border=NA)
  }
  
  # finite+infinite inner (broad)
  p; title(paste("broad (inner area for all values)",syn_vals[4]),cex.main=cex_main,font.main=1)
  l1;l2;l3
  if (!is_empty(mk2_only_interpol_concs_viab_above50)){
    lines(mk2_only_interpol_concs_viab_above50, pred_ec50_chek1_broad,col="orange",lwd=3);
    vert_line_area<-seq(y[length(y)],pred_ec50_chek1_broad[length(pred_ec50_chek1_broad)], length.out=100)
    polygon(c(x[x>=min(mk2_only_interpol_concs_viab_above50,na.rm=T)],rep(max(mk2_only_interpol_concs_viab_above50,na.rm=T),length(vert_line_area)),
              rev(mk2_only_interpol_concs_viab_above50)),
            c(y[x>=min(mk2_only_interpol_concs_viab_above50,na.rm=T)], vert_line_area,rev(pred_ec50_chek1_broad)), col="orange", border=NA)
  }
}

############
# scatter plots comparing different area-based synergy measures

function_synergy_area_compare_plot <- function(df_synergies,list_cols,
                                               cex_val,cex_lab_val,cex_axis_val,cex_main_val,mtext_string,mtext_val){
  for (counter in 1:length(list_cols) ) {
    xlab_str=str_replace(colnames(df_synergies)[list_cols[[counter]][1]+1],"synergy_area_","")
    ylab_str=str_replace(colnames(df_synergies)[list_cols[[counter]][2]+1],"synergy_area_","")
    x=df_synergies[,list_cols[[counter]][1]+1]; y=df_synergies[,list_cols[[counter]][2]+1]
    cor_val=round(cor(x,y,use="pairwise.complete.obs"),2)
    title_str=paste(xlab_str, "vs", ylab_str, "cor=", cor_val)
    xlim_vals=c(quantile(x,0.05,na.rm=T),quantile(x,0.95,na.rm=T)); ylim_vals=c(quantile(y,0.05,na.rm=T),quantile(y,0.95,na.rm=T))
    plot(x,y,xlab=xlab_str, ylab=ylab_str, cex=cex_val, cex.lab=cex_lab_val,cex.axis=cex_axis_val,xlim=xlim_vals,ylim=ylim_vals)
    title(title_str,cex.main=cex_main_val,font.main=1)
    mtext(mtext_string,side=3,line=-3,outer=TRUE,cex=mtext_val,font.main=1)
  }
}

##########################
##########################

function_plot_name <- function(list_subplots,filename) {
  if (max(unlist(list_subplots))>7) {
    filename<-str_replace(filename,".eps", "_predchek1.eps")
  } else {filename<-str_replace(filename,".eps", "_predmk2.eps")}
  filename
}

function_comparison_all_synergy_measures <- function(list_subplots,syn_measures,mar_vals,subplot_vals,cex_lab_val,cex_axis_val){
  list_comp_syn <- list_subplots
  par(mfrow=subplot_vals,mar=mar_vals); all_synergy_measures_per_cell_line <- syn_measures
  for (counter in 1:length(list_comp_syn)) {
    plot(all_synergy_measures_per_cell_line[,list_comp_syn[[counter]][1]+1],
         all_synergy_measures_per_cell_line[,list_comp_syn[[counter]][2]+1],
         xlim=as.numeric(lims[list_comp_syn[[counter]][1],]), ylim=as.numeric(lims[list_comp_syn[[counter]][2],]),
         xlab=str_replace(rownames(lims)[list_comp_syn[[counter]][1]],"synergy_area_",""),
         ylab=str_replace(rownames(lims)[list_comp_syn[[counter]][2]],"synergy_area_",""), 
         cex=3, cex.lab=cex_lab_val,cex.axis=cex_axis_val)
    title_str=paste("corr=",round(cor(all_synergy_measures_per_cell_line[,list_comp_syn[[counter]][1]+1], 
                      all_synergy_measures_per_cell_line[,list_comp_syn[[counter]][2]+1],use="pairwise.complete.obs"),2),sep = "")
    title(title_str,cex.main=4,font.main=1)
  }
  # filename
}

###################################

function_single_cell_line_plotter <- function(file_path, vector){

for (counter in vector) {
start_cell_line=stop_cell_line=counter
# calculate curves and areas for one cell line
list_pars_indiv_cell_line <- function_synergy_calculator(start_cell_line,stop_cell_line,viab_df_tidy_all,chek1_conc,mk2_conc,
                                                         n_interpol,fine_chek1_conc,fine_mk2_conc,mesh_fine)
# VALUES!!
l_pred_mk2<-list_pars_indiv_cell_line[[1]]; l_pred_chek1<-list_pars_indiv_cell_line[[2]]
ec50_chek1_mk2_values<-list_pars_indiv_cell_line[[3]]; ec50_chek1_mk2_values=list_pars_indiv_cell_line[[3]]

# PLOT FOR ONE CELL LINE with 4 SYN MEASURES
height_width_vals<-c(12,18); cex_lab_val<-2.5;cex_main<-3;cex_axis_val=1.5;margin_params<-c(5,5,5,2)
# MK2 predicted
# SAVE
postscript(paste(file_path,unique(viab_df_tidy_all$cell_line)[start_cell_line],"_pred_mk2_synergy_areas.eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
# PLOTTING FUNCTION
ymax<-15; cex_lab_val<-2.5; cex_main<-3; cex_axis_val=1.5; margin_params<-c(5,5,5,2)
xmax<-max(l_pred_mk2[[length(l_pred_mk2)-1]])*1.1; if (is.infinite(xmax)) { xmax=5.5}
function_synergy_area_plot_pred_mk2(height_width_vals, start_cell_line,xmax,ymax,
                                    l_pred_mk2,ec50_chek1_mk2_values,cex_lab_val,cex_main,cex_axis_val,margin_params)
dev.off()

# CHEK1 is predicted
postscript(paste(file_path,unique(viab_df_tidy_all$cell_line)[start_cell_line],"_pred_chek1_synergy_areas.eps",sep=""),
           height=height_width_vals[1],width=height_width_vals[2],onefile=FALSE,horizontal=FALSE,paper='special')
# PLOTTING FUNCTION
xmax <- max(l_pred_chek1[[length(l_pred_chek1)-1]])*1.1; ymax<-5.5;
function_synergy_area_plot_pred_chek1(height_width_vals,start_cell_line,xmax,ymax,
                                      l_pred_chek1,ec50_chek1_mk2_values,cex_lab_val,cex_main,cex_axis_val,margin_params)
dev.off()
} # for loop
}