#!/usr/bin/env Rscript


# change case of d or nd
perftable = function (list_predict, list_real, verbose = 0){
  
  nb_value = length (list_real)
  i = 1
  tp = 0.0
  fp = 0.0
  tn = 0.0
  fn = 0.0
  while(i <= nb_value ){
    if (list_predict[i]==1){
      if (list_predict[i] == list_real[i]){
        tp = tp + 1
      }else {
        fp = fp + 1
      }
    }else{
      if (list_predict[i] == list_real[i]){
        tn = tn + 1
      }else {
        fn = fn + 1
      }
    }
    i = i + 1
  }
  if(verbose == 1){
    print (paste ("TP : ", tp, sep = ""))
    print (paste ("TN : ", tn, sep = ""))
    print (paste ("FP : ", fp, sep = ""))
    print (paste ("FN : ", fn, sep = ""))
  }
  tableval = c(tp,tn,fp,fn)
  return (tableval)
}



accuracy = function (tp, tn, fp, fn){
  return ((tp + tn)/(tp + fp + tn +fn))
}

precision = function (tp, fp){
  return (tp/(tp + fp))
}

recall = function (tp, fn){
  return (tp/(tp + fn))
}

specificity = function (tn, fp){
  sp = tn/(tn + fp)
  if(is.na(sp)){
    sp = 0
  }
  return (sp)
}

sensibility = function (tp, fn){
  se = tp/(tp + fn)
  if(is.na(se)){
    se = 0
  }
  return (se)
}


BCR = function (tp, tn, fp, fn){
  return (0.5*(tp/(tp+fn) + tn/(tn+fp)))
  
}

MCC = function (tp, tn, fp, fn){
  numerator = tp*tn-fp*fn
  denumerator = (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)
  mcc = numerator / sqrt(denumerator)
  if(is.na(mcc)){
    mcc = 0
  }
  return (mcc)
}


qualityPredictList = function (test_vector, real_vector){
  v_predict = perftable (test_vector, real_vector, 1)
  print (paste ("accuracy : ", accuracy(v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("precision : ",precision(v_predict[1], v_predict[3]), sep = ""))
  #print (paste ("recall : ", recall(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("sensibility : ", sensibility(v_predict[1], v_predict[4]), sep = ""))
  print (paste ("specificity : ", sensibility(v_predict[2], v_predict[3]), sep = ""))
  print (paste ("BCR (balanced classification rate) : ", BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("BER (balanced error rate) : ", 1 - BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  print (paste ("MCC (Matthew) : ", MCC (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
  return (v_predict)
}



################
#     MAIN     #
################
args <- commandArgs(TRUE)
p_perf = args[1]

#p_perf = "c://Users/aborr/research/ILS/HERG/comparison_study/pred/pred_all"

d_perf = read.csv(p_perf, sep = "\t")
d_perf = na.omit(d_perf)

l_pred = qualityPredictList(d_perf$pred, d_perf$aff)

