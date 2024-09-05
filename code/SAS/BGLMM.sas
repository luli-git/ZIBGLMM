proc import out=simulationdata
    datafile="/home/u63419104/simulations_data.csv" 
    dbms=csv
    replace;
    getnames=YES;
run;
 
ods exclude all;                 
ods output FitStatistics=Output1;   
ods output ParameterEstimates=Output3;   
ods output AdditionalEstimates=Output5; 
title "Probit model";
title2 "Zero-Inflated Probit Model for Two Groups";
data filtered_data;
    set simulationdata;
    if metaanalysisid < 100;   /*For illustration, you can use the first 100, to conduct full experiment, change it to if metaanalysisid > 0 */
run;
proc nlmixed data=filtered_data fd df=1000 gtol=1e-10 tech=CONGRA;
 parms beta1=-2, sig1=-1, beta2=-2, sig2=-1 fz=0;
 where metaanalysisid not in (454,2491,12671);
 
 by metaanalysisid;
 rho=(exp(2*fz)-1)/(exp(2*fz)+1);
 beta1i = beta1 + mu1 ;
 pred1i = exp(beta1i)/(1+exp(beta1i));
 beta2i = beta2 + mu2 ;
 pred2i = exp(beta2i)/(1+exp(beta2i));
 C = (16*sqrt(3)/(15*3.1415));
 P1 = 1/(1+exp(-beta1/sqrt(1+C*C*exp(sig1*2)))); * marginal probablity in group 1 in none ZI part;
 P2 = 1/(1+exp(-beta2/sqrt(1+C*C*exp(sig2*2)))); * marginal probablity in group 2 in none ZI part;
 logll = y1*log(pred1i)+(n1-y1)*log(1-pred1i) + y2*log(pred2i)+(n2-y2)*log(1-pred2i) ;
 model yy ~ general(logll);
 random mu1 mu2 ~normal([0, 0], [exp(2*sig1), rho*exp(sig1+sig2), exp(2*sig2)]) subject=site;
 estimate "logRR" log(P2/P1);
 estimate "RR" P2/P1;
run;
ods exclude none;                  
 
PROC EXPORT DATA=Output1
    OUTFILE='/home/u63419104/bglmm_1.csv'
    DBMS=CSV
    REPLACE;
RUN;

proc print data=Output3 noobs;
  var metaanalysisid Parameter Estimate Lower Upper;
run;
PROC EXPORT DATA=Output3
    OUTFILE='/home/u63419104/bglmm_2.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=Output5
    OUTFILE='/home/u63419104/bglmm_3.csv'
    DBMS=CSV
    REPLACE;
RUN;
