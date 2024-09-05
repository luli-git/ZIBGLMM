proc import out=simulationdata
    datafile="/home/u63419104/simulations_data.csv" 
    dbms=csv
    replace;
    getnames=YES;
run;
 
/* ods exclude all;                  */
ods output FitStatistics=Output1;   
ods output ParameterEstimates=Output3;   
ods output AdditionalEstimates=Output5; 
title "Probit model";
title2 "Zero-Inflated Probit Model for Two Groups";
data filtered_data;
    set simulationdata;
    if metaanalysisid < 100;  /*For illustration, you can use the first 100, to conduct full experiment, change it to if metaanalysisid > 0 */
run;
proc nlmixed data=filtered_data fd df=1000 gtol=1e-10 tech=CONGRA;
 parms logitpi=-1 beta1=-2, sig1=-1, beta2=-2, sig2=-1 fz=0;
 where metaanalysisid not in (36,192,528,535,626,719,820,874,903,1042,1101,1156,1196,1229,1249,1262,1337,
 1338,1447,1472,1499,1861,1880,1932,1951,2021,2104,2105,2260,2502,2537,2554,2596,2626,2709,2722,2727,
 2754,2781,2793,2803,2903,3138,3743,4033,4382,4443,4651,4669,4688,4730,4768,4868,4999,5005,5008,5050,
 5103,5141,5161,5227,5239,5286,5295,5309,5337,5358,5392,5450,5482,5515,5530,5559,5696,5895,5956,5978,
 7080,7112,7121,7131,7189,7252,7481,7530,7595,7639,7690,7971,8003,8047,8077,8092,8191,8211,8237,8263,
 8306,8312,8315,8318,8321,8424,8428,8440,8474,8504,8560,8669,8717,8764,8773,8775,8797,8829,8874,8898,
 8910,9015,9113,9152,9236,9338,9436,9464,9521,9644,9874,9886,9905,9976,10159,10186,10206,10325,10412,
 10439,10514,10578,10722,10733,10750,10799,10801,10839,10866,10869,10887,10946,11064,11078,11091,11111,
 11166,11167,11255,11261,11304,11341,11447,11476,11477,11519,11525,11624,11638,11648,11675,11690,11739,
 11796,11820,11838,11890,11977,11999,12082,12124,12373,12393,12474,12493,12593,12611,12709,12733,13059,
 13061,13087,13168,13427,13479,13539,13571,13601,13609,13697,13796,13843,13889,13957,14063,14066,14097,
 14124,14145,14198,14211,14232,14270,14281,14316,14349,14371,14393,14402,14437,14452,14456,14480,14516,
 14546,14585,14598,14618,14697,14790,14795,14810,14828,14829,14832,14845,14890,14922,14937,15088,15197,
 15234,15651,15830,15924,16039,16044,16271,16365,16367,16371,16372,16416,16460,16461,16623,16657,16714,
 16795,16817,16924,17015,17030,17062,17082,17097,17143,17148,17165,17170,17226,17264,17282,17283,17291,
 17298,17344,17423,17448,17468,17496,17502,17533,17570,17599,17616,17687,17701,17708,17709,17760,17767,
 17817,17865,17891,17900,17915,17952,17962);

 by metaanalysisid; 
 pi = exp(logitpi)/(1+exp(logitpi));
 rho=(exp(2*fz)-1)/(exp(2*fz)+1);
 beta1i = beta1 + mu1 ;
 pred1i = exp(beta1i)/(1+exp(beta1i));
 beta2i = beta2 + mu2 ;
 pred2i = exp(beta2i)/(1+exp(beta2i));
 C = (16*sqrt(3)/(15*3.1415));
 P1 = 1/(1+exp(-beta1/sqrt(1+C**2*exp(sig1*2)))); * marginal probablity in group 1 in none ZI part;
 P2 = 1/(1+exp(-beta2/sqrt(1+C**2*exp(sig2*2)))); * marginal probablity in group 2 in none ZI part;
 if y1 = 0 and y2 = 0 then do;
  	logll = log(pi+(1-pi)*exp(n1*log(1-pred1i)+n2*log(1-pred2i) ) );
/*   	logll = log(pi+(1-pi)* (1-pred1i)**n1 * (1-pred2i)**n2  ); */
/* 15  */
  	if pi < 0.000001 then logll = n1*log(1-pred1i)+n2*log(1-pred2i);
 end;
 if y1 > 0 or y2>0 then logll = log(1-pi)+ y1*log(pred1i)+(n1-y1)*log(1-pred1i) 
                     + y2*log(pred2i)+(n2-y2)*log(1-pred2i) ;
 model yy ~ general(logll);
 random mu1 mu2 ~normal([0, 0], [exp(2*sig1), rho*exp(sig1+sig2), exp(2*sig2)]) subject=site;
 estimate "Zero-Inflation Proportion pi" pi;
 estimate "logRR" log(P2/P1);
 estimate "RR" P2/P1;
run;
/* ods exclude none;                   */
/* proc print data=Output5 noobs; */
/* run; */
PROC EXPORT DATA=Output1
    OUTFILE='/home/u63419104/zibglmm_1.csv'
    DBMS=CSV
    REPLACE;
RUN;

/* proc print data=Output3 noobs; */
/*   var metaanalysisid Parameter Estimate Lower Upper; */
run;
PROC EXPORT DATA=Output3
    OUTFILE='/home/u63419104/zibglmm_2.csv'
    DBMS=CSV
    REPLACE;
RUN;

PROC EXPORT DATA=Output5
    OUTFILE='/home/u63419104/zibglmm_3.csv'
    DBMS=CSV
    REPLACE;
RUN;
