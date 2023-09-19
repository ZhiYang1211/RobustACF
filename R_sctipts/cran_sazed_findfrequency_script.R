library(AER)
library(fpp)
library(fpp2)
library(astsa)
library(expsmooth)
library(fma)
library(TSA)
library(sazedR)
library(forecast)

data_list = c(BondYield='AER', DutchSales='AER', UKNonDurables='AER',
              birth='astsa', cmort='astsa',
              flu='astsa', gas='astsa', hor='astsa', part='astsa',
              prodn='astsa', qinfl='astsa', qintr='astsa', rec='astsa',
              so2='astsa', soi='astsa', sunspotz='astsa', tempr='astsa',
              unemp='astsa', UnempRate='astsa',
              bonds='expsmooth', cangas='expsmooth', enplanements='expsmooth',
              frexport='expsmooth', mcopper='expsmooth', ukcars='expsmooth',
              usgdp='expsmooth', utility='expsmooth',
              vehicles='expsmooth', visitors='expsmooth',
              airpass='fma', beer='fma', bricksq='fma', condmilk='fma', dole='fma',
              elec='fma', fancy='fma', hsales='fma', hsales2='fma',
              invent15='fma',labour='fma', milk='fma', motion='fma',
              pigs='fma', plastics='fma', pollution='fma', qelec='fma',
              qsales='fma', shampoo='fma', ukdeaths='fma',
              usdeaths='fma', uselec='fma', writing='fma',
              cafe='fpp', euretail='fpp',
              a10='fpp2', ausbeer='fpp2', auscafe='fpp2', ausgdp='fpp2',
              austourists='fpp2', debitcards='fpp2', elecequip='fpp2',
              gasoline='fpp2', h02='fpp2', hyndsight='fpp2', qauselec='fpp2',
              qcement='fpp2', qgas='fpp2', usmelec='fpp2',
              airmiles='TSA', beersales='TSA',
              co2='TSA', flow='TSA', hours='TSA', milk='TSA', JJ='TSA', oilfilters='TSA',
              prescrip='TSA', prey.eq='TSA', retail='TSA', tempdub='TSA',
              wages='TSA', winnebago='TSA')

#data_list = c(BondYield='AER', DutchSales='AER', UKNonDurables='AER')

# get periods detected by sazed and findfrequency
freqs = c()
sazed_freqs = c()
findf_freqs = c()
for (data.name in names(data_list)) {
    pkg <- data_list[[data.name]]
    data(list=c(data.name), package=pkg)
    data <- data.name %>% get()
    freqs <- append(freqs, frequency(data))
    x <- as.vector(data)
    
    sazed_fre = sazed(x)
    findf_fre = findfrequency(x)
    sazed_freqs <- append(sazed_freqs, sazed_fre)
    findf_freqs <- append(findf_freqs, findf_fre)
}

# cal error
error_measure_num<-function(periods_detect, periods_real, errors){
  if (length(periods_detect) != length(periods_real)){
    return(-1)
  }
  abs_diff <- abs(floor(periods_detect) - floor(periods_real))
  abs_diff_ratio <- abs_diff / floor(periods_real)
  num <- length(periods_detect)
  errors_sort <- sort(errors)
  num_error = length(errors_sort)
  error_distr <- integer(num_error+1)
  for (i in 1:num){
    for (j in 1:num_error){
      if (abs_diff_ratio[i] <= errors_sort[j]){
        error_distr[j] = error_distr[j] + 1
        break
      }
    }
    if (abs_diff_ratio[i] > errors_sort[num_error]){
      error_distr[num_error+1] = error_distr[num_error+1] + 1
    }
  }

  return(error_distr)
}


errors = c(0, 0.05, 0.1, 0.15)

sazed_error_distr = error_measure_num(sazed_freqs, freqs, errors)
findf_error_distr = error_measure_num(findf_freqs, freqs, errors)




  