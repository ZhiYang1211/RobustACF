library(AER)
library(fpp)
library(fpp2)
library(astsa)
library(expsmooth)
library(fma)
library(TSA)

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
              a10='fpp2', ausbeer='fpp2', auscafe='fpp2', #ausgdp='fpp2',
              austourists='fpp2', debitcards='fpp2', elecequip='fpp2',
              gasoline='fpp2', h02='fpp2', hyndsight='fpp2', qauselec='fpp2',
              qcement='fpp2', qgas='fpp2', usmelec='fpp2',
              airmiles='TSA', beersales='TSA',
              co2='TSA', flow='TSA', hours='TSA', milk='TSA', JJ='TSA', oilfilters='TSA',
              prescrip='TSA', prey.eq='TSA', retail='TSA', tempdub='TSA',
              wages='TSA', winnebago='TSA')
for (data.name in names(data_list)) {
     pkg <- data_list[[data.name]]
     data(list=c(data.name), package=pkg)
     data <- data.name %>% get()
     freq <- as.character(frequency(data))
     x <- as.vector(data)
     df <- data.frame(id=c(1:length(x)), value=x)
     filename <- paste0(data.name, '_', pkg, '_FREQ_', freq, '.csv')
     write.csv(df, paste0('C:/Users/i326292/Desktop/crandata/', filename), row.names=FALSE)
}
