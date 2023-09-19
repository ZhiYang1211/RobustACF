library(sazedR)
library(forecast)


# seasonal_no_trend.csv
data_csv2<-read.csv('data\\seasonal_no_trend.csv')
ts2 = ts(data_csv2['value'])
sazed(ts2)
findfrequency(ts2)

# seasonal_no_trend.csv, outlier
data_csv2<-read.csv('data\\seasonal_no_trend.csv')
ts2 = ts(data_csv2['value'])
ts2[71] = 2
sazed(ts2)
findfrequency(ts2)

# seasonal.csv
data_csv2<-read.csv('data\\seasonal.csv')
ts2 = ts(data_csv2['value'])
sazed(ts2)
findfrequency(ts2)

# seasonal.csv, outliers
data_csv2<-read.csv('data\\seasonal.csv')
ts2 = ts(data_csv2['value'])
ts2[76] = 5
ts2[77] = 5
ts2[78] = 5
ts2[161] = -2
ts2[162] = -2
sazed(ts2)
findfrequency(ts2)

# tsdl_271.csv, T = 365
data_csv2<-read.csv('data\\tsdl_271.csv')
ts2 = ts(data_csv2['value'])
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_271.csv, T = 365, two outliers
data_csv2<-read.csv('data\\tsdl_271.csv')
ts2 = ts(data_csv2['value'])
ts2[601] = -40
ts2[1101] = 30
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_331.csv, T = 365
data_csv2<-read.csv('data\\tsdl_331.csv')
ts2 = ts(data_csv2['value'])
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_331.csv, T = 365, two outliers
data_csv2<-read.csv('data\\tsdl_331.csv')
ts2 = ts(data_csv2['value'])
ts2[201] = -30
ts2[801] = 30
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_91_weekly.csv, T = 52
data_csv2<-read.csv('data\\tsdl_91_weekly.csv')
ts2 = ts(data_csv2['value'])
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_91_weekly.csv, T = 52, three outliers
data_csv2<-read.csv('data\\tsdl_91_weekly.csv')
ts2 = ts(data_csv2['value'])
ts2[61] = 10
ts2[401] = 35
ts2[301] = 30
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_92_weekly.csv, T = 52
data_csv2<-read.csv('data\\tsdl_92_weekly.csv')
ts2 = ts(data_csv2['value'])
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_92_weekly.csv, T = 52, an outlier
data_csv2<-read.csv('data\\tsdl_92_weekly.csv')
ts2 = ts(data_csv2['value'])
ts2[111] = 2.5
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_244.csv, T = 12
data_csv2<-read.csv('data\\tsdl_244.csv')
ts2 = ts(data_csv2['value'])
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_244.csv, T = 12, an outlier
data_csv2<-read.csv('data\\tsdl_244.csv')
ts2 = ts(data_csv2['value'])
ts2[26] = 300
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_358.csv, T = 12
data_csv2<-read.csv('data\\tsdl_358.csv')
ts2 = ts(data_csv2['value'])
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)

# tsdl_358.csv, T = 12, two outliers
data_csv2<-read.csv('data\\tsdl_358.csv')
ts2 = ts(data_csv2['value'])
ts2[16] = 1000
ts2[66] = 1100
x = 1:length(ts2)
plot(x, ts2, type = 'l')+
  grid()
sazed(ts2)
findfrequency(ts2)
