This dataset is obtained by extracting all open-source univariate time series
of packages listed in the "Time Series Data" section of the CRAN Task View on Time
Series Analysis, a comprehensive collection of R packages on time series analysis.

The time series included here span across a wide variety of application
domains, ranging from economic indicators, such as employment rates or retail sales,
to environmental measurements, such as pollution levels or the number of sunspots.

* Only time series of the R object class "ts" are considered.
* Time-series data without specified frequency value are abandoned
* Time-series data with specified frequency value yet without any real seasonal pattern are also abandoned
* Total number of time-series data is 81, all in csv format(1st column - time index, 2nd column - value)
* File naming convention : [DataName]_[R-pkgName]_FREQ_[Frequency].csv