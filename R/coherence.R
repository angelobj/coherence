#' Magnitude Square Coherence Estimation
#'
#' Computes the magnitude square coherence between two signals using the `gsignal::mscohere` function.
#' This function estimates the coherence between two signals based on the Hann window and given parameters.
#'
#' @param x A data frame with two numeric columns representing the two signals to analyze.
#' @param window Integer specifying the length of the Hann window for segmenting the signals. Default is 300.
#' @param overlap Numeric value specifying the number of overlapping samples between adjacent segments.
#'               Should be a fraction of the window length, typically between 0 and `window`. Default is calculated as `cohOverlap * length(hann(window)) / samplingFreq`.
#' @param samplingFreq Sampling frequency of the signals in Hz (samples per second). Default is 1000.
#' @param nfft Integer specifying the number of points used in the FFT computation. Should be a power of 2 and at least as large as the window length. Default is 2000.
#'
#' @return A list containing:
#'   \item{f}{Frequency vector}
#'   \item{Cxy}{Coherence estimate between the two signals}
#'
#' @export
coherence <- function(x, window = 300, overlap = 0.5, samplingFreq = 1000, nfft = 2000) {
  if (!requireNamespace("gsignal", quietly = TRUE)) {
    stop("The 'gsignal' package is required but not installed. Install it using install.packages('gsignal').")
  }
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (ncol(x) != 2) stop("x must have two columns")

  gsignal::mscohere(
    x %>% as.matrix(),
    window = gsignal::hann(window),
    overlap = overlap * window,  # Now using a proportion of the window
    nfft = nfft,  # Ensure it is passed correctly
    fs = samplingFreq
  )
}

#' Fisher z-transformation
#'
#' Fisher's z-transformation converts Pearson's r to the normally distributed variable z.
#' The formula for the transformation is:
#' \deqn{z_r = \tanh^{-1}(r) = \frac{1}{2} \log \left( \frac{1+r}{1-r} \right)}
#'
#' @param x Numeric vector or matrix containing Pearson correlation coefficients (r) to be transformed.
#' @return A numeric vector or matrix with transformed values.
#' @export
fisherz <- function(x) {
  # Check if DescTools is installed
  if (!requireNamespace("DescTools", quietly = TRUE)) {
    stop("The 'DescTools' package is required but not installed. Install it using install.packages('DescTools').")
  }

  # Perform Fisher z-transformation
  DescTools::FisherZ(x)
}


#' Time Frequency Representation (based on Dr. Christopher Laine matlab code).
#' @param x Vector
#' @param freqs Frequencies
#' @param Fs Sampling Frequency
#' @param width Number of cycles
#' @examples
#' x<-sin(seq(0,pi,by=0.1));plot(x)
#' tfr(x,freqs=seq(8,50,by=2))
#' @export
tfr<-function(x=NULL,freqs=seq(8,50,by=2),Fs=1000,width=7){
  #if(!"gsignal" %in% rownames(installed.packages())){stop("Package 'gsignal' must be installed")}
  if(is.null(x)) stop("Data must be provided")
  if(class(x)%in%c('matrix','data.frane')) stop('Data must be a data.frame or matrix')
  if(class(x)=='matrix') warning('Data must be a data.frame and is converted in the process.')
  if(class(x)=='matrix') x<-as.data.frame(x)
  complex_tfr<-sapply(freqs,function(f){
    print(paste0("Time frequency representation of ",f," (",round(length(f)/which(f==freqs)*100),"%)"))
    time<-(1:ceiling(width/f*Fs))/Fs;
    gaussvect<-gsignal::gausswin(length(time))
    gaussvect<-gaussvect/sum(gaussvect);
    sinewave<-(exp(2*sqrt(as.complex(-1))*pi*f*time));
    morlet<-sinewave*gaussvect
    complex_tfr<-(2)*gsignal::conv(x,morlet,'same');
  },USE.NAMES = T)
  return(complex_tfr)
}

#' Phase angle between complex valued vectors (based on Dr. Christopher Laine matlab code).
#' @param x complex valued vector
#' @param y complex valued vector
#' @examples
#' x<-c(2 +   2i,3 +   2i,4 +   7i,5 +   7i)
#' y<-c(3 +   2i,4 +   2i,5 +   7i,6 +   7i)
#' msc(x,y)
#' @export
msc<-function(x,y){
  mscoh<-(abs(sum(x*Conj(y))^2))/(sum(abs(x)^2)*sum(abs(y)^2))
  coh<-(abs(sum(x*Conj(y))))/sqrt((sum(abs(x)^2)*sum(abs(y)^2)))
  icoh<-Im(((sum(x*Conj(y))))/sqrt((sum(abs(x)^2)*sum(abs(y)^2))))
  phase_angle<-Arg(((sum(x*Conj(y))))/sqrt((sum(abs(x)^2)*sum(abs(y)^2))));
  return(mscoh)
}

#' Coherence (rho) to Z-score conversion(based on Dr. Christopher Laine matlab code).
#' @param x Raw coherence values
#' @param datalen The total length (in samples) of the original signals that coherence was calculated on (samples, not segments or seconds)
#' @param window Vector for tapering segments: rectwin(1024) or gausswin(2048);
#' @param overlap Fraction of the window to overlap, e.g. 0 or 0.5
#' @export
coh2z<-function(x,datalen,window,overlap){
  winlen=length(window);
  Lvec<-if(overlap==0){floor(datalen/winlen)}else{
    L1=floor(datalen/winlen);
    L2 = floor( (L1-1)/(1-overlap))+1;
    k1=floor(overlap*winlen);
    k2=ceiling((1-overlap)*winlen);
    wprime_numerator=sum(window[1:k1]*window[(k2+1):length(window)]);
    wprime_denominator= sum(window^2);
    wprime=(wprime_numerator/wprime_denominator)^2;
    w=1/(1+2*wprime);
    Lvec=w*L2;
    Fz = atanh(sqrt(x));
    Zscores = Fz/sqrt(1/(2*Lvec));
    return(Zscores)
  }
}
