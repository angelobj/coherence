# git push -u origin main

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
other<-function(){}
