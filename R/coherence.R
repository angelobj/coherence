#' Compute Magnitude-Squared Coherence based on 'gsignal' package
#'
#' Computes the magnitude-squared coherence between two signals using the `gsignal::mscohere()` function.
#' This function estimates how well the two signals correspond at each frequency.
#'
#' @param x A data frame with exactly two numeric columns representing the signals to estimate magnitude.squared coherence from.
#' @param window Integer specifying the length of the Hann window for segmenting the signals. Default is 300.
#' @param overlap Numeric value (between 0 and 1) specifying the fraction of overlap between adjacent segments. Default is 0.5 (50% overlap).
#' @param samplingFreq Numeric value specifying the sampling frequency in Hz. Default is 1000.
#' @param nfft Integer specifying the number of FFT points used in computation. Default is 2000.
#' @param plot Logical indicating whether to generate a plot of the coherence spectrum. Default is `TRUE`.
#' @param title Optional character string specifying the plot title. If `NULL`, no title is displayed.
#' @param xlab Optional character string specifying the x-axis label. Default is `"Frequencies (Hz)"`.
#' @param ylab Optional character string specifying the y-axis label. Default is `"Magnitude-Squared Coherence"`.
#' @param out Character string specifying the output format:
#'   \item{"data"}{(Default) Returns a data frame containing the coherence results.}
#'   \item{"all"}{Returns a list containing both the coherence data and the plot.}
#'
#' @return Depending on the value of `out`:
#'   \item{data}{A data frame with:
#'     \item{Freqs}{Frequencies in Hz}
#'     \item{Cxy}{Magnitude-squared coherence values at each frequency}
#'   }
#'   \item{all}{A list containing:
#'     \item{coh}{Coherence results as a data frame}
#'     \item{plot}{Coherence plot (if `plot = TRUE`)}
#'   }
#'
#' @examples
#' # Generate simulated signals
#' signals <- simulateSignals(samplingFreq = 1000, duration = 8, target_freq = 10, noise_level = 0.5, plot = FALSE)
#'
#' # Compute coherence and return data only
#' coh_result <- coherence(signals[, c("signal1", "signal2")], samplingFreq = 1000, out = "data", plot = FALSE)
#' head(coh_result)
#'
#' # Compute coherence and return both data and plot
#' coh_result <- coherence(signals[, c("signal1", "signal2")], samplingFreq = 1000, out = "all", plot = TRUE)
#' print(coh_result$plot)
#'
#' # Custom plot labels
#' coherence(signals[, c("signal1", "signal2")], samplingFreq = 1000, plot = TRUE,
#'           title = "Coherence Analysis", xlab = "Frequency (Hz)", ylab = "Coherence Magnitude")
#'
#' @export
coherence <- function(x, window = 300, overlap = 0.5, samplingFreq = 1000, nfft = 2000,plot=T,
                      title=NULL,xlab=NULL,ylab=NULL,
                      out='data') {
  if (!requireNamespace("gsignal", quietly = TRUE)) {
    stop("The 'gsignal' package is required but not installed. Install it using install.packages('gsignal').")
  }
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (ncol(x) != 2) stop("x must have two columns")

  coh_data<<-gsignal::mscohere(
    as.matrix(x),
    window = gsignal::hann(window),
    overlap=overlap*length(gsignal::hann(window))/samplingFreq,
    nfft = nfft,
    fs = samplingFreq
  ) %>% do.call(bind_cols,.)  %>% as.data.frame() %>% setNames(c('Freqs','Cxy'))

  out_plot <- if (plot) {

    ggplot2::ggplot(coh_data, ggplot2::aes(x = Freqs, y = Cxy)) +
      ggplot2::geom_line(color = "blue", linewidth = 1) +  # Ensuring proper mapping of data
      ggplot2::labs(
        title = if (!is.null(title)) title else NULL,
        x = if (!is.null(xlab)) xlab else "Frequencies (Hz)",
        y = if (!is.null(ylab)) ylab else "Magnitude-Squared Coherence"
      ) +
      ggplot2::theme_minimal(base_size = 15) +  # Clean plot style
      ggplot2::theme(
        legend.position = "none",
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black")
      )
  } else {
    NULL
  }

  # Ensure the plot is displayed when `plot = TRUE`
  if (!is.null(out_plot)) print(out_plot)

  # Return results
  if(is.null(out)||out=='data'){
    return(coh_data)
  }else if(!is.null(out)&&out=='all'){
    return(list(coh = coh_data, plot = out_plot))
  }else{
    stop("out must be 'data' or 'all'")
  }
}

#' Generate Simulated Signals and Perform FFT Analysis
#'
#' Simulates two signals: one as a pure sine wave and the other as a noisy version of the sine wave.
#' Computes the Fast Fourier Transform (FFT) of both signals and optionally generates plots of the time and frequency domain representations.
#'
#' @param samplingFreq Numeric value specifying the sampling frequency of the signals in Hz. Default is 1000 Hz.
#' @param duration Numeric value specifying the duration of the signals in seconds. Default is 8 seconds.
#' @param target_freq Numeric value specifying the target frequency of the sine wave in Hz. Default is 10 Hz.
#' @param noise_level Numeric value specifying the standard deviation of the Gaussian noise added to the second signal. Default is 0.5.
#' @param plot Logical indicating whether to generate plots of the signals. Default is TRUE.
#' @param filter_freqs Optional numeric value specifying a frequency threshold (in Hz). If provided, only frequency components above this value will be included in the FFT analysis.
#' @param seed Integer specifying the random seed for reproducibility. Default is 1.
#' @param out Character string specifying the output format:
#'   \item{"data"}{(Default) Returns only the simulated signals as a data frame.}
#'   \item{"all"}{Returns a list containing the simulated signals, FFT analysis, and plots.}
#'
#' @return Depending on the value of `out`:
#'   \item{data}{A data frame containing the simulated signals (time, signal1, signal2).}
#'   \item{all}{A list containing:
#'     \item{signals}{Data frame with the simulated signals}
#'     \item{fft}{FFT results for both signals}
#'     \item{plot}{Time and frequency domain plots (if `plot = TRUE`)}
#'   }
#'
#' @examples
#' # Generate signals and return only the data
#' result <- simulateSignals()
#'
#' # Generate signals and return all outputs
#' result <- simulateSignals(out = "all")
#'
#' @export
simulateSignals <- function(samplingFreq = 1000, duration = 8, target_freq = 10, noise_level = 0.5, plot = TRUE,
                            filter_freqs = NULL,seed=1,out='data') {
  set.seed(seed)  # Ensures reproducibility

  # Ensure necessary packages are available
  required_pkgs <- c("ggplot2", "magrittr", "tidyr", "ggpubr", "dplyr")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(paste("The following required packages are missing:", paste(missing_pkgs, collapse = ", "),
               "\nInstall them using install.packages()."))
  }

  # Time vector
  Time <- seq(0, duration, by = 1 / samplingFreq)

  # Generate the base signal (pure sine wave at target frequency)
  signal1 <- sin(2 * pi * target_freq * Time)

  # Generate the second signal (sine wave + broadband noise)
  signal2 <- signal1 + rnorm(length(Time), mean = 0, sd = noise_level)  # Ensures correlation only at target frequency
  signal1<-signal1+rnorm(length(Time), mean = 0, sd = noise_level)

  # Store in a data frame
  signals <- data.frame(Time, signal1, signal2)

  ### Compute FFT for both signals ###
  fft_data <- lapply(c('signal1', 'signal2'), function(x) {
    fft_analysis(signals[[x]], samplingFreq, x, filter_freqs)
  }) %>% do.call(rbind, .) %>%
    mutate(Signal = factor(case_when(
      Signal == 'signal1' ~ 'Signal 1',
      Signal == 'signal2' ~ 'Signal 2'
    ), levels = c('Signal 1', 'Signal 2'), ordered = TRUE))


  ### Generate Plots ###
  out_plot <- if (plot) {

    # Prepare long-format data for time-domain plotting
    long_signals <- signals %>%
      tidyr::pivot_longer(!Time, names_to = 'Signal', values_to = 'Amplitude') %>%
      mutate(Signal = factor(case_when(
        Signal == 'signal1' ~ 'Signal 1',
        Signal == 'signal2' ~ 'Signal 2'
      ), levels = c('Signal 2', 'Signal 1'), ordered = TRUE))

    # Time-domain plot
    time_plot <- ggplot2::ggplot() +
      ggplot2::geom_line(aes(x = Time, y = Amplitude, colour = Signal),
                         data = long_signals, alpha = 0.9) +  # Corrected data source
      ggplot2::scale_color_manual(values = c('Signal 1' = 'red', 'Signal 2' = '#0f5dd1')) +
      ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::theme(
        legend.position = 'top',
        legend.title = element_blank(),
        strip.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        text = element_text(size = 20)
      )

    # Frequency-domain (FFT) plot
    fft_plot <- ggplot2::ggplot() +
      ggplot2::geom_line(aes(x = Freqs, y = Amplitude, colour = Signal), data = fft_data, alpha = 0.9) +
      ggplot2::scale_color_manual(values = c('Signal 1' = 'red', 'Signal 2' = '#0f5dd1')) +
      ggplot2::guides(color = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::theme(
        strip.text.x = element_blank(),
        strip.text = element_text(color = "black", size = 15),
        legend.position = 'top',
        legend.title = element_blank(),
        strip.background = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"),
        text = ggplot2::element_text(size = 20)
      ) + ggplot2::facet_wrap(~Signal, scales = 'free_y')

    # Combine plots
    ggpubr::ggarrange(time_plot, fft_plot, common.legend = TRUE, nrow = 1)
  } else { NULL }

  if(!is.null(out_plot))
    print(out_plot)

  # Return results
  if(is.null(out)||out=='data'){
    return(signals)
  }else if(!is.null(out)&&out=='all'){
    return(list(signals = signals, fft = fft_data, plot = out_plot))
  }else{
    stop("out must be 'data' or 'all'")
  }
}
#' Perform Fast Fourier Transform (FFT) Analysis for easy plotting
#'
#' Computes the Fast Fourier Transform (FFT) of a given signal and returns its frequency spectrum.
#'
#' @param signal Numeric vector representing the time-domain signal to analyze.
#' @param samplingFreq Numeric value specifying the sampling frequency of the signal in Hz.
#' @param label Character string representing the label for the signal (e.g., "Signal 1", "Signal 2").
#' @param filter_freqs Optional numeric value specifying a frequency threshold (in Hz). If provided, only frequency components above this value will be returned.
#'
#' @return A data frame containing:
#'   \item{Freqs}{Frequencies in Hz}
#'   \item{Amplitude}{Magnitude of the FFT at each frequency}
#'   \item{Signal}{Label for the signal}
#'
#' @export
fft_analysis <- function(signal, samplingFreq, label, filter_freqs = NULL) {
  N <- length(signal)  # Number of samples

  fft_values <- abs(fft(signal)) / N  # Normalize the FFT output
  Amplitude <- fft_values[1:(N/2)]  # Keep only the first half (positive frequencies)
  Freqs <- seq(0, samplingFreq / 2, length.out = length(Amplitude))  # Match Amplitude length

  # Create FFT data frame
  fft_data <- data.frame(Freqs, Amplitude, Signal = rep(label, times = length(Amplitude)))

  # Apply frequency filtering if filter_freqs is provided
  if (!is.null(filter_freqs) && is.numeric(filter_freqs) && filter_freqs > 0) {
    fft_data <- fft_data %>% dplyr::filter(Freqs > filter_freqs)
  }

  return(fft_data)
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
