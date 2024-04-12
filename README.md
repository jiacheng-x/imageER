# imageER
This package provides implementation of methods introduced by Qiu (1998, 2004, 2009) with a focus on Qiu (1998). The paper introduced a method to detect edges/jumps in image and use the edge information to denoise/recover image, such that the edge locations would not be blurred. Additionally, **imageER** provides a function to select tuning parameter used during edge detection and denoising.

## Installation

To install **imageER**, run the following command in your R console:
```R
library(devtools)
install_github("jiacheng-x/imageER")
```

## Example
This is a basic example which shows you how to use available functions on single channel square images
```R
library(imageER)
#Generate a surface with a circular jump in the middle
n <- 100
x <- y <- (1:n)/n
f <- matrix(0, n , n)

for (i in 1:n){
  for (j in 1:n){
    if ((x[i]-0.5)^2+(y[j]-0.5)^2 > 0.25^2){
      f[i,j] <- -2*(x[i]-0.5)^2-2*(y[j]-0.5)^2
    } else{
      f[i,j] <- -2*(x[i]-0.5)^2-2*(y[j]-0.5)^2+1
    }
  }
}

#Adding normal noise
sigma <- 1/2 * sd(as.vector(f))
set.seed(1234)
noise <- matrix(rnorm(n*n,0,sigma), n, n)
z <- f + noise

#Detecting Edge
edge_map <- edge_detect(z, k=5, h=4, alpha_n=0.05)


#Recover/Denoise Surface
recover_surface(z, k1=5, edge=edge_map)

#Select optimal parameters
k <- c(3, 5, 7, 9, 11)
k1 <- c(5, 10, 15, 20, 25)

par_select(k, k1, z, h = 3, alpha_n = 0.05, alpha_jp = 0.05, omega = 0.5)
```

When you have single channel non-suqare images, use `detection_recover_blocking()` instead.

```R
#Detection and recover
result = detection_recover_blocking(z, h=4, k=5, k1=5, block_size=64)
```

## Reference
- Qiu, Peihua. “Discontinuous Regression Surfaces Fitting.” The Annals of Statistics 26, no. 6 (1998): 2218–45.
- Qiu, Peihua. “The Local Piecewisely Linear Kernel Smoothing Procedure for Fitting Jump Regression Surfaces.” Technometrics 46, no. 1 (February 1, 2004): 87–98. [https://doi.org/10.1198/004017004000000149](https://doi.org/10.1198/004017004000000149).

- Qiu, Peihua. “Jump-Preserving Surface Reconstruction from Noisy Data.” Annals of the Institute of Statistical Mathematics 61, no. 3 (September 1, 2009): 715–51. [https://doi.org/10.1007/s10463-007-0166-9](https://doi.org/10.1007/s10463-007-0166-9).

## Acknowledgments

I would like to extend my thanks to Prof. Peihua Qiu, who provided invaluable assistance when I first started on implementing his method. His encouragement and permission to develop and share this package on GitHub have been fundamental to this project.

Additionally, I am deeply grateful to my advisor Prof. Richard Charnigo for his guidance during the development of the tuning parameter selection criterion. His insightful suggestions during the debugging process have significantly enhanced the functionality and reliability of this package.


