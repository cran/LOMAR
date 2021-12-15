#' pssk
#' 
#' Compute the persistence scale-space kernel on persistence diagrams. 
#' Reference: Jan Reininghaus, Stefan Huber, Ulrich Bauer, and Roland Kwitt. A stable multi-scale kernel for topological machine learning. In Proceedings of the IEEE conference on computer vision and pattern recognition (CVPR), pages 4741–4748, 2015.
#'
#' @param Dg1 a persistence diagram as a n1 x 3 matrix where each row is a topological feature
#' and the columns are dimension, birth and death of the feature.
#' @param Dg2 another persistence diagram as a n2 x 3 matrix
#' @param sigma kernel bandwidth
#' @param dimensions vector of the dimensions of the topological features to consider, if NULL (default) use all available dimensions
#' @return kernel value
#' @examples 
#' D1 <- matrix(c(0,0,0,1,1,0,0,0,1.5, 3.5,2,2.5,3, 4, 6), ncol = 3, byrow = FALSE)
#' D2 <- matrix(c(0,0,1,1,0, 0, 1.2, 2, 1.4, 3.2,4.6,6.5), ncol = 3, byrow = FALSE)
#' K <- pssk(Dg1 = D1, Dg2 = D2, sigma = 1)
#' @export

pssk <- function(Dg1 = NULL, Dg2 = NULL, sigma = NULL, dimensions = NULL) {
  if (sigma <= 0) {
    stop("ERROR: Parameter sigma must be strictly positive.\n")
  }
  if(!is.null(dimensions)) {
    Dg1 <- Dg1[Dg1[,1] %in% dimensions,, drop = FALSE]
    Dg2 <- Dg2[Dg2[,1] %in% dimensions,, drop = FALSE]
  }
  dims <- unique(c(Dg1[,1], Dg2[,1]))
  K <- 0
  for(d in dims) {
    PD1 <- Dg1[Dg1[,1] == d, 2:3, drop = FALSE]
    PD2 <- Dg2[Dg2[,1] == d, 2:3, drop = FALSE]
    n1 <- nrow(PD1)
    n2 <- nrow(PD2)
    if(n1==0 || n2==0) { next }
    k <- 0
    # Note: nested loops 25% faster than using apply()
    for(i in 1:n1) {
      p <- PD1[i,]
      for(j in 1:n2) {
        q <- PD2[j,]
        qbar <- rev(q)
        d1 <- (p[1]-q[1])^2+(p[2]-q[2])^2
        d2 <- (p[1]-qbar[1])^2+(p[2]-qbar[2])^2
        k <- k + exp(-d1/(8*sigma)) -exp(-d2/(8*sigma))
      }
    }
    k <- k/(8*pi*sigma)
    K <- K + k
  }
  return(K)
}

#' sliced_Wd
#' 
#' Compute sliced Wasserstein distance or kernel. 
#' Reference: Mathieu Carriere, Marco Cuturi, and Steve Oudot. Sliced Wasserstein kernel for persistence diagrams. In Proceedings of the 34th International Conference on Machine Learning, volume 70 of Proceedings of Machine Learning Research, pages 664–673, 2017.
#'
#' @param Dg1 a persistence diagram as a n1 x 3 matrix where each row is a topological feature
#' and the columns are dimension, birth and death of the feature.
#' @param Dg2 another persistence diagram as a n2 x 3 matrix
#' @param M number of slices (default: 10)
#' @param sigma kernel bandwidth (default: 1)
#' @param dimensions  vector of the dimensions of the topological features to consider, if NULL (default) use all available dimensions
#' @param return.dist logical (default: FALSE). Whether to return the kernel or distance value.
#' @return kernel or distance value
#' @examples 
#' D1 <- matrix(c(0,0,0,1,1,0,0,0,1.5, 3.5,2,2.5,3, 4, 6), ncol = 3, byrow = FALSE)
#' D2 <- matrix(c(0,0,1,1,0, 0, 1.2, 2, 1.4, 3.2,4.6,6.5), ncol = 3, byrow = FALSE)
#' K <- sliced_Wd(Dg1 = D1, Dg2 = D2, M = 10, sigma = 1, return.dist = TRUE)
#' @export

sliced_Wd <- function(Dg1, Dg2, M = 10, sigma = 1, dimensions = NULL, return.dist = FALSE) {
 
  if(!is.null(dimensions)) {
    Dg1 <- Dg1[Dg1[,1] %in% dimensions,, drop = FALSE]
    Dg2 <- Dg2[Dg2[,1] %in% dimensions,, drop = FALSE]
  }
  dims <- unique(c(Dg1[,1], Dg2[,1]))
  SW <- 0
  for(d in dims) {
    PD1 <- Dg1[Dg1[,1] == d, 2:3, drop = FALSE]
    PD2 <- Dg2[Dg2[,1] == d, 2:3, drop = FALSE]
    n1 <- nrow(PD1)
    n2 <- nrow(PD2)
    if(n1==0 || n2==0) { next }
    # Project diagram points on the diagonal and add to the other diagram
    L <- c(2,2) # Vector representing the diagonal
    p.PD1 <- matrix(0, ncol = 2, nrow = n1)
    for(i in 1:n1) {
      p <- PD1[i,]
      s <- (L %*% p)/(L%*%L)
      p.PD1[i,] <- s[1,1] * L
    }
    p.PD2 <- matrix(0, ncol = 2, nrow = n2)
    for(i in 1:n2) {
      p <- PD2[i,]
      s <- (L %*% p)/(L%*%L)
      p.PD2[i,] <- s[1,1] * L
    }
    PD1 <- rbind(PD1, p.PD2)
    PD2 <- rbind(PD2, p.PD1)
    
    theta <- -pi/2
    s <- pi/M
    K <- 0
    V1 <- numeric(n1+n2)
    V2 <- numeric(n1+n2)
    for(i in 1:M) {
      L <- c(cos(theta), sin(theta))
      V1 <- apply(PD1, 1, function(x) { x %*% L })
      V1 <- sort(V1)
      V2 <- apply(PD2, 1, function(x) { x %*% L })
      V2 <- sort(V2)
      K <- K + s * (sum(abs(V1 - V2)))
      theta <- theta + s
    }
    K <- K/pi
    SW <- SW + K
  }
  if (!return.dist) {
    return(exp(-SW/sigma))
  } else {
    return(SW)
  }
}

#' get_persistence_diagrams
#' 
#' Compute persistence diagrams for a list of point sets.
#' By default, compute persistent homology from the Vietoris-Rips filtration.
#' If use.dtm is TRUE, compute instead the persistent homology of the sublevel
#' set of the distance to measure evaluated over a grid.
#' 
#' @importFrom foreach %dopar%
#' @param point.sets list of point sets, each as a data frame with columns x,y,z
#' @param maxdimension maximum dimension of the homological features to be computed
#' @param maxscale limit of the Vietoris-Rips filtration
#' @param use.dtm logical (default: FALSE), whether to use the distance to measure function
#' @param m0 parameter for the dtm function
#' @param grid.by vector of space between points of the grid for the dtm function along each dimension
#' @param ncpu number of parallel threads to use for computation
#' @return a list of persistence diagrams as n x 3 matrices. Each row is a topological feature
#'  and the columns are dimension, birth and death of the feature.
#' @examples 
#' PS <- list(data.frame(x = c(2.4,-6.9,4.6,-0.7,-3.3,-4.9,-3.5,-3.5,4.2,-7),
#'                       y = c(5.7,1.9,4.8,3.4,-3,-2.1,7.2,1.8,6.1,-1.6),
#'                       z = c(2.7,-0.1,-0.7,-0.6,0.4,-1.5,-0.6,-0.9,2.2,0.7)),
#'            data.frame(x = c(0,0,3.1,-5.6,-5,-7.4,-0.7,-7.7,-6.7,4,4.2,0.2,5.8,3.9,3.9),
#'                       y = c(6.3,-6.1,-3.5,4.6,-4.1,0.3,8.8,-2.3,2.9,3.7,-1.4,-3.9,5.5,-1.2,-6.7),
#'                       z = c(-1.5,1.7,-0.4,-1.4,1.8,1.7,-0.9,-1.8,-0.5,1.7,1.3,0.5,-1.4,1.6,-0.1)))
#' Diags <- get_persistence_diagrams(point.sets = PS, maxdimension = 1, maxscale = 5, ncpu = 1)
#' @export

get_persistence_diagrams <- function(point.sets = NULL, 
                                     maxdimension = NULL, 
                                     maxscale = NULL, 
                                     use.dtm = FALSE,
                                     m0 = NULL,
                                     grid.by = NULL,
                                     ncpu = 1) {
  Diag <- list()
  i <- NULL
  # Set up cluster for parallelization
  cluster <- parallel::makeCluster(ncpu)
  doParallel::registerDoParallel(cluster)
  n <- length(point.sets)
  if(!use.dtm) {
    Diag <- foreach::foreach(i = c(1:n), .packages = 'TDA') %dopar% {
      TDA::ripsDiag(X = point.sets[[i]],
                    maxdimension,
                    maxscale,
                    library = "GUDHI",
                    printProgress = FALSE)[["diagram"]]
    }
  } else {
    point.sets <- lapply(point.sets, function(x) {scale(x, center = TRUE, scale = FALSE)})
    min.dims <- apply(sapply(point.sets, function(x) {apply(x, 2, min)}), 1, min)
    max.dims <- apply(sapply(point.sets, function(x) {apply(x, 2, max)}), 1, max)
    lims <- NULL
    for(l in 1:length(max.dims)) {
      lims <- cbind(lims, c(min.dims[l], max.dims[l]))
    }
    Diag <- foreach::foreach(i = c(1:n), .packages = 'TDA') %dopar% {
      TDA::gridDiag(X = point.sets[[i]], 
                    FUN = TDA::dtm, 
                    m0 = m0, 
                    lim = lims,
                    by = grid.by, 
                    sublevel = TRUE, 
                    library = "GUDHI", 
                    printProgress = FALSE )$diagram
    }
  }
  on.exit(parallel::stopCluster(cluster))
  return(Diag)
}

#' get_kernel_matrix
#' 
#' Compute kernel/distance matrix between persistence diagrams.
#' 
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#' @param Diag list of persistence diagrams as n x 3 matrices
#' @param method which kernel or distance to compute. One of sWd (for sliced Wasserstein kernel) or pssk (for the persistence scale-space kernel)
#' @param sigma kernel bandwidth
#' @param return.dist logical (default: FALSE) for method sWd, whether to return the sliced Wasserstein distance matrix instead of the kernel.
#' @param M number of slices for the sliced Wasserstein kernel
#' @param dimensions vector of the dimensions of the topological features to consider, if NULL (default) use all available dimensions
#' @param ncpu number of parallel threads to use for computation
#' @return a matrix
#' @examples 
#' PS <- list(data.frame(x = c(2.4,-6.9,4.6,-0.7,-3.3,-4.9,-3.5,-3.5,4.2,-7),
#'                       y = c(5.7,1.9,4.8,3.4,-3,-2.1,7.2,1.8,6.1,-1.6),
#'                       z = c(2.7,-0.1,-0.7,-0.6,0.4,-1.5,-0.6,-0.9,2.2,0.7)),
#'            data.frame(x = c(0,0,3.1,-5.6,-5,-7.4,-0.7,-7.7,-6.7,4,4.2,0.2,5.8,3.9,3.9),
#'                       y = c(6.3,-6.1,-3.5,4.6,-4.1,0.3,8.8,-2.3,2.9,3.7,-1.4,-3.9,5.5,-1.2,-6.7),
#'                       z = c(-1.5,1.7,-0.4,-1.4,1.8,1.7,-0.9,-1.8,-0.5,1.7,1.3,0.5,-1.4,1.6,-0.1)),
#'            data.frame(x = c(-9.8,-5.2,12.5,2.5,4.5,1.3,-0.2,0.4,9.3,-1.4,0.5,-1.1,-7.7),
#'                       y = c(-4.2,1.5,-0.5,12,-3,-7.2,10.9,6.7,-1.3,10,6.7,-6.2,2.9),
#'                       z = c(3.4,-3.8,-1.4,1.8,3.5,2.5,2.6,-4.8,-3.8,3.9,4.1,-3.6,-4)))
#' Dgs <- get_persistence_diagrams(point.sets = PS, maxdimension = 1, maxscale = 5, ncpu = 1)
#' K <- get_kernel_matrix(Diag = Dgs, method = 'sWd', dimensions = c(0,1), M = 10, sigma = 5)
#' @export

get_kernel_matrix <- function(Diag = NULL, method = c("sWd", "pssk"), dimensions = NULL, return.dist = FALSE, M = NULL, sigma = NULL, ncpu = 1) {
  n <- length(Diag)
  K <- matrix(NA, nrow = n, ncol = n)
  i <- j <- NULL
  cluster <- parallel::makeCluster(ncpu)
  doParallel::registerDoParallel(cluster)
  if(method == "sWd") {
    # Compute pairwise sliced Wasserstein distances
    K <- foreach::foreach(j = c(1:n), .combine='cbind') %:%
      foreach::foreach(i = c(1:n), .combine = 'c') %dopar% {
        Di <- Diag[[i]]
        Dj <- Diag[[j]]
        sliced_Wd(Di, Dj, M, sigma, dimensions, return.dist)
      }
  } else if (method == "pssk") {
    K <- foreach::foreach(j = c(1:n), .combine='cbind') %:%
      foreach::foreach(i = c(1:n), .combine = 'c') %dopar% {
        Di <- Diag[[i]]
        Dj <- Diag[[j]]
        pssk(Di, Dj, sigma, dimensions)
      }
  } else {
    stop(paste0("Unknown kernel method: ", method,".\n"))
  }
  on.exit(parallel::stopCluster(cluster))
  return(K)  
}
