#' The logit transformation
#'
#' This function takes the expected value mu (in the range of 0 to 1) and transforms
#' it to the linear predictor yita in the real number domain.
#'
#' @param x A numeric value for the expected value mu
#'
#' @return A numeric value for the linear predictor yita
#'
#' @examples
#' \dontrun{logit(0.5)}
#'
#' @export
logit <- function(x) {
  log(x / (1 - x))
}

#' Conversion of genoSnpMatrix to a genotype matrix
#'
#' This function takes a genoSnpMatrix in SnpMatrix format and converts it to a
#' genotype matrix with 0, 1, 2 encoding. By default, it standardizes the matrix
#' using the formula (x - 2*p) / sqrt(2*p*(1-p)), where x is the encoded numeric
#' value in each column and p is the minor allele frequency (maf) in each column.
#'
#' @param genoSnpMatrix a genotype matrix in SnpMatrix format
#' @param standardization logical, indicating whether to standardize the matrix
#'
#' @return a genotype matrix with 0, 1, 2 encoding and optionally standardized
#'
#' @export
#'
#' @examples
#' \dontrun{
#' load("./inst/PhenoGenoMap.RData")
#' G <- snp_Data$genotypes
#' Ga <- define_Ga(G)}
#'
define_Ga <- function(genoSnpMatrix, standardization = TRUE) {
  # Check if genoSnpMatrix is of SnpMatrix class
  if (!inherits(genoSnpMatrix, "SnpMatrix")) {
    stop("Input genoSnpMatrix must be of SnpMatrix class.")
  }

  # Convert genoSnpMatrix to a numeric matrix
  Ga <- as(genoSnpMatrix, "numeric")

  # Function to fill missing values with mode
  fill_na <- function(genotype_matrix) {
    for (col in 1:ncol(genotype_matrix)) {
      if (anyNA(genotype_matrix[, col])) {
        mode_value <- as.numeric(names(which.max(table(genotype_matrix[, col]))))
        genotype_matrix[is.na(genotype_matrix[, col]), col] <- mode_value
      }
    }
    return(genotype_matrix)
  }

  # Fill missing values in Ga
  Ga <- fill_na(Ga)

  # Function to flip 0s to 2s and 2s to 0s
  flip.matrix <- function(x) {
    zero2 <- which(x == 0)
    two0 <- which(x == 2)
    x[zero2] <- 2
    x[two0] <- 0
    return(x)
  }

  # Apply flip.matrix to Ga
  Ga <- flip.matrix(Ga)

  # Standardize Ga if standardization is TRUE
  if (standardization) {
    calculate_Xs <- function(Ga) {
      n <- nrow(Ga)
      m <- ncol(Ga)

      p <- apply(Ga, 2, function(marker) {
        A_freq <- sum(marker) / (2 * n)

        return(A_freq)
      })

      #Z <- matrix(0, nrow = n, ncol = m)

      # Generate Z matrix
      for (j in 1:m) {
        if (all(Ga[, j] == 0)) {
          next  # If the elements of a column are all 0, skip the processing of this column
        }
        pj <- p[j]
        Ga[, j] <- (Ga[, j] - 2 * pj)/ sqrt(2*(pj * (1 - pj)))
      }

      return(Ga)
    }
    Ga <- calculate_Xs(Ga)
  }

  return(Ga)
}

#' 1D scanning for additive and additive-by-environment SNPs
#'
#' This function conducts a one-dimensional scan on SNPs based on the given parameters.
#'
#' @param phenotype a vector of phenotype values
#' @param covariates a matrix of covariates (PCA)
#' @param environment a vector of environmental factors
#' @param weights a vector of weights
#' @param data a data frame containing the input data
#' @param threshold the threshold for precise calculation and estimation
#' @param theta the theta parameter when family is "negative.binomial"
#' @param genodata a genotype matrix
#' @param map the map file
#' @param family the distribution family
#' @param parallel logical, parallel computation
#' @param num_threads the number of threads for parallel computation
#' @param filename the filename for output
#'
#' @return A list containing the 1D scan results, including chromosome, SNP name,
#'         position, score value and p-value.
#' @export
#'
#' @examples
#' \dontrun{rae <- single.snp.score(phenotype = y,covariates = cov, environment = env,
#'                          genodata = Ga,map = map,family = "binomial")}
#'
single.snp.score <- function(phenotype = phenotype, covariates = NULL, environment = NULL, weights = NULL,
                             data = NULL,threshold = 0.01, theta = NULL, genodata = genodata,map = map,
                             family = "gaussian", parallel = FALSE,num_threads = 2,filename = NULL) {
  weight <- numeric(2)
  # Process phenotype, covariates, and environment inputs
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("The 'data' input must be a data frame.")
    }
    if(!phenotype %in% colnames(data)){
      stop("The 'phenotype' input must be a single column name present in the 'data' data frame.")
    }
    if (is.character(phenotype)) {
      phenotype <- data[[phenotype]]
    }
    if (!all(!is.na(phenotype) & phenotype != "" & phenotype != ".")) {
      stop("Phenotype contains NA, empty values, or dots.")
    }
    if (!is.null(covariates)) {
      if (is.character(covariates)) {
        covariates_names <- strsplit(covariates, "\\+")[[1]]
        covariates <- as.matrix(data[, covariates_names, drop = FALSE])
      }
      if (!all(!is.na(covariates) & covariates != "" & covariates != ".")) {
        stop("covariates contains NA, empty values, or dots.")
      }
      cov <- as.matrix(cbind(1, covariates))
    } else {
      cov <- matrix(1, nrow = length(phenotype), ncol = 1)
    }
    if (!is.null(environment)) {
      if(!environment %in% colnames(data)){
        stop("The 'environment' input must be a single column name present in the 'data' data frame.")
      }
      environment <- data[[environment]]
      if (!all(!is.na(environment) & environment != "" & environment != ".")) {
        stop("Environment contains NA, empty values, or dots.")
      }
      envir <- length(unique(environment))
    }else{
      envir <- 1
    }
    if (!is.null(weights)) {
      if(!weights %in% colnames(data)){
        stop("The 'weights' input must be a single column name present in the 'data' data frame.")
      }
      weight <- data[[weights]]
      if (!all(!is.na(weight) & weight != "" & weight != ".")) {
        stop("Weights contains NA, empty values, or dots.")
      }
    }
  } else {
    if (!(is.vector(phenotype) || (is.matrix(phenotype) && ncol(phenotype) == 1))) {
      stop("Error: phenotype must be a vector or a matrix with ncol equal to 1.")
    }

    phenotype <- phenotype

    if (!all(!is.na(phenotype) & phenotype != "" & phenotype != ".")) {
      stop("Phenotype contains NA, empty values, or dots.")
    }
    if (!is.null(covariates)) {
      if (!is.matrix(covariates)) {
        covariates <- as.matrix(covariates)
      }
      if (nrow(covariates) != length(phenotype)) {
        print(nrow(covariates))
        stop("Covariates length does not match the length of phenotype.")
      }

      if (!all(!is.na(covariates) & covariates != "" & covariates != ".")) {
        stop("Covariates contains NA, empty values, or dots.")
      }
      cov <- cbind(1, covariates)
    } else {
      cov <- matrix(1, nrow = length(phenotype), ncol = 1)
    }
    if (!is.null(environment)) {
      if (!(is.vector(environment) || (is.matrix(environment) && ncol(environment) == 1))) {
        stop("Error: environment must be a vector or a matrix with ncol equal to 1.")
      }
      environment <- environment
      if (length(environment) != length(phenotype)) {
        stop("Environment length does not match the length of phenotype.")
      }

      if (!all(!is.na(environment) & environment != "" & environment != ".")) {
        stop("Environment contains NA, empty values, or dots.")
      }
      envir <- length(unique(environment))
    }else{
      envir <- 1
    }
    if (!is.null(weights)) {
      if (!(is.vector(weights) || (is.matrix(weights) && ncol(weights) == 1))) {
        stop("Error: weights must be a vector or a matrix with ncol equal to 1.")
      }
      weight <- weights
      if (length(weight) != length(phenotype)) {
        stop("Weights length does not match the length of phenotype.")
      }

      if (!all(!is.na(weight) & weight != "" & weight != ".")) {
        stop("Weights contains NA, empty values, or dots.")
      }
    }
  }
  if ( is(genodata, "SnpMatrix")) {
    stop("The 'genodata' input cannot be in SnpMatrix format.")
  }
  # Check genodata dimensions
  if (nrow(genodata) != (length(phenotype)/envir)) {
    stop("The number of rows in genodata must match the length of phenotype.")
  }

  if (ncol(genodata) != nrow(map)) {
    stop("The number of columns in genodata must match the number of rows in map.")
  }

  # Check map column names
  required_colnames <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")
  if (!all(colnames(map) == required_colnames)) {
    stop("The column names of map must be: 'chromosome', 'snp.name', 'cM', 'position', 'allele.1', 'allele.2'. Please rename the columns accordingly.")
  }

  # Check if the 'snp.name' column in 'map' matches the column names of 'genodata' in both length and order
  if (!identical(map$snp.name, colnames(genodata))) {
    stop("The 'snp.name' column in 'map' does not match the column names of 'genodata' in terms of order.")
  }

  all_families <- c("gaussian", "binomial", "poisson", "negative.binomial")
  if (!family %in% all_families) {
    stop("Invalid family parameter. Supported families are 'gaussian', 'poisson',
         'binomial' and 'negative.binomial'.")
  }

  if (family == "gaussian") {
    phenotype <- scale(phenotype)
  }

  if (family == "binomial") {
    if (any(phenotype > 1)) {
      if (is.null(weights)) {
        stop("Phenotype data for Binomial family must have weights specified for counts greater than 1.")
      }else{
        phenotype <- phenotype / weight
      }
    } else if (any(phenotype < 0) || any(phenotype > 1)) {
      stop("Phenotype data for Binomial family must only contain values between 0 and 1.")
    }

    # If weights are NULL, generate weight vector
    if (is.null(weights) && any(phenotype >= 0 & phenotype <= 1)) {
      if (!all(phenotype %in% c(0, 1))) {
        weight <- rep(100, length(phenotype))
        weights <- weight
      }
    }

  }

  if (family == "poisson" || family == "negative.binomial") {
    if (any(phenotype < 0)) {
      stop("Phenotype data for '", family, "' family must be non-negative.")
    }
  }
  if (family == "negative.binomial") {
    if (is.null(theta)) {
      ml <- glm.nb(phenotype ~ covariates)
      theta <- ml$theta
    }else{
      theta <- theta
    }

  }else{
    theta <- 0.0
  }

  ra <- glm_NR(phenotype,weight, cov, family,theta)
  ramu <- Ca_Mu(phenotype,weight,  cov, family, ra)
  rawii <- Ca_Wii(phenotype,weight,  family, ramu,ra,theta)
  #Calculate ra1 using score_test function
  if(!is.null(covariates)){
    tXcW <-  t(covariates) * rawii
    if (qr(tXcW %*% covariates)$rank == ncol(tXcW %*% covariates)) {
      XtWXc_inv <- solve(tXcW %*% covariates)
    } else {
      XtWXc_inv <- ginv(tXcW %*% covariates)
    }
    proj <- covariates %*% XtWXc_inv %*% tXcW
  }else{
    proj <-matrix(0,ncol(cov),ncol(cov))
  }
  RcppParallel::setThreadOptions(numThreads = num_threads)
  ra1 <- score_test(phenotype, ramu, rawii, genodata,proj,ra,family, envir, threshold,theta,epistasis_test = FALSE, parallel = parallel)

  #Create matrices and store p_values and x2 values
  snp_name <- as.matrix(colnames(genodata))
  chromosome <- as.matrix(map[,1])
  position <- as.matrix(map[,4])
  x2_values <- as.matrix(ra1$x2)
  p_values <- as.matrix(ra1$p_value)

  gwasp1 <- cbind(snp_name, chromosome, position, x2_values, p_values)
  colnames(gwasp1) <- c("SNP", "CHR", "BP", "X2", "P")

  # If a file name is specified, output to the file
  if (!is.null(filename)) {
    write.table(gwasp1, file = filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  # Set row names of matrices as column names of genodata
  rownames(p_values) <- colnames(genodata)
  rownames(x2_values) <- colnames(genodata)
  colnames(p_values) <- paste0("P.", envir, "df")
  colnames(x2_values) <- paste0("Chi.squared.", envir, "df")
  # Return the matrices
  if (family == "negative.binomial") {
    return(list(gwas1D = gwasp1 ,p_values = p_values, x2_values = x2_values,theta = theta))
  }else{
    return(list(gwas1D = gwasp1 ,p_values = p_values, x2_values = x2_values))
  }

}

#' 2D scanning for additive-additive epistasis and epistasis-by-environment SNPs
#'
#' This function conducts a two-dimensional scan for on a SNP matrix based on the
#' given parameters and significant SNPs from a 1D scanning.
#'
#'
#' @param significantASNP a vector of significant SNP names from a 1D scanning
#' @param phenotype a vector of phenotype values
#' @param covariates a matrix of covariates (PCA)
#' @param environment a vector of environmental factors
#' @param weights a vector of weights
#' @param data a data frame containing the input data
#' @param threshold the threshold for precise calculation and estimation
#' @param theta the theta parameter when family is "negative.binomial"
#' @param genodata a genotype matrix
#' @param map the map file
#' @param family the distribution family
#' @param parallel logical, parallel computation
#' @param num_threads the number of threads for parallel computation
#' @param filename the filename for output
#'
#' @return A list containing the 1D scan results for additive-additive epistasis,
#'        including cnromosome, SNP name, position,p-values, and score values.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sig <- c("SNP25","SNP67")
#' raae <- epi.snp.score(significantASNP = sig, phenotype = y, covariates = cov,
#'                      environment = env,genodata = Ga,map = map,family = "binomial",
#'                      parallel = T, num_threads = 4,filename = "gwasp2.txt")}
#'
epi.snp.score <- function(significantASNP = NULL, phenotype = phenotype, covariates = NULL, environment = NULL,
                          weights = NULL, data = NULL,threshold = 0.01,theta = NULL,genodata = genodata,
                          map = map,family = "gaussian", parallel = FALSE, num_threads = 2, filename = NULL) {

  if (!is.null(significantASNP)) {
    if (!(is.vector(significantASNP) || (is.matrix(significantASNP) && ncol(significantASNP) == 1))) {
      stop("Error: If significantASNP is provided, it must be a vector or a matrix with ncol equal to 1.")
    }
  }

  weight <- numeric(2)
  # Process phenotype, covariates, and environment inputs
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("The 'data' input must be a data frame.")
    }
    if(!phenotype %in% colnames(data)){
      stop("The 'phenotype' input must be a single column name present in the 'data' data frame.")
    }
    if (is.character(phenotype)) {
      phenotype <- data[[phenotype]]
    }
    if (!all(!is.na(phenotype) & phenotype != "" & phenotype != ".")) {
      stop("Phenotype contains NA, empty values, or dots.")
    }
    if (!is.null(covariates)) {
      if (is.character(covariates)) {
        covariates_names <- strsplit(covariates, "\\+")[[1]]
        covariates <- as.matrix(data[, covariates_names, drop = FALSE])
      }
      if (!all(!is.na(covariates) & covariates != "" & covariates != ".")) {
        stop("covariates contains NA, empty values, or dots.")
      }
      cov <- as.matrix(cbind(1, covariates))
    } else {
      cov <- matrix(1, nrow = length(phenotype), ncol = 1)
    }
    if (!is.null(environment)) {
      if(!environment %in% colnames(data)){
        stop("The 'environment' input must be a single column name present in the 'data' data frame.")
      }
      environment <- data[[environment]]
      if (!all(!is.na(environment) & environment != "" & environment != ".")) {
        stop("Environment contains NA, empty values, or dots.")
      }
      envir <- length(unique(environment))
    }else{
      envir <- 1
    }
    if (!is.null(weights)) {
      if(!weights %in% colnames(data)){
        stop("The 'weights' input must be a single column name present in the 'data' data frame.")
      }
      weight <- data[[weights]]
      if (!all(!is.na(weight) & weight != "" & weight != ".")) {
        stop("Weights contains NA, empty values, or dots.")
      }
    }
  } else {
    if (!(is.vector(phenotype) || (is.matrix(phenotype) && ncol(phenotype) == 1))) {
      stop("Error: phenotype must be a vector or a matrix with ncol equal to 1.")
    }

    phenotype <- phenotype

    if (!all(!is.na(phenotype) & phenotype != "" & phenotype != ".")) {
      stop("Phenotype contains NA, empty values, or dots.")
    }
    if (!is.null(covariates)) {
      if (!is.matrix(covariates)) {
        covariates <- as.matrix(covariates)
      }
      if (nrow(covariates) != length(phenotype)) {
        print(nrow(covariates))
        stop("Covariates length does not match the length of phenotype.")
      }

      if (!all(!is.na(covariates) & covariates != "" & covariates != ".")) {
        stop("Covariates contains NA, empty values, or dots.")
      }
      cov <- cbind(1, covariates)
    } else {
      cov <- matrix(1, nrow = length(phenotype), ncol = 1)
    }
    if (!is.null(environment)) {
      if (!(is.vector(environment) || (is.matrix(environment) && ncol(environment) == 1))) {
        stop("Error: environment must be a vector or a matrix with ncol equal to 1.")
      }
      environment <- environment
      if (length(environment) != length(phenotype)) {
        stop("Environment length does not match the length of phenotype.")
      }

      if (!all(!is.na(environment) & environment != "" & environment != ".")) {
        stop("Environment contains NA, empty values, or dots.")
      }
      envir <- length(unique(environment))
    }else{
      envir <- 1
    }
    if (!is.null(weights)) {
      if (!(is.vector(weights) || (is.matrix(weights) && ncol(weights) == 1))) {
        stop("Error: weights must be a vector or a matrix with ncol equal to 1.")
      }
      weight <- weights
      if (length(weight) != length(phenotype)) {
        stop("Weights length does not match the length of phenotype.")
      }

      if (!all(!is.na(weight) & weight != "" & weight != ".")) {
        stop("Weights contains NA, empty values, or dots.")
      }
    }
  }
  if ( is(genodata, "SnpMatrix")) {
    stop("The 'genodata' input cannot be in SnpMatrix format.")
  }

  # Check genodata dimensions
  if (nrow(genodata) != (length(phenotype)/envir)) {
    stop("The number of rows in genodata must match the length of phenotype.")
  }

  if (ncol(genodata) != nrow(map)) {
    stop("The number of columns in genodata must match the number of rows in map.")
  }

  # Check map column names
  required_colnames <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")
  if (!all(colnames(map) == required_colnames)) {
    stop("The column names of map must be: 'chromosome', 'snp.name', 'cM', 'position', 'allele.1', 'allele.2'. Please rename the columns accordingly.")
  }

  # Check if the 'snp.name' column in 'map' matches the column names of 'genodata' in both length and order
  if (!identical(map$snp.name, colnames(genodata))) {
    stop("The 'snp.name' column in 'map' does not match the column names of 'genodata' in terms of order.")
  }

  all_families <- c("gaussian", "binomial", "poisson", "negative.binomial")
  if (!family %in% all_families) {
    stop("Invalid family parameter. Supported families are 'gaussian', 'poisson',
         'binomial' and 'negative.binomial'.")
  }

  if (family == "gaussian") {
    phenotype <- scale(phenotype)
  }
  # if (family == "binomial") {
  #   if (is.null(weights)) {
  #     if (any(phenotype < 0) || any(phenotype > 1)) {
  #       stop("Phenotype data for Binomial family must only contain values between 0 and 1.")
  #     }
  #   } else {
  #     if (any(phenotype < 0)) {
  #       stop("Phenotype data for 'binomial' family must be non-negative.")
  #     }
  #   }
  # }
  if (family == "binomial") {
    if (any(phenotype > 1)) {
      if (is.null(weights)) {
        stop("Phenotype data for Binomial family must have weights specified for counts greater than 1.")
      }else{
        phenotype <- phenotype / weight
      }
    } else if (any(phenotype < 0) || any(phenotype > 1)) {
      stop("Phenotype data for Binomial family must only contain values between 0 and 1.")
    }

    # If weights are NULL, generate weight vector
    if (is.null(weights) && any(phenotype >= 0 & phenotype <= 1)) {
      if (!all(phenotype %in% c(0, 1))) {
        weight <- rep(100, length(phenotype))
        weights <- weight
      }
    }
  }
  if (family == "poisson" || family == "negative.binomial") {
    if (any(phenotype < 0)) {
      stop("Phenotype data for '", family, "' family must be non-negative.")
    }
  }

  if (!is.null(significantASNP)) {
    Xc <- matrix(nrow = nrow(genodata), ncol = length(significantASNP))

    for (i in 1:length(significantASNP)) {
      cn <- significantASNP[i]
      ci <- which(colnames(genodata) == cn)
      Xc[, i] <- genodata[, ci]
    }
    Xcc <- matrix(nrow = (nrow(genodata)*envir), ncol = (length(significantASNP)*envir))
    Xcc <- generateSubmatrix(Xc, envir)
    cov <- cbind(cov,Xcc)
    if(!is.null(covariates)){
      covariates <- cbind(covariates,Xcc)
    }
  }
  if (family == "negative.binomial") {
    if (is.null(theta)) {
      ml <- glm.nb(phenotype ~ cov)
      theta <- ml$theta
    }else{
      theta <- theta
    }

  }else{
    theta <- 0.0
  }

  ra <- glm_NR(phenotype,weight, cov, family,theta)
  ramu <- Ca_Mu(phenotype,weight, cov, family, ra)
  rawii <- Ca_Wii(phenotype,weight, family, ramu,ra,theta)

  if(!is.null(covariates)){

    tXcW <-  t(covariates) * rawii
    if (qr(tXcW %*% covariates)$rank == ncol(tXcW %*% covariates)) {
      XtWXc_inv <- solve(tXcW %*% covariates)
    } else {
      XtWXc_inv <- ginv(tXcW %*% covariates)
    }
    proj <- covariates %*% XtWXc_inv %*% tXcW
  }else{
    proj <- matrix(0,ncol(cov),ncol(cov))
  }
  #Calculate ra, ramu, rawii using provided functions
  #Calculate ra1 using score_test function
  RcppParallel::setThreadOptions(numThreads = num_threads)
  ra1 <- score_test(phenotype, ramu, rawii, genodata,proj, ra,family, envir,threshold,theta, epistasis_test = TRUE, parallel = parallel)

  #Create matrices and store p_values and x2 values
  x2_values <- as.matrix(ra1$x2)
  p_values <- as.matrix(ra1$p_value)
  gwasp2 <- matrix(NA,nrow = length(x2_values),ncol = 8)

  colna <- colnames(genodata)
  m <- ncol(genodata)
  #SNP <- character(m*(m-1)/2)

  col_idx <- 1
  for (i in 1:(m-1)) {
    for (j in (i+1):m) {
      gwasp2[col_idx,1] <- colna[i]
      gwasp2[col_idx,2] <- map[i,1]
      gwasp2[col_idx,3] <- map[i,4]
      gwasp2[col_idx,4] <- colna[j]
      gwasp2[col_idx,5] <- map[j,1]
      gwasp2[col_idx,6] <- map[j,4]
      #SNP[col_idx] <- paste0(colna[i],"_",colna[j])
      col_idx <- col_idx + 1
    }
  }
  gwasp2[,7] <- x2_values
  gwasp2[,8] <- p_values
  colnames(gwasp2) <- c("SNP1", "CHR1", "BP1","SNP2", "CHR2", "BP2", "X2", "P")

  if (!is.null(filename)) {
    write.table(gwasp2, file = filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  }

  # Set row names of matrices as column names of genodata
  #rownames(p_values) <- SNP
  #rownames(x2_values) <- SNP
  colnames(p_values) <- paste0("P.", envir, "df")
  colnames(x2_values) <- paste0("Chi.squared.", envir, "df")
  # Return the matrices
  if (family == "negative.binomial") {
    return(list(gwas2D = gwasp2, p_values = p_values, x2_values = x2_values,theta = theta))
  }else{
    return(list(gwas2D = gwasp2, p_values = p_values, x2_values = x2_values))
  }

}

#' Filtering and parameter estimation for significant SNPs in 1D Scaning
#'
#' This function filters the results from an additive 1D scaning to obtain significant
#' SNP markers and allows for parameter estimation.
#'
#' @param phenotype a vector of phenotype values
#' @param covariates a matrix of covariates
#' @param environment a vector of environmental
#' @param weights a vector of weights
#' @param data a data frame containing the input data
#' @param theta the theta parameter
#' @param genodata a genotype matrix
#' @param map the map file
#' @param Pvalues a vector of 1D p-values
#' @param Bonferr bonferroni correction factor, 0.5/M (M is the number of SNPs)
#' @param family the distribution family
#' @param lasso logical, LASSO method
#' @param selection logical, forward selection method
#' @param estimate logical, estimating parameters for significant markers
#' @param GEI_random logical, gene-environment interaction effects as random effects
#'
#' @return A list containing information about 1D significant SNPs, covariates coefficients,
#'         genotype coefficients, GEI coefficients, GEI VarCorr, and model AIC.
#'
#' @export
#'
#' @examples
#'  \dontrun{rae1 <- single.snp.estimate(phenotype = y, covariates = cov, environment = env,
#'                             genodata = Ga,map = map,Pvalues = p1, Bonferr = 0.5,
#'                             family = "binomial",selection = T,estimate = T)}
#'
single.snp.estimate <- function(phenotype = phenotype, covariates = NULL, environment = NULL, weights = NULL,data = NULL,
                                theta = NULL, genodata = genodata,map = map,Pvalues = Pvalues, Bonferr = 0.5,
                                family = "gaussian",lasso = FALSE, selection = FALSE,estimate = FALSE, GEI_random = TRUE) {

  weight <- numeric(2)
  # Process phenotype, covariates, and environment inputs
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("The 'data' input must be a data frame.")
    }
    if(!phenotype %in% colnames(data)){
      stop("The 'phenotype' input must be a single column name present in the 'data' data frame.")
    }
    if (is.character(phenotype)) {
      phenotype <- data[[phenotype]]
    }
    if (!all(!is.na(phenotype) & phenotype != "" & phenotype != ".")) {
      stop("Phenotype contains NA, empty values, or dots.")
    }
    if (!is.null(covariates)) {
      if (is.character(covariates)) {
        covariates_names <- strsplit(covariates, "\\+")[[1]]
        covariates <- as.matrix(data[, covariates_names, drop = FALSE])
      }
      if (!all(!is.na(covariates) & covariates != "" & covariates != ".")) {
        stop("Covariates contains NA, empty values, or dots.")
      }
      cov <- as.matrix(cbind(1, covariates))
    } else {
      cov <- matrix(1, nrow = length(phenotype), ncol = 1)
    }

    if (!is.null(environment)) {
      if(!environment %in% colnames(data)){
        stop("The 'environment' input must be a single column name present in the 'data' data frame.")
      }
      environment <- data[[environment]]
      if (!all(!is.na(environment) & environment != "" & environment != ".")) {
        stop("Environment contains NA, empty values, or dots.")
      }
      envir <- length(unique(environment))
      environment <- factor(environment)
      ename <-  levels(environment)
    }else{
      envir <- 1
    }
    if (!is.null(weights)) {
      if(!weights %in% colnames(data)){
        stop("The 'weights' input must be a single column name present in the 'data' data frame.")
      }
      weight <- data[[weights]]
      if (!all(!is.na(weight) & weight != "" & weight != ".")) {
        stop("Weights contains NA, empty values, or dots.")
      }
    }
  } else {
    if (!(is.vector(phenotype) || (is.matrix(phenotype) && ncol(phenotype) == 1))) {
      stop("Error: phenotype must be a vector or a matrix with ncol equal to 1.")
    }

    phenotype <- phenotype
    if (!all(!is.na(phenotype) & phenotype != "" & phenotype != ".")) {
      stop("Phenotype contains NA, empty values, or dots.")
    }
    if (!is.null(covariates)) {
      if (!is.matrix(covariates)) {
        covariates <- as.matrix(covariates)
      }
      if (nrow(covariates) != length(phenotype)) {
        stop("Covariates length does not match the length of phenotype.")
      }

      if (!all(!is.na(covariates) & covariates != "" & covariates != ".")) {
        stop("Covariates contains NA, empty values, or dots.")
      }
      cov <- cbind(1, covariates)
    } else {
      cov <- matrix(1, nrow = length(phenotype), ncol = 1)
    }
    if (!is.null(environment)) {
      if (!(is.vector(environment) || (is.matrix(environment) && ncol(environment) == 1))) {
        stop("Error: environment must be a vector or a matrix with ncol equal to 1.")
      }
      environment <- environment
      if (length(environment) != length(phenotype)) {
        stop("Environment length does not match the length of phenotype.")
      }

      if (!all(!is.na(environment) & environment != "" & environment != ".")) {
        stop("Environment contains NA, empty values, or dots.")
      }
      environment <- factor(environment)
      envir <- length(unique(environment))
      ename <-  levels(environment)
    }else{
      envir <- 1
    }
    if (!is.null(weights)) {
      if (!(is.vector(weights) || (is.matrix(weights) && ncol(weights) == 1))) {
        stop("Error: weights must be a vector or a matrix with ncol equal to 1.")
      }
      weight <- weights
      if (length(weight) != length(phenotype)) {
        stop("Weights length does not match the length of phenotype.")
      }

      if (!all(!is.na(weight) & weight != "" & weight != ".")) {
        stop("Weights contains NA, empty values, or dots.")
      }
    }
  }

  if ( is(genodata, "SnpMatrix")) {
    stop("The 'genodata' input cannot be in SnpMatrix format.")
  }

  all_families <- c("gaussian", "binomial", "poisson", "negative.binomial")
  if (!family %in% all_families) {
    stop("Invalid family parameter. Supported families are 'gaussian', 'poisson',
         'binomial' and 'negative.binomial'.")
  }

  # Check genodata dimensions
  if (nrow(genodata) != (length(phenotype)/envir)) {
    stop("The number of rows in genodata must match the length of phenotype.")
  }

  if (ncol(genodata) != nrow(map)) {
    stop("The number of columns in genodata must match the number of rows in map.")
  }

  # Check map column names
  required_colnames <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")
  if (!all(colnames(map) == required_colnames)) {
    stop("The column names of map must be: 'chromosome', 'snp.name', 'cM', 'position', 'allele.1', 'allele.2'. Please rename the columns accordingly.")
  }

  # Check if the 'snp.name' column in 'map' matches the column names of 'genodata' in both length and order
  if (!identical(map$snp.name, colnames(genodata))) {
    stop("The 'snp.name' column in 'map' does not match the column names of 'genodata' in terms of order.")
  }

  phenotype1 <- phenotype
  if (family == "gaussian") {
    phenotype1 <- scale(phenotype)
  }
  if (family == "binomial") {
    if (any(phenotype > 1)) {
      if (is.null(weights)) {
        stop("Phenotype data for Binomial family must have weights specified for counts greater than 1.")
      }else{
        phenotype <- phenotype / weight
        phenotype1 <- phenotype
      }
    } else if (any(phenotype < 0) || any(phenotype > 1)) {
      stop("Phenotype data for Binomial family must only contain values between 0 and 1.")
    }

    # If weights are NULL, generate weight vector
    if (is.null(weights) && any(phenotype >= 0 & phenotype <= 1)) {
      if (!all(phenotype %in% c(0, 1))) {
        weight <- rep(100, length(phenotype))
        weights <- weight
      }
    }
  }
  if (family == "poisson" || family == "negative.binomial") {
    if (any(phenotype < 0)) {
      stop("Phenotype data for '", family, "' family must be non-negative.")
    }
  }
  if (family == "negative.binomial") {
    if (is.null(theta)) {

      ml <- glm.nb(phenotype ~ cov)
      theta <- ml$theta

    }else{
      theta <- theta
    }

  }else{
    theta <- 0.0
  }
  #Pvalues <- Pvalues
  m <- ncol(genodata)
  #Pvalues <- p.adjust(Pvalues, method = "fdr")
  thr <- 0.5 / m
  thr1 <- Bonferr / m
  # Record indexes smaller than the FDR threshold
  significant_l_indices <- which(Pvalues <= thr1)
  if (length(significant_l_indices) == 0) {
    miss <- NULL
    message("No significant markers found.")
    return(miss)
  }

  if(lasso){
    if (family == "negative.binomial" ) {
      message("Lasso method is not applicable for negative binomial models. Performing Bonferroni correction instead.")
      fin <- significant_l_indices
    }else{
      if(length(significant_l_indices) >= 2 ){
        gs <- as.matrix(genodata[, significant_l_indices])
        if(envir == 1){

          colnames(gs)  <- significant_l_indices

        }else{

          gs <- generateSubmatrix(gs,envir)
          gro <- rep(significant_l_indices,envir)
          colnames(gs)  <- gro

        }
        if (family == "binomial" && !is.null(weights) ){
          phe1 <- phenotype * weight
          phe2 <- weight - phe1
          pheno <- cbind(phe1,phe2)
          cvobs <-  cv.glmnet(data.matrix(gs), as.matrix(pheno),family = family)
          obs <-  glmnet(data.matrix(gs), as.matrix(pheno),family = family,lambda = cvobs$lambda.min )
        }else{
          cvobs <-  cv.glmnet(data.matrix(gs), as.matrix(phenotype), family = family)
          obs <-  glmnet(data.matrix(gs), as.matrix(phenotype),  family = family,lambda = cvobs$lambda.min )
        }

        fin <- as.numeric(rownames(obs$beta)[as.numeric(obs$beta)>0])
        fin <- unique(fin)
      }else{
        fin <- significant_l_indices
      }
    }

  }

  if (selection) {
    # Generate threshold
    thr2 <- (Bonferr * 0.01) / m
    significant_s_indices <- which(Pvalues <= thr2)
    if (length(significant_s_indices) == 0) {
      thr2 <- (Bonferr * 0.1) / m
      significant_s_indices <- which(Pvalues <= thr2)
      if (length(significant_s_indices) == 0) {
        fin <- significant_l_indices
      }else if (identical(significant_s_indices, significant_l_indices)) {
        fin <- significant_l_indices
      } else {
        significant_l_indices <- significant_l_indices - 1
        significant_s_indices <- significant_s_indices - 1

        # Generate index vectors by forward selection method
        fin <- forward_selection(phenotype1,weight, cov, genodata, family, envir, thr,theta, significant_s_indices, significant_l_indices)
        fin <- fin + 1
      }
    }else if (identical(significant_s_indices, significant_l_indices)) {
      fin <- significant_l_indices
    } else {
      significant_l_indices <- significant_l_indices - 1
      significant_s_indices <- significant_s_indices - 1

      fin <- forward_selection(phenotype1, weight,cov, genodata, family, envir,thr, theta, significant_s_indices, significant_l_indices)
      fin <- fin + 1
    }

  }
  if( !selection){
    fin <- significant_l_indices
  }
  colsnp <- colnames(genodata)[fin]
  gwassig1 <- matrix(NA,nrow =length(fin) ,ncol = 4)
  gwassig1[,1] <- colsnp
  gwassig1[,2] <- map[fin,1]
  gwassig1[,3] <- map[fin,4]
  gwassig1[,4] <- Pvalues[fin]
  colnames( gwassig1) <- c("SNP", "CHR", "BP","P")
  if(estimate){
    n <- nrow(genodata)
    allga <- matrix(NA, nrow = n * envir, ncol = length(fin))
    for (i in 1:envir) {
      allga[((i - 1) * n + 1):(i * n), ] <- genodata[, fin]
    }

    coly <- c("y")
    # Create data frame
    if (is.null(covariates) && is.null(environment)) {
      quan <- as.data.frame(cbind(phenotype, allga))
      colnames(quan) <- c(coly, colsnp)
    } else if (is.null(covariates)) {
      quan <- as.data.frame(cbind(phenotype, environment, allga))
      cole <- c("environment")
      colnames(quan) <- c(coly,cole, colsnp)
    } else if (is.null(environment)) {
      colx <- colnames(covariates)

      quan <- as.data.frame(cbind(phenotype, covariates, allga))
      colnames(quan) <- c(coly,colx, colsnp)
    } else {
      colx <- colnames(covariates)
      cole <- c("environment")
      quan <- as.data.frame(cbind(phenotype, covariates, environment, allga))
      colnames(quan) <- c(coly, colx, cole, colsnp)
    }
    if (!is.null(weights)){
      weight <- as.matrix(weight)
      colw <- c("frequency")
      colnames(weight) <- colw
      quan <- cbind(quan,weight)
    }
    colx <- if (!is.null(covariates)) colnames(covariates) else 1

    # When envir is 1, treat covariates (if it exists) and allga as fixed effects, and use the glm function to calculate
    if (envir == 1) {
      formula <- as.formula(paste(coly, "~",paste(colx, collapse = "+"), "+", paste(colsnp, collapse = "+")))

      if (family == "binomial" && !is.null(weights) ){
        # formula <- as.formula(paste("cbind(",coly,",",colw,"-",coly,")", "~",paste(colx, collapse = "+"), "+", paste(colsnp, collapse = "+")))
        result <- glm(formula, data = quan,weights = frequency, family = family)
      }else if (family == "negative.binomial") {
        result <- glm(formula, data = quan, family = negative.binomial(theta))
      }else{
        result <- glm(formula, data = quan, family = family)
      }
      estimates <- coef(summary(result))
      cov_coe <- estimates[!rownames(estimates) %in% colsnp,]
      if (nrow(cov_coe) < length(colx)) {
        message("There are NAs in the estimated coefficients of covariates. The corresponding variables have been removed from the results.")
      }

      fe_coe <- matrix(NA,nrow = length(colsnp),ncol = ncol(estimates))
      rownames(fe_coe) <- colsnp
      colnames(fe_coe) <- colnames(estimates)
      ma1 <- colsnp %in% rownames(estimates)
      ma2 <- rownames(estimates) %in% colsnp
      fe_coe[ma1,] <- estimates[ma2,]
      #fe_coe <- estimates[ncol(cov)+1:ncol(estimates),]
      # std_errors <- sqrt(diag(vcov(result)))
      # z_values <- estimates / std_errors
      # coefficients <- data.frame(estimate = estimates, std.error = std_errors, z_value = z_values)
      aic <- AIC(result)
      return(list(gwas1D = gwassig1,covariates_coefficients = cov_coe ,genotype_coefficients = fe_coe, AIC = aic))
    } else {
      if (!GEI_random) {

        gf <- genodata[, fin]
        gf <- generateSubmatrix(gf, envir)

        rann <- c()
        for (e in ename) {
          for (snp in colsnp) {
            combination <- paste0(snp, "_", e)
            rann <- c(rann, combination)
          }
        }
        colnames(gf) <- rann

        # Add the gf matrix to the quan data frame, treat covariates (if it exists) and allga as fixed effects, and calculate it with the glm function
        quan <- cbind(quan, gf)
        formula <- as.formula(paste(coly, "~", paste(colx, collapse = "+"), "+", paste(colsnp, collapse = "+"),
                                    "+",paste(rann, collapse = "+")))
        if (family == "binomial" && !is.null(weights) ){
          result <- glm(formula, data = quan,weights = frequency, family = family)
        }else if (family == "negative.binomial") {
          result <- glm(formula, data = quan, family = negative.binomial(theta))
        }else{
          result <- glm(formula, data = quan, family = family)
        }


        fixed_effects <- coef(summary(result))
        cov_coe <- fixed_effects[!rownames(fixed_effects) %in% colsnp,]
        if (nrow(cov_coe) < length(colx)) {
          message("There are NAs in the estimated coefficients of covariates. The corresponding variables have been removed from the results.")
        }

        fe_coe <- matrix(NA,nrow = length(colsnp),ncol = ncol(fixed_effects))
        rownames(fe_coe) <- colsnp
        colnames(fe_coe) <- colnames(fixed_effects)
        ma1 <- colsnp %in% rownames(fixed_effects)
        ma2 <- rownames(fixed_effects) %in% colsnp
        fe_coe[ma1,] <- fixed_effects[ma2,]

        gei_coe <- matrix(NA,nrow = length(rann),ncol = ncol(fixed_effects))
        colnames(gei_coe) <- colnames(fixed_effects)
        rownames(gei_coe) <- rann
        ma3 <- rann %in% rownames(fixed_effects)
        ma4 <- rownames(fixed_effects) %in% rann
        gei_coe[ma3,] <- fixed_effects[ma4,]
        # fixed_effects <- coef(summary(result))
        # cov_coe <- fixed_effects[1:ncol(cov),]
        # fe_coe <- fixed_effects[colsnp,]
        # gei_coe <- matrix(NA,nrow = length(rann),ncol = ncol(fixed_effects))
        # colnames(gei_coe) <- colnames(fixed_effects)
        # rownames(gei_coe) <- rann
        # gei_hf <- fixed_effects[(ncol(cov)+length(colsnp)+1):nrow(fixed_effects),]
        # gei_coe[1:nrow(gei_hf),] <- gei_hf
        #rownames(estimates) <- colnames(result$coefficients)
        aic <- AIC(result)
        return(list(gwas1D = gwassig1,covariates_coefficients = cov_coe ,genotype_coefficients = fe_coe,GEI_coefficients = gei_coe, AIC = aic))
      } else {

        # Treat covariates (if it exists) and snp as fixed effects, and the interaction between snp and environment as random effects, and use the glmer function to calculate
        formula <- as.formula(paste(coly, "~",paste(colx, collapse = "+"), "+", paste(colsnp, collapse = "+"),"+",
                                    paste(paste0("(0 + ", colsnp, "|", cole, ")", collapse = "+"))))

        if(family == "gaussian"){
          result <- lmer(formula, data = quan)
        }else if(family == "binomial" && !is.null(weights) ){
          result <- glmer(formula, data = quan,weights = frequency, family = family,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
        }else if (family == "negative.binomial") {
          result <- glmer(formula, data = quan, family = negative.binomial(theta),control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)) )
        }else{
          result <- glmer(formula, data = quan, family = family,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
        }
        fixed_effects <- coef(summary(result))

        cov_coe <- fixed_effects[!rownames(fixed_effects) %in% colsnp,]
        if (nrow(cov_coe) < length(colx)) {
          message("There are NAs in the estimated coefficients of covariates. The corresponding variables have been removed from the results.")
        }

        fe_coe <- matrix(NA,nrow = length(colsnp),ncol = ncol(fixed_effects))
        rownames(fe_coe) <- colsnp
        colnames(fe_coe) <- colnames(fixed_effects)
        ma1 <- colsnp %in% rownames(fixed_effects)
        ma2 <- rownames(fixed_effects) %in% colsnp
        fe_coe[ma1,] <- fixed_effects[ma2,]

        # colsnp_effects <- fixed_effects[colsnp,]
        # other_effects <- fixed_effects[1:ncol(cov),]
        random_effects <- VarCorr(result)
        #rownames(estimates) <- colnames(result$coefficients)
        rane <- ranef(result)
        rann <- c()

        for (snp in colsnp) {
          for (e in ename) {

            combination <- paste0(snp, ":", e)
            rann <- c(rann, combination)
          }
        }
        ranv<- c()

        for (env in names(rane)) {
          for (snp in names(rane[[env]])) {

            ranv <- c(ranv, rane[[env]][[snp]])
          }
        }
        names(ranv) <- rann
        aic <- AIC(result)
        return(list(gwas1D = gwassig1,covariates_coefficients = cov_coe, genotype_coefficients = fe_coe,
                    GEI_coefficients = ranv, GEI_varcorr = random_effects, AIC = aic))
      }
    }

  }else{

    return(list(Sr.No = fin,gwas1D = gwassig1))

  }
}

#' Filtering and parameter estimation for significant SNPs in 2D Scaning
#'
#' This function filters the results from an A-A epistatic 2D scaning to obtain
#' significant SNP markers and allows for parameter estimation.
#'
#' @param significantASNP a vector of significant SNP names from 1D scanning
#' @param phenotype a vector of phenotype values
#' @param covariates a matrix of covariates (PCA)
#' @param environment a vector of environmental factors
#' @param weights a vector of weights
#' @param data a data frame containing the input data
#' @param theta the theta parameter
#' @param genodata a genotype matrix
#' @param map the map file
#' @param Pvalues a vector of 2D p-values
#' @param Bonferr bonferroni correction factor, 0.5/(M * (M - 1)/2)
#' @param family the distribution family
#' @param lasso logical, LASSO method
#' @param selection logical, forward selection method
#' @param estimate logical, estimate parameters for significant markers
#' @param GEI_random logical, gene-environment interaction effects as random effects
#'
#' @return A list containing information about significant SNP markers from 1D and 2D scans,
#'         covariates coefficients, genotype coefficients, GEI coefficients, GEI VarCorr, and AIC.
#'
#' @export
#'
#' @examples
#'  \dontrun{
#' sig <- c("SNP25","SNP67")
#' raae1 <- epi.snp.estimate(significantASNP = sig,phenotype = y, covariates = cov,
#'                           environment = env, genodata = Ga,map = map,Pvalues = p2,
#'                           Bonferr = 0.5,family = "binomial",selection = T,estimate = T)}
#'
epi.snp.estimate <- function(significantASNP = NULL,phenotype = phenotype, covariates = NULL, environment = NULL,
                             weights = NULL,data = NULL,theta = NULL, genodata = genodata,map = map,Pvalues = Pvalues,
                             Bonferr = 0.5, family = "gaussian",lasso = FALSE,selection = FALSE,
                             estimate = FALSE,GEI_random = TRUE) {

  if (!is.null(significantASNP)) {
    if (!(is.vector(significantASNP) || (is.matrix(significantASNP) && ncol(significantASNP) == 1))) {
      stop("Error: If significantASNP is provided, it must be a vector or a matrix with ncol equal to 1.")
    }
  }

  weight <- numeric(2)
  # Process phenotype, covariates, and environment inputs
  if (!is.null(data)) {
    if (!is.data.frame(data)) {
      stop("The 'data' input must be a data frame.")
    }
    if(!phenotype %in% colnames(data)){
      stop("The 'phenotype' input must be a single column name present in the 'data' data frame.")
    }
    if (is.character(phenotype)) {
      phenotype <- data[[phenotype]]
    }
    if (!all(!is.na(phenotype) & phenotype != "" & phenotype != ".")) {
      stop("Phenotype contains NA, empty values, or dots.")
    }
    if (!is.null(covariates)) {
      if (is.character(covariates)) {
        covariates_names <- strsplit(covariates, "\\+")[[1]]
        covariates <- as.matrix(data[, covariates_names, drop = FALSE])
      }
      if (!all(!is.na(covariates) & covariates != "" & covariates != ".")) {
        stop("Covariates contains NA, empty values, or dots.")
      }
      cov <- as.matrix(cbind(1, covariates))
    } else {
      cov <- matrix(1, nrow = length(phenotype), ncol = 1)
    }

    if (!is.null(environment)) {
      if(!environment %in% colnames(data)){
        stop("The 'environment' input must be a single column name present in the 'data' data frame.")
      }
      environment <- data[[environment]]
      if (!all(!is.na(environment) & environment != "" & environment != ".")) {
        stop("Environment contains NA, empty values, or dots.")
      }
      envir <- length(unique(environment))
      environment <- factor(environment)
      ename <-  levels(environment)
    }else{
      envir <- 1
    }
    if (!is.null(weights)) {
      if(!weights %in% colnames(data)){
        stop("The 'weights' input must be a single column name present in the 'data' data frame.")
      }
      weight <- data[[weights]]
      if (!all(!is.na(weight) & weight != "" & weight != ".")) {
        stop("Weights contains NA, empty values, or dots.")
      }
    }
  } else {
    if (!(is.vector(phenotype) || (is.matrix(phenotype) && ncol(phenotype) == 1))) {
      stop("Error: phenotype must be a vector or a matrix with ncol equal to 1.")
    }

    phenotype <- phenotype

    if (!all(!is.na(phenotype) & phenotype != "" & phenotype != ".")) {
      stop("Phenotype contains NA, empty values, or dots.")
    }
    if (!is.null(covariates)) {
      if (!is.matrix(covariates)) {
        covariates <- as.matrix(covariates)
      }
      if (nrow(covariates) != length(phenotype)) {
        stop("Covariates length does not match the length of phenotype.")
      }

      if (!all(!is.na(covariates) & covariates != "" & covariates != ".")) {
        stop("Covariates contains NA, empty values, or dots.")
      }
      cov <- cbind(1, covariates)
    } else {
      cov <- matrix(1, nrow = length(phenotype), ncol = 1)
    }
    if (!is.null(environment)) {
      if (!(is.vector(environment) || (is.matrix(environment) && ncol(environment) == 1))) {
        stop("Error: environment must be a vector or a matrix with ncol equal to 1.")
      }
      environment <- environment
      if (length(environment) != length(phenotype)) {
        stop("Environment length does not match the length of phenotype.")
      }

      if (!all(!is.na(environment) & environment != "" & environment != ".")) {
        stop("Environment contains NA, empty values, or dots.")
      }
      environment <- factor(environment)
      envir <- length(unique(environment))
      ename <-  levels(environment)
    }else{
      envir <- 1
    }
    if (!is.null(weights)) {
      if (!(is.vector(weights) || (is.matrix(weights) && ncol(weights) == 1))) {
        stop("Error: weights must be a vector or a matrix with ncol equal to 1.")
      }
      weight <- weights
      if (length(weight) != length(phenotype)) {
        stop("Weights length does not match the length of phenotype.")
      }

      if (!all(!is.na(weight) & weight != "" & weight != ".")) {
        stop("Weights contains NA, empty values, or dots.")
      }
    }
  }
  #lencov <- ncol(cov)
  if ( is(genodata, "SnpMatrix")) {
    stop("The 'genodata' input cannot be in SnpMatrix format.")
  }

  # Check genodata dimensions
  if (nrow(genodata) != (length(phenotype)/envir)) {
    stop("The number of rows in genodata must match the length of phenotype.")
  }

  if (ncol(genodata) != nrow(map)) {
    stop("The number of columns in genodata must match the number of rows in map.")
  }

  # Check map column names
  required_colnames <- c("chromosome", "snp.name", "cM", "position", "allele.1", "allele.2")
  if (!all(colnames(map) == required_colnames)) {
    stop("The column names of map must be: 'chromosome', 'snp.name', 'cM', 'position', 'allele.1', 'allele.2'. Please rename the columns accordingly.")
  }

  # Check if the 'snp.name' column in 'map' matches the column names of 'genodata' in both length and order
  if (!identical(map$snp.name, colnames(genodata))) {
    stop("The 'snp.name' column in 'map' does not match the column names of 'genodata' in terms of order.")
  }

  all_families <- c("gaussian", "binomial", "poisson", "negative.binomial")
  if (!family %in% all_families) {
    stop("Invalid family parameter. Supported families are 'gaussian', 'poisson',
         'binomial' and 'negative.binomial'.")
  }
  phenotype1 <- phenotype
  if (family == "gaussian") {
    phenotype1 <- scale(phenotype)
  }
  if (family == "binomial") {
    if (any(phenotype > 1)) {
      if (is.null(weights)) {
        stop("Phenotype data for Binomial family must have weights specified for counts greater than 1.")
      }else{
        phenotype <- phenotype / weight
        phenotype1 <- phenotype
      }
    } else if (any(phenotype < 0) || any(phenotype > 1)) {
      stop("Phenotype data for Binomial family must only contain values between 0 and 1.")
    }

    # If weights are NULL, generate weight vector
    if (is.null(weights) && any(phenotype >= 0 & phenotype <= 1)) {
      if (!all(phenotype %in% c(0, 1))) {
        weight <- rep(100, length(phenotype))
        weights <- weight
      }
    }
  }
  if (family == "poisson" || family == "negative.binomial") {
    if (any(phenotype < 0)) {
      stop("Phenotype data for '", family, "' family must be non-negative.")
    }
  }
  #gwassig1 <- ifelse(!is.null(significantASNP), matrix(NA, nrow = length(significantASNP), ncol = 3), NULL)
  n <- nrow(genodata)
  if (!is.null(significantASNP)) {
    Xc <- matrix(nrow = n, ncol = length(significantASNP))

    for (i in 1:length(significantASNP)) {
      cn <- significantASNP[i]
      ci <- which(colnames(genodata) == cn)
      Xc[, i] <- genodata[, ci]
    }
    Xcc <- matrix(nrow = (n*envir), ncol = (length(significantASNP)*envir))
    Xcc <- generateSubmatrix(Xc, envir)

    cov <- cbind(cov,Xcc)
    gwassig1 <- matrix(NA, nrow = length(significantASNP), ncol = 3)
    significantASNP <- as.matrix(significantASNP)
    gwassig1[,1] <- significantASNP
    gwassig1[,2] <- map[match(significantASNP, map$snp.name),1]
    gwassig1[,3] <- map[match(significantASNP, map$snp.name),4]
    colnames(gwassig1) <- c("SNP", "CHR", "BP")

  }

  if (family == "negative.binomial") {
    if (is.null(theta)) {

      ml <- glm.nb(phenotype ~ cov)
      theta <- ml$theta

    }else{
      theta <- theta
    }

  }else{
    theta <- 0.0
  }
  m <- ncol(genodata)
  num <- m * (m - 1) / 2
  #Pvalues <- Pvalues
  #Pvalues <- p.adjust(Pvalues, method = "fdr")
  thr <- 0.5 / num
  thr1 <- Bonferr / num

  significant_l_indices <- which(Pvalues <= thr1)
  if (length(significant_l_indices) == 0) {
    miss <- NULL
    message("No significant markers found.")
    return(miss)
  }


  number <- matrix(NA, nrow = num, ncol = 2)
  col_idx <- 1

  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      number[col_idx, 1] <- i
      number[col_idx, 2] <- j
      col_idx <- col_idx + 1
    }
  }
  Gaa <- matrix(NA, nrow = n, ncol = length(significant_l_indices))
  SNP<- character(length(significant_l_indices))
  col_idx <- 1

  for (index in significant_l_indices) {
    i <- number[index, 1]
    j <- number[index, 2]

    SNP[col_idx] <- paste0(colnames(genodata)[i], "_", colnames(genodata)[j])
    Gaa[, col_idx] <- genodata[, i] * genodata[, j]

    col_idx <- col_idx + 1
  }
  colnames(Gaa) <- SNP

  if(lasso){
    if (family == "negative.binomial" ) {
      message("Lasso method is not applicable for negative binomial models. Performing Bonferroni correction instead.")
      fin <- 1:length(significant_l_indices)
    }else{
      if(length(significant_l_indices) >= 2 ){
        fin_l <- 1:length(significant_l_indices)
        gs <- as.matrix(Gaa)
        if(envir == 1){

          colnames(gs)  <- fin_l

        }else{

          gs <- generateSubmatrix(gs,envir)
          gro <- rep(fin_l,envir)
          colnames(gs)  <- gro
        }
        if (family == "binomial" && !is.null(weights) ){
          phe1 <- phenotype * weight
          phe2 <- weight - phe1
          pheno <- cbind(phe1,phe2)
          cvobs <-  cv.glmnet(data.matrix(gs), as.matrix(pheno),family = family)
          obs <-  glmnet(data.matrix(gs), as.matrix(pheno),family = family,lambda = cvobs$lambda.min )

        }else{
          cvobs <-  cv.glmnet(gs, phenotype, family = family)
          obs <-  glmnet(gs, phenotype,  family = family,lambda = cvobs$lambda.min )
        }
        fin <- as.numeric(rownames(obs$beta)[as.numeric(obs$beta)>0])
        fin <- unique(fin)
      }else{
        fin <- 1:length(significant_l_indices)
      }
    }

  }

  if (selection) {

    thr2 <- (Bonferr * 0.01) / num
    significant_s_indices <- which(Pvalues <= thr2)
    if (length(significant_s_indices) == 0) {
      thr2 <- (Bonferr * 0.1) / num
      significant_s_indices <- which(Pvalues <= thr2)
      if (length(significant_s_indices) == 0) {
        fin <- 1:length(significant_l_indices)
      }else if (identical(significant_s_indices, significant_l_indices)) {
        fin <- 1:length(significant_l_indices)
      } else {
        fin_l <- 1:length(significant_l_indices)
        fin_s <- match(significant_s_indices, significant_l_indices)
        fin_l <- fin_l - 1
        fin_s <- fin_s - 1

        fin <- forward_selection(phenotype1,weight, cov, Gaa,family, envir, thr,theta, fin_s, fin_l)
        fin <- fin + 1
      }
    }else if (identical(significant_s_indices, significant_l_indices)) {
      fin <- 1:length(significant_l_indices)
    } else {
      fin_l <- 1:length(significant_l_indices)
      fin_s <- match(significant_s_indices, significant_l_indices)
      fin_l <- fin_l - 1
      fin_s <- fin_s - 1

      fin <- forward_selection(phenotype1,weight, cov, Gaa,family, envir,thr, theta, fin_s, fin_l)
      fin <- fin + 1
    }

  }

  if(!selection){
    fin <- 1:length(significant_l_indices)
  }

  sli <- significant_l_indices[fin]
  gwassig2 <- matrix(NA,nrow = length(fin),ncol = 7)
  col_idx <- 1
  for (index in sli) {
    i <- number[index, 1]
    j <- number[index, 2]

    gwassig2[col_idx,1] <- colnames(genodata)[i]
    gwassig2[col_idx,2] <- map[i,1]
    gwassig2[col_idx,3] <- map[i,4]
    gwassig2[col_idx,4] <- colnames(genodata)[j]
    gwassig2[col_idx,5] <- map[j,1]
    gwassig2[col_idx,6] <- map[j,4]

    col_idx <- col_idx + 1
  }
  gwassig2[,7] <- Pvalues[sli]
  colnames( gwassig2) <- c("SNP1","CHR1", "BP1","SNP2","CHR2", "BP2","P")

  if (!is.null(significantASNP)) {
    Xcn <- cbind(Xc,Gaa[, fin])
    allga <- matrix(NA, nrow = n * envir, ncol = ncol(Xcn))
    #allga <- matrix(NA, nrow = nrow(genodata) * envir, ncol = length(fin) )
    #Xcn <- matrix(NA, nrow = nrow(genodata) * envir,ncol =length(significantASNP))
    for (i in 1:envir) {
      allga[((i - 1) * n + 1):(i * n), ] <- Xcn
      #Xcn[((i - 1) * nrow(genodata) + 1):(i * nrow(genodata)), ] <- Xc
      #allga[((i - 1) * nrow(genodata) + 1):(i * nrow(genodata)), ] <- XcnGaa[, fin]
    }
    #Xc <- rbind(Xc,Xc)
    #allga <- cbind(Xc,allga)
    #allga <- cbind(Xcn,allga)
    colsnp <- SNP[fin]
    colsnp <- c(significantASNP,colsnp)
  }else{
    allga <- matrix(NA, nrow = n * envir, ncol = length(fin))
    for (i in 1:envir) {
      allga[((i - 1) * n + 1):(i * n), ] <- Gaa[, fin]
    }
    colsnp <- SNP[fin]
  }
  if(estimate){
    coly <- c("y")

    if (is.null(covariates) && is.null(environment)) {
      quan <- as.data.frame(cbind(phenotype, allga))
      colnames(quan) <- c(coly, colsnp)
    } else if (is.null(covariates)) {
      quan <- as.data.frame(cbind(phenotype, environment, allga))
      cole <- c("environment")
      colnames(quan) <- c(coly,cole, colsnp)
    } else if (is.null(environment)) {
      colx <- colnames(covariates)

      quan <- as.data.frame(cbind(phenotype, covariates, allga))
      colnames(quan) <- c(coly,colx, colsnp)
    } else {
      colx <- colnames(covariates)
      cole <- c("environment")
      quan <- as.data.frame(cbind(phenotype, covariates, environment, allga))
      colnames(quan) <- c(coly, colx, cole, colsnp)
    }
    if (!is.null(weights)){
      weight <- as.matrix(weight)
      colw <- c("frequency")
      colnames(weight) <- colw
      quan <- cbind(quan,weight)
    }
    colx <- if (!is.null(covariates)) colnames(covariates) else 1

    if (envir == 1) {
      formula <- as.formula(paste(coly, "~",paste(colx, collapse = "+"), "+", paste(colsnp, collapse = "+")))

      if (family == "binomial" && !is.null(weights) ){
        # formula <- as.formula(paste("cbind(",coly,",",colw,"-",coly,")", "~",paste(colx, collapse = "+"), "+", paste(colsnp, collapse = "+")))
        result <- glm(formula, data = quan,weights = frequency, family = family)
      }else if (family == "negative.binomial") {
        result <- glm(formula, data = quan, family = negative.binomial(theta))
      }else{
        result <- glm(formula, data = quan, family = family)
      }
      estimates <- coef(summary(result))
      cov_coe <- estimates[!rownames(estimates) %in% colsnp,]
      if (nrow(cov_coe) < length(colx)) {
        message("There are NAs in the estimated coefficients of covariates. The corresponding variables have been removed from the results.")
      }
      #coe2 <- estimates[rownames(estimates) %in% colsnp,]

      fe_coe <- matrix(NA,nrow = length(colsnp),ncol = ncol(estimates))
      rownames(fe_coe) <- colsnp
      colnames(fe_coe) <- colnames(estimates)
      ma1 <- colsnp %in% rownames(estimates)
      ma2 <- rownames(estimates) %in% colsnp
      fe_coe[ma1,] <- estimates[ma2,]
      #fe_coe <- estimates[rownames(estimates) %in% colsnp,]
      #fe_coe <- estimates[colsnp,]

      aic <- AIC(result)
      if (!is.null(significantASNP)) {
        return(list(gwas1D = gwassig1,gwas2D = gwassig2,covariates_coefficients = cov_coe ,genotype_coefficients = fe_coe, AIC = aic))
      }else{
        return(list(gwas2D = gwassig2,covariates_coefficients = cov_coe ,genotype_coefficients = fe_coe, AIC = aic))
      }
    } else {
      if (!GEI_random) {

        gf <- as.matrix(allga[1:n, ])
        gf <- generateSubmatrix(gf, envir)

        rann <- c()
        for (e in ename) {
          for (snp in colsnp) {
            combination <- paste0(snp, "_", e)
            rann <- c(rann, combination)
          }
        }
        colnames(gf) <- rann

        quan <- cbind(quan, gf)
        formula <- as.formula(paste(coly, "~", paste(colx, collapse = "+"), "+", paste(colsnp, collapse = "+"),
                                    "+",paste(rann, collapse = "+")))
        if (family == "binomial" && !is.null(weights) ){
          result <- glm(formula, data = quan,weights = frequency, family = family)
        }else if (family == "negative.binomial") {
          result <- glm(formula, data = quan, family = negative.binomial(theta))
        }else{
          result <- glm(formula, data = quan, family = family)
        }

        fixed_effects <- coef(summary(result))
        cov_coe <- fixed_effects[!rownames(fixed_effects) %in% colsnp,]
        if (nrow(cov_coe) < length(colx)) {
          message("There are NAs in the estimated coefficients of covariates. The corresponding variables have been removed from the results.")
        }
        #cov_coe <- fixed_effects[1:lencov,]

        fe_coe <- matrix(NA,nrow = length(colsnp),ncol = ncol(fixed_effects))
        rownames(fe_coe) <- colsnp
        colnames(fe_coe) <- colnames(fixed_effects)
        ma1 <- colsnp %in% rownames(fixed_effects)
        ma2 <- rownames(fixed_effects) %in% colsnp
        fe_coe[ma1,] <- fixed_effects[ma2,]
        #fe_coe <- fixed_effects[colsnp,]

        gei_coe <- matrix(NA,nrow = length(rann),ncol = ncol(fixed_effects))
        colnames(gei_coe) <- colnames(fixed_effects)
        rownames(gei_coe) <- rann
        ma3 <- rann %in% rownames(fixed_effects)
        ma4 <- rownames(fixed_effects) %in% rann
        gei_coe[ma3,] <- fixed_effects[ma4,]
        #gei_hf <- fixed_effects[(lencov+length(colsnp)+1):nrow(fixed_effects),]
        #gei_coe[1:nrow(gei_hf),] <-  gei_hf

        #rownames(estimates) <- colnames(result$coefficients)
        aic <- AIC(result)
        if (!is.null(significantASNP)) {
          return(list(gwas1D = gwassig1,gwas2D = gwassig2,covariates_coefficients = cov_coe ,genotype_coefficients = fe_coe,GEI_coefficients = gei_coe, AIC = aic))
        }else{
          return(list(gwas2D = gwassig2,covariates_coefficients = cov_coe ,genotype_coefficients = fe_coe,GEI_coefficients = gei_coe, AIC = aic))
        }
      } else {

        formula <- as.formula(paste(coly, "~",paste(colx, collapse = "+"), "+", paste(colsnp, collapse = "+"),"+",
                                    paste(paste0("(0 + ", colsnp, "|", cole, ")", collapse = "+"))))

        if(family == "gaussian"){
          result <- lmer(formula, data = quan)
        }else if(family == "binomial" && !is.null(weights) ){
          result <- glmer(formula, data = quan,weights = frequency, family = family,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
        }else if (family == "negative.binomial") {
          result <- glmer(formula, data = quan, family = negative.binomial(theta),control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
        }else{
          result <- glmer(formula, data = quan, family = family,control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
        }
        fixed_effects <- coef(summary(result))

        cov_coe <- fixed_effects[!rownames(fixed_effects) %in% colsnp,]
        if (nrow(cov_coe) < length(colx)) {
          message("There are NAs in the estimated coefficients of covariates. The corresponding variables have been removed from the results.")
        }
        #cov_coe <- fixed_effects[1:lencov,]

        fe_coe <- matrix(NA,nrow = length(colsnp),ncol = ncol(fixed_effects))
        rownames(fe_coe) <- colsnp
        colnames(fe_coe) <- colnames(fixed_effects)
        ma1 <- colsnp %in% rownames(fixed_effects)
        ma2 <- rownames(fixed_effects) %in% colsnp
        fe_coe[ma1,] <- fixed_effects[ma2,]
        #colsnp_effects <- fixed_effects[colsnp,]

        random_effects <- VarCorr(result)
        #rownames(estimates) <- colnames(result$coefficients)
        rane <- ranef(result)
        rann <- c()

        for (snp in colsnp) {
          for (e in ename) {

            combination <- paste0(snp, ":", e)
            rann <- c(rann, combination)
          }
        }
        ranv<- c()

        for (env in names(rane)) {
          for (snp in names(rane[[env]])) {

            ranv <- c(ranv, rane[[env]][[snp]])
          }
        }
        names(ranv) <- rann
        aic <- AIC(result)
        e_vector <- significant_l_indices[fin]
        if (!is.null(significantASNP)) {
          return(list(gwas1D = gwassig1,Sr.No = e_vector,gwas2D = gwassig2,covariates_coefficients = cov_coe, genotype_coefficients = fe_coe,
                      GEI_coefficients = ranv, GEI_varcorr = random_effects, AIC = aic))
        }else{
          return(list(Sr.No = e_vector,gwas2D = gwassig2,covariates_coefficients = cov_coe, genotype_coefficients = fe_coe,
                      GEI_coefficients = ranv, GEI_varcorr = random_effects, AIC = aic))
        }

      }
    }

  }else{
    e_vector <- significant_l_indices[fin]
    if (!is.null(significantASNP)) {
      return(list(gwas1D = gwassig1,Sr.No = e_vector,gwas2D = gwassig2))
    }else{
      return(list(Sr.No = e_vector,gwas2D = gwassig2))
    }
  }

}
