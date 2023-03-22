interpol.TPS <- function (x, y)
{
  distance <- Norm (x - y)
  if (length (x) == 2)
  {
    if (distance < .Machine $ double.eps)
      return (0)
    else
      return (- log (distance ^ 2) * (distance ^ 2))
  }
  if (length (x) == 3)
    return (distance)
}

interpol.L1.TPS <- function (y, X)
  return (aaply (X, 1, interpol.TPS, y = y))

interpol.L2.TPS <- function (X, Y)
  return (aaply (Y, 1, interpol.L1.TPS, X = X))

#' TPS transform
#'
#' Calculates the Thin Plate Spline transform between a reference shape and a target shape
#'
#' @param target.shape Target shape
#' @param reference.shape Reference shape
#'
#' @returns A list with the transformation parameters and a function that gives
#' the value of the TPS function at each point for numerical differentiation
#'
#' @author Guilherme Garcia
#' @export
#' @importFrom stats spline
TPS <- function (target.shape, reference.shape)
{
  ## calculates thin plate splines deformations between reference and target shapes
  theta.p <- function (p)
    return (b + (A %*% p) + t (W) %*% interpol.L1.TPS (p, Q))
  Q <- reference.shape
  P <- target.shape
  k <- nrow (Q)
  m <- ncol (Q)
  R <- interpol.L2.TPS (Q, Q)
  # print (R)
  one.k <- array (1, c (k, 1))
  zero.11 <- array (0, c (1, 1))
  zero.m1 <- array (0, c (m, 1))
  zero.mm <- array (0, c (m, m))
  L <- rbind (R, t (one.k), t (Q))
  L <- cbind (L, rbind (one.k, zero.11, zero.m1))
  L <- cbind (L, rbind (Q, t (zero.m1), zero.mm))
  P.zero <- rbind (P, t (zero.m1), zero.mm)
  spline <- solve (L, P.zero)
  W <- spline [1:k, ]
  b <- t (t (spline [1+k, ]))
  A <- t (spline [(2:(m+1))+k, ])
  return (list ('W' = W, 'A' = A,
                'b' = b, 'Q' = Q, 'P' = P, 'theta.p' = theta.p))
}


#' Midline rotate
#'
#' Returns the rotation matrix that aligns a specimen sagital line
#' to plane y = 0 (2D) or z = 0 (3D)
#'
#' @param X shape array
#' @param midline rows for the midline landmarks
#'
#' @returns Rotation matrix
#'
#' @author Guilherme Garcia
#' @export
Rotate2MidlineMatrix <- function (X, midline)
{
  ncl <- ncol (X)
  Xm <- na.omit (X [midline, ])
  Mm <- matrix (apply (Xm, 2, mean), byrow = TRUE, nrow = nrow (X), ncol = ncl)
  Xc <- X - Mm
  W <- na.omit (Xc [midline, ])
  RM <-svd (var (W))$v
  return (RM)
}

#' Local Jacobian calculation
#'
#' Calculates jacobians for a given interpolation in a set of points
#' determined from tesselation (as centroids of each tetrahedron defined, for now...)
#'
#' @param spline Thin plate spline calculated by the TPS function
#' @param tesselation matrix of landmarks.
#' @param ... Additional arguments to some function
#' @note Jacobians are calculated on the row centroids of the tesselation matrix.
#'
#' @returns array of jacobians calculated at the centroids
#'
#' @author Guilherme Garcia
#' @export
#' @importFrom numDeriv jacobian
JacobianArray <- function (spline, tesselation, ...)
{
  with (spline,
        {
          Q.tetra <- Q [tesselation, ]
          dim (Q.tetra) <- c (dim (tesselation), ncol (Q))
          Q.centroids <- apply (Q.tetra, 1, colMeans)
          jacobs <- apply (Q.centroids, 2, jacobian, func = theta.p, ...)
          dim (jacobs) <- c (ncol (Q), ncol (Q), ncol (Q.centroids))
          return (jacobs)
        })
}

#' Centered jacobian residuals
#'
#' Calculates mean jacobian matrix for a set of jacobian matrices
#' describing a local aspect of shape deformation for a given set of volumes,
#' returning log determinants of deviations from mean jacobian (Woods, 2003).
#'
#' @param jacobArray Arrays of Jacobian calculated in the JacobianArray function
#'
#' @returns array of centered residual jacobians
#'
#' @author Guilherme Garcia
#' @author Diogo Melo
#' @references Woods, Roger P. 2003. “Characterizing Volume and Surface
#' Deformations in an Atlas Framework: Theory, Applications, and Implementation.” NeuroImage 18 (3):769-88.
#' @export
#' @importFrom expm logm
Center2MeanJacobianFast <- function (jacobArray)
{
  logm.single <- function (Ai, inv.Mk) return (logm (Ai %*% inv.Mk))
  A <- jacobArray
  N <- dim (A) [3]
  Mk <- MeanMatrix(A, tol = .Machine$double.eps)

  inv.Mk <- solve (Mk)
  centered.now <- apply (A, 3, logm.single, inv.Mk = inv.Mk)

  centered.now <- array (centered.now, c(nrow (A), nrow (A), N))
  log.det <- apply (centered.now, 3, function (x) return (sum (diag (x))))
  return (log.det)
}

#' Local Shape Variables
#'
#' Calculates the local shape variables of a set of landmarks using the sequence:
#' - TPS transform between all shapes and the mean shape
#' - Jacobian of the TPS transforms at the centroid of rows of the landmarks in the tesselation argument
#' - Mean center the Jacobians using the Karcher Mean
#' - Take the determinant of the centered jacobians
#'
#' @param gpa Procustes aligned landmarks.
#' @param cs Centoid sizes
#' @param landmarks unaligned landmarks. Ignored if both gpa and cs are passed.
#' @param tesselation matrix of rows of the landmarks. The centroid of each row
#' is used to mark the position of the jacobians
#' @param run_parallel Logical. If computation should be paralleled. Use with
#' caution, can make things worse. Requires that at parallel back-end like doMC
#' be registered
#' @returns List with TPS functions, jacobian matrices, local shape variables, mean shape, centroid sizes and individual IDs
#' @author Guilherme Garcia
#' @author Diogo Melo
#' @export
#' @importFrom Morpho ProcGPA cSize arrMean3
#' @importFrom plyr alply laply aaply
LocalShapeVariables <- function (gpa = NULL, cs = NULL, landmarks = NULL, tesselation,
                         run_parallel = FALSE)
{
  if(is.null(gpa) & is.null(cs)){
    cat("Running GPA")
    gpa <- ProcGPA(landmarks, CSinit = TRUE)
    cs <- apply(landmarks, 3, cSize)
    gpa <- gpa$rotated
    dimnames(gpa)[[3]] = dimnames(landmarks)[[3]]
  }
  mshape <- arrMean3(gpa)
  dimnames(mshape) <- dimnames(gpa)[1:2]
  tps <- alply (gpa, 3, TPS,
                reference.shape = mshape,
                .parallel = run_parallel, .progress = "text")
  print ('tps done')
  jacobs <- laply (tps, JacobianArray, tesselation = tesselation,
                   .parallel = run_parallel, .progress = "text")
  jacobs <- aperm (jacobs, c(2, 3, 1, 4), resize = TRUE)
  print ('jacobs done')
  local <- aaply (jacobs, 4, Center2MeanJacobianFast,
                  .parallel = run_parallel, .progress = "text")
  local <- t (local)
  return (list ('tps' = tps,
                'jacobians' = jacobs,
                'local' = local,
                'reference' = mshape,
                'cs' = cs,
                'info' = dimnames(gpa)[[3]]))
}
