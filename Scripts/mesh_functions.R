# mesh_functions.R




get_spde_matrices <- function(x) {
  x <- x$spde[c("c0", "g1", "g2")]
  names(x) <- c("M0", "M1", "M2") # legacy INLA names needed!
  x
}

# from TMB examples repository:
make_anisotropy_spde <- function(spde, anistropy = TRUE) {
  if (anistropy) {
    inla_mesh <- spde$mesh
    Dset <- 1:2
    inla_mesh <- spde$mesh
    # Triangle info
    TV <- inla_mesh$graph$tv # Triangle to vertex indexing
    V0 <- inla_mesh$loc[TV[, 1], Dset] # V = vertices for each triangle
    V1 <- inla_mesh$loc[TV[, 2], Dset]
    V2 <- inla_mesh$loc[TV[, 3], Dset]
    E0 <- V2 - V1 # E = edge for each triangle
    E1 <- V0 - V2
    E2 <- V1 - V0
    # Calculate Areas
    TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
    Tri_Area <- rep(NA, nrow(E0))
    for (i in seq_len(length(Tri_Area))) Tri_Area[i] <- TmpFn(E0[i, ], E1[i, ]) / 2
    
    ret <- list(
      n_s = spde$mesh$n,
      n_tri = nrow(TV),
      Tri_Area = Tri_Area,
      E0 = E0,
      E1 = E1,
      E2 = E2,
      TV = TV - 1,
      G0 = spde$spde$c0,
      G0_inv = as(Matrix::diag(1 / Matrix::diag(spde$spde$c0)), "TsparseMatrix")
    )
  } else {
    ret <- list(
      n_s = 0L, n_tri = 0L, Tri_Area = rep(0, 1), E0 = matrix(0, 1),
      E1 = matrix(0, 1), E2 = matrix(0, 1), TV = matrix(0, 1),
      G0 = Matrix::Matrix(0, 1, 1, doDiag = FALSE), G0_inv = Matrix::Matrix(0, 1, 1, doDiag = FALSE)
    )
  }
  ret
}