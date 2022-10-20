#' Title
#'
#' @param ts
#' @param method
#' @param mss
#' @param pt
#' @param is_Euclidean_distance
#' @param distance
#' @param block_size
#' @importFrom stats dist
#' @import Rcpp
#'
#' @return
#' @export
#' @useDynLib CP, .registration = TRUE
#' @examples
CP = function(ts,
              method = 'origin',
              mss = 3,
              pt = 0.05,
              is_Euclidean_distance = TRUE,
              distance = NULL,
              block_size = 1)
{
  if (is.matrix(ts) != TRUE)
    ts = as.matrix(ts)

  ts_len = dim(ts)[1]
  ts_d = dim(ts)[2]
  ts_copy = as.vector(t(ts))
  change_point_list = rep(0, ts_len)
  change_point_num = c(0)

  if (method == "origin")
  {
    method_char = 'o'
  } else if (method == "gss")
  {
    method_char = 'g'
  } else if (method == "powell")
  {
    method_char = 'p'
  } else
  {
    stop("Method should be \'origin\',\'gss\' or \'powell\'")
  }

  if (mss %% 1 != 0)
  {
    stop("min segment size should be int")
  }

  if (ts_len <= 2 * mss)
  {
    stop("time series length should be more than 2 * min_segment_size")
  }

  if ((pt <= 0) | (pt >= 1))
  {
    stop("pt should more than 0 and less than 1")
  }

  if (!is_Euclidean_distance %in% c(TRUE, FALSE))
  {
    stop("is_Euclidean_distance should be bool")
  }




  if (is_Euclidean_distance == TRUE) {
    dist_copy = dist(ts,
                     method = "euclidean",
                     upper = TRUE,
                     diag = TRUE)
    dist_copy = as.vector(t(as.matrix(dist_copy)))
  } else {
    if (is.null(distance))
    {
      stop("distance should be given if is_Euclidean_distance is False")
    } else {
      if ((dim(distance)[1] == ts_len) & (dim(distance)[2] == ts_len)) {
        dist_copy = as.vector(t(as.matrix(distance)))
      } else {
        stop(paste0("distance should have shape (", ts_len, ", ", ts_len, ")"))
      }
    }
  }

  change_point_set = cp_run(
    ts_copy,
    ts_len,
    ts_d,
    mss,
    pt,
    is_Euclidean_distance,
    dist_copy,
    block_size,
    change_point_list,
    change_point_num,
    method_char
  )

  return (change_point_set[-1])

}
