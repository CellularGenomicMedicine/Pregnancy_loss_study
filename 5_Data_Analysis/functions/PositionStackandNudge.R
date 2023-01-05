position_stack_and_nudge <- function(x = 0, y = 0, vjust = 1, reverse = FALSE) {
  
  ggproto(NULL, PositionStackAndNudge,
          x = x,
          y = y,
          vjust = vjust,
          reverse = reverse
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @noRd
PositionStackAndNudge <- ggproto("PositionStackAndNudge", PositionStack,
                                 x = 0,
                                 y = 0,
                                 
                                 setup_params = function(self, data) {
                                   c(
                                     list(x = self$x, y = self$y),
                                     ggproto_parent(PositionStack, self)$setup_params(data)
                                   )
                                 },
                                 
                                 compute_layer = function(self, data, params, panel) {
                                   # operate on the stacked positions (updated in August 2020)
                                   data = ggproto_parent(PositionStack, self)$compute_layer(data, params, panel)
                                   
                                   x_orig <- data$x
                                   y_orig <- data$y
                                   # transform only the dimensions for which non-zero nudging is requested
                                   if (any(params$x != 0)) {
                                     if (any(params$y != 0)) {
                                       data <- transform_position(data, function(x) x + params$x, function(y) y + params$y)
                                     } else {
                                       data <- transform_position(data, function(x) x + params$x, NULL)
                                     }
                                   } else if (any(params$y != 0)) {
                                     data <- transform_position(data, function(x) x, function(y) y + params$y)
                                   }
                                   data$nudge_x <- data$x
                                   data$nudge_y <- data$y
                                   data$x <- x_orig
                                   data$y <- y_orig
                                   
                                   data
                                 },
                                 
                                 compute_panel = function(self, data, params, scales) {
                                   ggproto_parent(PositionStack, self)$compute_panel(data, params, scales)
                                 }
)