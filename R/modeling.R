# Modeling functions.

# Generate a single kinetic model ------

#' Generate a kinetic object model.
#'
#' @description Generates a mixed-effect n-order model and embeds it
#' in a kinetic class object.
#' @details Technically, a wrapper around \code{\link[lme4]{lmer}}
#' (if family is 'gaussian') and \code{\link[lme4]{glmer}} functions.
#' @param data a data frame.
#' @param response name of the dependent variable.
#' @param time name of the time variable.
#' @param ID name of the pairing/grouping or block variable.
#' @param family modeling family.
#' @param model order.
#' @param ... extra arguments passed to \code{\link[lme4]{lmer}} or
#' \code{\link[lme4]{glmer}}.
#' @return a kinetic class object.
#' @export

  model_kinetic <- function(data,
                            response,
                            time = 'time',
                            ID = 'ID',
                            family = 'gaussian',
                            order = 2, ...) {

    ## entry control

    if(!any(class(data) == 'data.frame')) {

      stop('A data frame needed as data argument', call. = FALSE)

    }

    if(!any(c(response, time, ID) %in% names(data))) {

      stop('A variable absent from the input data frame.', call. = FALSE)

    }

    order <- as.integer(order)

    ## model formula

    if(order == 0) {

      mod_formula <- paste(response, paste0('(1|', ID, ')'), sep = '~')

    } else {

      mod_formula <- paste0(response, '~', time)

      if(order == 1) {

        mod_formula <- paste0(mod_formula, '+', paste0('(1|', ID, ')'))

      } else {

        mod_formula <- paste0(mod_formula, '+',
                              paste0('I(', time, '^', 2:order, ')', collapse = '+'), '+',
                              paste0('(1|', ID, ')'))

      }

    }

    mod_formula <- as.formula(mod_formula)

    ## modeling

    if(family == 'gaussian') {

      model <- lme4::lmer(formula = mod_formula,
                          data = data, ...)

    } else {

      model <- lme4::glmer(formula = mod_formula,
                           data = data,
                           family = family, ...)

    }

    kinet::kinetic(model = model,
                   response = response,
                   ID = ID,
                   time = time,
                   family = family,
                   order = order)

  }

# Kinetic object list ------

#' Generate a list of kinetic objects.
#'
#' @description Generates a list of kinetic objects for a vector of response
#' variable names. A wrapper around \code{\link{model_kinetic}} enabling for
#' model construction in parallel.
#' @param response a vector with the response variable names.
#' @param .parallel logical, should the model construction be done in parallel?
#' @param .errors logical, should modeling failures be returned?
#' @inheritParams model_kinetic
#' @return a list with kinetic class objects, one for each response variable.
#' @export

  model_kinetic_lst <- function(data,
                                responses = NULL,
                                time = 'time',
                                ID = 'ID',
                                family = 'gaussian',
                                order = 2,
                                .parallel = FALSE,
                                .errors = FALSE, ...) {

    ## entry control

    stopifnot(is.logical(.parallel))
    stopifnot(is.logical(.errors))

    ## Benchmarking

    start_time <- Sys.time()
    message(paste('Creating', length(responses), 'kinetic models'))
    on.exit(message(paste('Elapsed', Sys.time() - start_time)))

    if(.parallel) {

      future::plan('multisession')

      kinet_lst <- furrr::future_map(responses,
                                     purrr::safely(kinet::model_kinetic),
                                     data = data,
                                     time = time,
                                     ID = ID,
                                     order = order,
                                     .options = furrr::furrr_options(seed = TRUE,
                                                                     packages = c('lme4',
                                                                                  'rlang')),
                                     ...)

      future::plan('sequential')

    } else {

      kinet_lst <- purrr::map(responses,
                              purrr::safely(kinet::model_kinetic),
                              data = data,
                              time = time,
                              ID = ID,
                              order = order, ...)

    }

    kinet_lst <- purrr::transpose(rlang::set_names(kinet_lst,
                                                   responses))

    kinet_lst <- purrr::map(kinet_lst, purrr::compact)

    if(.errors) {

      return(kinet_lst)

    } else {

      return(kinet_lst$result)

    }

  }

# END -----
