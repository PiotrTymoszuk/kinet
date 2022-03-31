# Model inference, fit goodness and likelihood-ratio test.

# Inference for a kinetic object list -----

#' Inference stats for the kinetic object list.
#'
#' @description Calculates inference stats for a list of kinetic objects,
#' optionally in parallel. See: \code{\link{coef.kinetic}} for details.
#' @param kinetic_objects a list of kinetic class objects.
#' @param .parallel logical, should the model construction be done in parallel?
#' @param .errors logical, should modeling failures be returned?
#' @param ... extra arguments passed to \code{\link{coef.kinetic}}.
#' @return a list of data frames with the inference stats.
#' @export

  inference_kinetic_list <- function(kinetic_objects,
                                     .parallel = FALSE,
                                     .errors = FALSE, ...) {

    ## entry control

    if(!is.list(kinetic_objects)) {

      stop('The function requires a list of kinetic objects.',
           call. = TRUE)

    }

    class_check <- purrr::map_lgl(kinetic_objects, kinet::is_kinetic)

    if(any(!class_check)) {

      stop('Non-kinetic class objects provided.',
           call. = TRUE)

    }

    stopifnot(is.logical(.parallel))
    stopifnot(is.logical(.errors))

    ## Benchmarking

    start_time <- Sys.time()
    message(paste('Inference for',
                  length(kinetic_objects),
                  'kinetic model objects'))
    on.exit(message(paste('Elapsed', Sys.time() - start_time)))

    ## Inference

    if(.parallel) {

      future::plan('multisession')

      coef_lst<-  furrr::future_map(kinetic_objects,
                                    purrr::safely(kinet:::coef.kinetic),
                                    .options = furrr::furrr_options(packages = c('dplyr',
                                                                                 'purrr',
                                                                                 'rlang',
                                                                                 'lme4',
                                                                                 'lmerTest'),
                                                                    seed = TRUE),
                                    ...)

      future::plan('sequential')

    } else {

      coef_lst<- purrr::map(kinetic_objects,
                            purrr::safely(kinet:::coef.kinetic), ...)

    }

    coef_lst <- purrr::transpose(coef_lst)

    coef_lst <- purrr::map(coef_lst, purrr::compact)

    if(.errors) {

      return(coef_lst)

    } else {

      return(coef_lst$result)

    }

  }

# Goodness of fit and normality of the residuals ------

#' Fit stats for a kinetic class object.
#'
#' @description Calculates fit errors, AIC, BIC and
#' raw pseudo-R squares (1 - MSE/Var(y)) for a kinetic model.
#' @param kinetic_object a kinetic class object.
#' @param ... additional arguments passed to \code{\link{residuals.kinetic}}.
#' @return a data frame with the n number of complete observations, AIC, BIC,
#' R squared, deviance, mae (mean absolute error), mse (mean squared error) and
#' rmse (root mean squared error).

  get_stats <- function(kinetic_object, ...) {

    stopifnot(kinet::is_kinetic(kinetic_object))

    resid_tbl <- residuals(kinetic_object, ...)

    tibble::tibble(n_complete = nobs(kinetic_object),
                   aic = AIC(kinetic_object),
                   bic = BIC(kinetic_object),
                   raw_rsq = 1 - (mean(resid_tbl$.resid^2)/var(resid_tbl$.outcome)),
                   deviance = deviance(kinetic_object, REML = FALSE),
                   mae = mean(abs(resid_tbl$.resid)),
                   mse = mean(resid_tbl$.resid^2),
                   rmse = sqrt(mean(resid_tbl$.resid^2)))

  }

#' Test for normality of the model residuals.
#'
#' @description Tests the normality of the model resioduals with Shapiro-Wilk
#' test.
#' @inheritParams get_stats
#' @return a tibble with the test results: test statistic and p value.

  normality <- function(kinetic_object, ...) {

    stopifnot(kinet::is_kinetic(kinetic_object))

    resids <- residuals(kinetic_object)$.resid

    tst_results <- shapiro.test(resids)

    tibble::tibble(type = 'normality',
                   test = 'Shapiro-Wilk test',
                   stat_name = 'W',
                   stat_value = tst_results[['statistic']],
                   df1 = NA,
                   df2 = NA,
                   p_value = tst_results[['p.value']])

  }

# LRT ----

#' Likelihood ratio test for a kinetic class object.
#'
#' @description Performs a likelihood ratio test for a kinetic class object.
#' Model terms are compared in the 0 - n order; additionally, the full model is
#' compared with the NULL model (order = 'global').
#' @param kinetic_object a kinetic class object.
#' @param ... additional arguments passed to \code{\link[stats]{anova}}.
#' @return a data frame with the likelihood, AIC, BIC values and the test
#' statistics.
#' @details The effect size statistic lambda is calculated with the formula
#' 2 * (logLik(model) - logLik(NULL model))
#' @export

  lrt <- function(kinetic_object, ...) {

    if(!kinet::is_kinetic(kinetic_object)) {

      stop('A kinetic class object is required.', call. = FALSE)

    }

    ## nested models

    nested_models <- purrr::map(kinetic_object$order:0,
                                kinet::model_kinetic,
                                data = model.frame(kinetic_object),
                                response = kinetic_object$response,
                                time = kinetic_object$time,
                                ID = kinetic_object$ID,
                                family = kinetic_object$family)

    nested_models <- rlang::set_names(nested_models,
                                      as.character(kinetic_object$order:0))

    ## LRT, particular terms

    lrt_res <-  purrr::map(1:length(nested_models),
                           function(x) try(anova(nested_models[[x + 1]]$model,
                                                 nested_models[[x]]$model),
                                           silent = TRUE))

    lrt_res <- purrr::map(lrt_res,
                          function(x) if(any(class(x) == 'try-error')) NULL else x)

    lrt_res <- purrr::compact(lrt_res)

    ## LRT, full versus null model

    lrt_res$global <- anova(nested_models[[length(nested_models)]]$model,
                            nested_models[[1]]$model)

    ## formatting the output

    output <- purrr::map_dfr(lrt_res[-length(lrt_res)],
                             ~data.frame(.x[2, ]))

    output <- rbind(output,
                    data.frame(lrt_res[[length(lrt_res)]][1, ]))

    output <- rbind(output,
                    data.frame(lrt_res$global[2, ]))

    output <- dplyr::mutate(output,
                            order = c(names(nested_models), 'global'),
                            order = factor(order),
                            response = kinetic_object$response)

    output <- tibble::as_tibble(dplyr::arrange(output, order))

    ## calculating the lambdas

    dplyr::mutate(output,
                  lambda = 2*(logLik - output$logLik[1]))

  }

#' Likelihood ration test for a list of kinetic objects.
#'
#' @description Performs likelihood ratio tests for a list of kinetic class
#' objects, optionally in parallel.
#' @param kinetic_objects a list with kinetic class objects.
#' @param .parallel logical, should the model construction be done in parallel?
#' @param .errors logical, should modeling failures be returned?
#' @param ... extra arguments passed to \code{\link[stats]{anova}}.
#' @return a list of data frames with the likelihood, AIC, BIC values and
#' the test statistics.
#' @details The effect size statistic lambda is calculated with the formula
#' 2 * (logLik(model) - logLik(NULL model))
#' @export

  lrt_list <- function(kinetic_objects,
                       .parallel = FALSE,
                       .errors = FALSE, ...) {

    ## entry control

    if(!is.list(kinetic_objects)) {

      stop('The function requires a list of kinetic objects.',
           call. = FALSE)

    }

    class_check <- purrr::map_lgl(kinetic_objects, kinet::is_kinetic)

    if(any(!class_check)) {

      stop('Non-kinetic class object provided.',
           call. = FALSE)

    }

    ## Benchmarking

    start_time <- Sys.time()
    message(paste('LRT check for',
                  length(kinetic_objects),
                  'kinetic model objects'))
    on.exit(message(paste('Elapsed', Sys.time() - start_time)))

    ## Computation

    if(.parallel) {

      future::plan('multisession')

      lrt_result_lst <- furrr::future_map(kinetic_objects,
                                          purrr::safely(lrt),
                                          ...,
                                          .options = furrr::furrr_options(packages = c('dplyr',
                                                                                       'tibble',
                                                                                       'rlang',
                                                                                       'lme4',
                                                                                       'lmerTest')))

      future::plan('sequential')

    } else {

      lrt_result_lst <- purrr::map(kinetic_objects,
                                   purrr::safely(lrt), ...)

    }

    lrt_result_lst <- purrr::transpose(lrt_result_lst)

    lrt_result_lst <- purrr::map(lrt_result_lst, purrr::compact)

    if(.errors) {

      return(lrt_result_lst)

    } else {

      return(lrt_result_lst$result)

    }

  }


# END ----
