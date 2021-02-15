
writeTempModel <- function (

  method_name

) {

  if (method_name == "berry") {

    model_string <- "
      ## Bayesian Hierarchical Model: Berry

      model {

        ## Prior
        #mu_raw   ~ dnorm(0, 1)
        #mu       <- mean_mu + tau * mu_raw
        mu       ~ dnorm(mean_mu, precision_mu)

        tau      ~ dnorm(0, precision_tau)I(0.001, )
        tau_prec <- pow(tau, -2)

        for (j in 1:J) {

          ## Prior
          theta[j] ~ dnorm(mu, tau_prec)

          ## Likelihood
          logit(p[j]) <- theta[j] + logit(p_t[j])
          r[j] ~ dbin(p[j], n[j])

        }

      }"

  } else if (method_name == "berry_mix") {

    model_string <- "
      ## Bayesian Hierarchical Model: Berry with mixture distribution for mu

      model {

        ## Prior
                  # Arguments of dnormmix must have length > 1!
        mu       ~ dnormmix(mean_mu, precision_mu, pi_mu)

        tau      ~ dnorm(0, precision_tau)I(0.001, )
        tau_prec <- pow(tau, -2)

        for (j in 1:J) {

          ## Prior
          theta[j] ~ dnorm(mu, tau_prec)

          ## Likelihood
          logit(p[j]) <- theta[j] + logit(p_t[j])
          r[j] ~ dbin(p[j], n[j])

        }

      }"

  } else if (method_name == "exnex") {

    model_string <- "
      ## Bayesian Hierarchical Model: ExNex
      ## Model taken from Neuenschwander, Beat, et al.
      ##  'Robust exchangeability designs for early phase clinical trials with multiple strata.'
      ##  Pharmaceutical statistics 15.2 (2016): 123-134.
      ##  https://doi.org/10.1002/pst.1730

      model {

        # prior distributions for EX-parameters
        for (jj in 1:Nexch) {
          mu[jj] ~ dnorm(mu_mean[jj], mu_prec[jj])
          prior_tau_prec[jj] <- pow(tau_HN_scale[jj], -2)
          tau[jj] ~ dnorm(0, prior_tau_prec[jj])I(0.001, )
          prec_tau[jj] <- pow(tau[jj], -2)
        }

        # log-odds parameters under EX
        for (jj in 1:Nexch) {
          for (j in 1:Nstrata) {
            re[jj,j] ~ dnorm(0, prec_tau[jj])
            LogOdds[jj, j] <- mu[jj] + re[jj, j]
          }
        }

        # log-odds parameters under NEX
        for (j in 1:Nstrata) {
          LogOdds[Nmix, j] ~ dnorm(nex_mean[j], nex_prec[j]) ## Deviation from Appendix of Neuenschwander:
        }                                                    ## Different parameters for each
                                                             ## Nex component possible.

        # latent mixture indicators:
        # exch_index: categorial 1, ..., Nmix = Nexch + 1
        # exch: Nstrata x Nmix matrix of 0/1 elements
        for (j in 1:Nstrata) {
          exch_index[j] ~ dcat(pMix[1:Nmix])
          for (jj in 1:Nmix) {
            exch[j, jj] <- equals(exch_index[j], jj)
          }
        }

        # pick theta
        for (j in 1:Nstrata) {
          theta[j] <- LogOdds[exch_index[j], j]
        }

        # likelihood part
        for (i in 1:Nstrata) {
          logit(p[i]) <- theta[i]
          # p_success[i] <- step(p[i] - p_cut) ## Deviation from Appendix of Neuenschwander:
          r[i] ~ dbin(p[i], n[i])              ## No posterior response rate assessment within JAGS.
        }

      }"

  } else if (method_name == "exnex_adj") {

    model_string <- "
      ## Bayesian Hierarchical Model: ExNex with Adjustment
      ## Model adapted from Neuenschwander, Beat, et al.
      ##  'Robust exchangeability designs for early phase clinical trials with multiple strata.'
      ##  Pharmaceutical statistics 15.2 (2016): 123-134.
      ##  https://doi.org/10.1002/pst.1730

      model {

        # prior distributions for EX-parameters
        for (jj in 1:Nexch) {
          mu[jj] ~ dnorm(mu_mean[jj], mu_prec[jj])
          prior_tau_prec[jj] <- pow(tau_HN_scale[jj], -2)
          tau[jj] ~ dnorm(0, prior_tau_prec[jj])I(0.001, )
          prec_tau[jj] <- pow(tau[jj], -2)
        }

        # log-odds parameters under EX
        for (jj in 1:Nexch) {
          for (j in 1:Nstrata) {
            re[jj,j] ~ dnorm(0, prec_tau[jj])
            LogOdds[jj, j] <- mu[jj] + re[jj, j]
          }
        }

        # log-odds parameters under NEX
        for (j in 1:Nstrata) {
          LogOdds[Nmix, j] ~ dnorm(nex_mean[j], nex_prec[j])
        }

        # latent mixture indicators:
        # exch_index: categorial 1, ..., Nmix = Nexch + 1
        # exch: Nstrata x Nmix matrix of 0/1 elements
        for (j in 1:Nstrata) {
          exch_index[j] ~ dcat(pMix[1:Nmix])
          for (jj in 1:Nmix) {
            exch[j, jj] <- equals(exch_index[j], jj)
          }
        }

        # pick theta
        for (j in 1:Nstrata) {
          theta[j] <- LogOdds[exch_index[j], j]
        }

        # likelihood part
        for (i in 1:Nstrata) {
          logit(p[i]) <- theta[i] + logit(p_target[i])
          r[i] ~ dbin(p[i], n[i])
        }

      }"

  } else {

    stop ("method_name must be one of berry, berry_mix, exnex, exnex_adj")

  }

  temp_model <- getTempModel(model_string)

  return (temp_model)

}

getTempModel <- function (

  temp_model_string

) {

  temp_model       <- tempfile()
  write_connection <- file(description = temp_model, open = "w")
  cat(temp_model_string, file = write_connection)
  close(write_connection)

  return (temp_model)

}
