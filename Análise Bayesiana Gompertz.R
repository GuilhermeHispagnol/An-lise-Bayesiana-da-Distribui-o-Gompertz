# Pacotes
library(stats4)
library(evd)
library(coda)
library(ggplot2)
library(grid)
library(gridExtra)

# Parâmetros do MCMC
R <- 105000
burn_in <-5000  
jump <- 5
num_samples <- R - burn_in

# Dados observados
x <- c(0.31, 2.18, 2.55, 3.63, 2.44, 2.07, 3.06, 3.00, 2.63, 1.08,
       2.38, 2.45, 1.16, 0.95, 2.29, 2.02, 1.60, 2.12, 1.19, 2.31,
       1.58, 0.59, 3.82, 2.68, 2.53, 1.01, 2.51, 1.86, 2.50, 3.51,
       2.07, 1.01, 1.12, 2.02, 0.91, 2.09, 1.27, 1.96, 2.67, 1.75,
       2.32, 2.87, 2.08, 2.80, 2.35, 0.28, 1.51, 2.43, 0.98, 2.12)

# Hiperparâmetros
ni_alfa <- 0.04
ni_beta <- 0.24

# Log-posterior da Gompertz
log_posterior <- function(a, b) {
  if (a <= 0 || b <= 0) return(-Inf)
  termo <- sum(exp(b * x) - 1)
  n <- length(x)
  n * log(a) + b * sum(x) - (a / b) * termo +
    dgamma(a, 0.01, 0.01, log = TRUE) +
    dgamma(b, 0.01, 0.01, log = TRUE)
}

# Função MCMC otimizada (sem thinning interno)
run_chain_optimized <- function(init_alfa, init_beta, R, burn_in) {
  samples_alfa <- numeric(R - burn_in)
  samples_beta <- numeric(R - burn_in)
  current_alfa <- init_alfa
  current_beta <- init_beta
  accept_alfa <- 0
  accept_beta <- 0
  
  for (i in 1:R) {
    # Proposta alfa
    prop_alfa <- rgamma(1, shape = current_alfa^2 / ni_alfa^2,
                        scale = ni_alfa^2 / current_alfa)
    log_ratio_alfa <- log_posterior(prop_alfa, current_beta) -
      log_posterior(current_alfa, current_beta) +
      dgamma(current_alfa, shape = prop_alfa^2 / ni_alfa^2,
             scale = ni_alfa^2 / prop_alfa, log = TRUE) -
      dgamma(prop_alfa, shape = current_alfa^2 / ni_alfa^2,
             scale = ni_alfa^2 / current_alfa, log = TRUE)
    if (log(runif(1)) < log_ratio_alfa) {
      current_alfa <- prop_alfa
      if (i > burn_in) accept_alfa <- accept_alfa + 1
    }
    
    # Proposta beta
    prop_beta <- rgamma(1, shape = current_beta^2 / ni_beta^2,
                        scale = ni_beta^2 / current_beta)
    log_ratio_beta <- log_posterior(current_alfa, prop_beta) -
      log_posterior(current_alfa, current_beta) +
      dgamma(current_beta, shape = prop_beta^2 / ni_beta^2,
             scale = ni_beta^2 / prop_beta, log = TRUE) -
      dgamma(prop_beta, shape = current_beta^2 / ni_beta^2,
             scale = ni_beta^2 / current_beta, log = TRUE)
    if (log(runif(1)) < log_ratio_beta) {
      current_beta <- prop_beta
      if (i > burn_in) accept_beta <- accept_beta + 1
    }
    
    if (i > burn_in) {
      idx <- i - burn_in
      samples_alfa[idx] <- current_alfa
      samples_beta[idx] <- current_beta
    }
  }
  
  return(list(
    alfa = samples_alfa,
    beta = samples_beta,
    accept_alfa = accept_alfa / (R - burn_in),
    accept_beta = accept_beta / (R - burn_in)
  ))
}

# Rodar 3 cadeias
set.seed(123)
chain1 <- run_chain_optimized(0.05, 0.8, R, burn_in)
chain2 <- run_chain_optimized(0.1, 1.0, R, burn_in)
chain3 <- run_chain_optimized(0.15, 1.2, R, burn_in)

# =========================================================
# Função para aplicar thinning dinamicamente
# =========================================================

apply_thinning <- function(chain1, chain2, chain3, jump){
  
  thin_chain <- function(samples, jump){
    idx <- seq(1, length(samples), by = jump)
    samples[idx]
  }
  
  alfa <- c(
    thin_chain(chain1$alfa, jump),
    thin_chain(chain2$alfa, jump),
    thin_chain(chain3$alfa, jump)
  )
  
  beta <- c(
    thin_chain(chain1$beta, jump),
    thin_chain(chain2$beta, jump),
    thin_chain(chain3$beta, jump)
  )
  
  list(alfa = alfa, beta = beta)
}

# =========================================================
# Aplicar thinning inicial
# =========================================================
# jump = 10 # aqui para o thining ser alterado
thin_results <- apply_thinning(chain1, chain2, chain3, jump)

alfa_thin <- thin_results$alfa
beta_thin <- thin_results$beta
#########################

# Diagnósticos
mcmc_chains <- mcmc.list(
  mcmc(cbind(alfa = chain1$alfa, beta = chain1$beta)),
  mcmc(cbind(alfa = chain2$alfa, beta = chain2$beta)),
  mcmc(cbind(alfa = chain3$alfa, beta = chain3$beta))
)

cat("Taxas de aceitação:\n")
cat("Cadeia 1 - alfa:", chain1$accept_alfa, "beta:", chain1$accept_beta, "\n")
cat("Cadeia 2 - alfa:", chain2$accept_alfa, "beta:", chain2$accept_beta, "\n")
cat("Cadeia 3 - alfa:", chain3$accept_alfa, "beta:", chain3$accept_beta, "\n")

print(summary(mcmc_chains))
print(gelman.diag(mcmc_chains))
print(effectiveSize(mcmc_chains))

# Tabela resumo
ic_alfa <- quantile(alfa_thin, c(0.025, 0.975))
ic_beta <- quantile(beta_thin, c(0.025, 0.975))

summary_table <- data.frame(
  Parâmetro = c("alfa", "beta"),
  Média = c(mean(alfa_thin), mean(beta_thin)),
  `Desvio Padrão` = c(sd(alfa_thin), sd(beta_thin)),
  `IC 2.5%` = c(ic_alfa[1], ic_beta[1]),
  `IC 97.5%` = c(ic_alfa[2], ic_beta[2])
)
print(summary_table)

# --- Traceplots + Densidades ---
df_trace <- data.frame(
  Iteration = seq(1, length(chain1$alfa)) + burn_in,
  alfa = chain1$alfa,
  beta = chain1$beta
)

p1 <- ggplot(df_trace, aes(x = Iteration, y = alfa)) +
  geom_line(color = "blue") +
  labs(title = "Traceplot de alfa", x = "Iteração", y = expression(alpha)) +
  theme_minimal()

p2 <- ggplot(df_trace, aes(x = Iteration, y = beta)) +
  geom_line(color = "red") +
  labs(title = "Traceplot de beta", x = "Iteração", y = expression(beta)) +
  theme_minimal()

p3 <- ggplot(data.frame(alfa = alfa_thin), aes(x = alfa)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Posterior de alfa", x = expression(alpha), y = "Densidade") +
  theme_minimal()

p4 <- ggplot(data.frame(beta = beta_thin), aes(x = beta)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Posterior de beta", x = expression(beta), y = "Densidade") +
  theme_minimal()

grid.arrange(p1, p2, p3, p4, ncol = 2)

# --- ACF ---
par(mfrow = c(1, 2))
acf(alfa_thin, main = "ACF de alfa", col = "blue", lwd = 2)
acf(beta_thin, main = "ACF de beta", col = "red", lwd = 2)
par(mfrow = c(1, 1))

cat("Taxas de aceitação:\n")
cat("Cadeia 1 - alfa:", chain1$accept_alfa, "beta:", chain1$accept_beta, "\n")
cat("Cadeia 2 - alfa:", chain2$accept_alfa, "beta:", chain2$accept_beta, "\n")
cat("Cadeia 3 - alfa:", chain3$accept_alfa, "beta:", chain3$accept_beta, "\n")

cat("\nResumo das cadeias:\n")
print(summary(mcmc_chains))

cat("\nGelman-Rubin diagnostic:\n")
print(gelman.diag(mcmc_chains, autoburnin = FALSE, multivariate = TRUE))

cat("\nEffective Sample Size:\n")
print(effectiveSize(mcmc_chains))

cat("\nGeweke diagnostic:\n")
print(lapply(mcmc_chains, geweke.diag))

cat("\nHeidelberger-Welch diagnostic:\n")
print(lapply(mcmc_chains, heidel.diag))

cat("\nRaftery-Lewis diagnostic:\n")
print(lapply(mcmc_chains, raftery.diag))

#------------- gráficos de ajuste

# --- Curva Gompertz ajustada ---
alfa_mean <- mean(alfa_thin)
beta_mean <- mean(beta_thin)
fdp_gompertz <- function(x) {
  alfa_mean * exp(beta_mean * x - alfa_mean / beta_mean * (exp(beta_mean * x) - 1))
}

dfx <- data.frame(x = x)
p5 <- ggplot(dfx, aes(x = x)) +
  geom_histogram(
    aes(y = after_stat(density), fill = "Histograma"),
    bins = 10,
    color = "black",
    alpha = 0.4,
    linewidth = 0.3
  ) +
  geom_density(
    aes(color = "Densidade Empírica", linetype = "Densidade Empírica"),
    size = 1,
    key_glyph = "path"
  ) +
  stat_function(
    fun = fdp_gompertz,
    aes(color = "Gompertz Ajustada", linetype = "Gompertz Ajustada"),
    size = 1,
    key_glyph = "path"
  ) +
  scale_fill_manual(name = NULL, values = c("Histograma" = "lightgreen")) +
  scale_color_manual(name = NULL, values = c(
    "Densidade Empírica" = "blue",
    "Gompertz Ajustada" = "red"
  )) +
  scale_linetype_manual(name = NULL, values = c(
    "Densidade Empírica" = "solid",
    "Gompertz Ajustada" = "dashed"
  )) +
  guides(
    fill = guide_legend(
      override.aes = list(
        fill = "lightgreen",
        color = "lightgreen",
        alpha = 0.4,
        linewidth = 0.3
      ),
      order = 1
    ),
    color = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  ) +
  labs(
    title = "Histograma com Curva Ajustada",
    x = "Tempo",
    y = "Densidade"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.spacing.x = unit(0.5, 'cm'),
    legend.key.width = unit(1.5, 'cm'),
    legend.key = element_rect(fill = "white", color = NA),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

grid.arrange(p5, ncol = 1)

# --- Log-verossimilhança ponto a ponto ---
n_samples <- length(alfa_thin)
log_lik_matrix <- matrix(NA, nrow = n_samples, ncol = length(x))

for (i in 1:n_samples) {
  a_i <- alfa_thin[i]
  b_i <- beta_thin[i]
  log_lik_matrix[i, ] <- log(a_i) + b_i * x - (a_i / b_i) * (exp(b_i * x) - 1)
}

# --- WAIC ---
lppd <- sum(log(colMeans(exp(log_lik_matrix))))
p_waic <- sum(apply(log_lik_matrix, 2, var))
waic <- -2 * (lppd - p_waic)
cat("\nWAIC =", round(waic, 2), "\n")

# --- DIC ---
loglik_mean_params <- sum(log(alfa_mean) + beta_mean * x -
                            (alfa_mean / beta_mean) * (exp(beta_mean * x) - 1))
mean_deviance <- -2 * mean(rowSums(log_lik_matrix))
dev_mean_params <- -2 * loglik_mean_params
dic <- 2 * mean_deviance - dev_mean_params
cat("DIC =", round(dic, 2), "\n")

#########

# ----------------------------------------------------------------------
# Gráficos adicionais: verossimilhança, priori e posteriori por parâmetro
# ----------------------------------------------------------------------

# 1. Definir funções auxiliares
# Log-verossimilhança (sem as prioris)
log_likelihood <- function(a, b) {
  if (a <= 0 || b <= 0) return(-Inf)
  termo <- sum(exp(b * x) - 1)
  n <- length(x)
  n * log(a) + b * sum(x) - (a / b) * termo
}

# Priori para alfa (Gamma(0.01, 0.01) - shape e rate)
dprior_alfa <- function(a) dgamma(a, shape = 0.01, rate = 0.01)
# Priori para beta (Gamma(0.01, 0.01))
dprior_beta <- function(b) dgamma(b, shape = 0.01, rate = 0.01)

# Médias posteriores (já calculadas)
alfa_mean <- mean(alfa_thin)
beta_mean <- mean(beta_thin)

# 2. Grid de valores para cada parâmetro (baseado nos intervalos observados)
grid_alfa <- seq(0.001, max(alfa_thin) * 1.2, length.out = 200)
grid_beta <- seq(0.001, max(beta_thin) * 1.2, length.out = 200)

# 3. Calcular log-verossimilhança condicional para alfa (fixando beta na média)
loglik_alfa <- sapply(grid_alfa, function(a) log_likelihood(a, beta_mean))
# Para plotar na escala original (verossimilhança), exponenciamos e normalizamos
lik_alfa <- exp(loglik_alfa - max(loglik_alfa))  # normalizada para máximo = 1

# 4. Calcular log-verossimilhança condicional para beta (fixando alfa na média)
loglik_beta <- sapply(grid_beta, function(b) log_likelihood(alfa_mean, b))
lik_beta <- exp(loglik_beta - max(loglik_beta))

# 5. Data frames para ggplot
df_alfa <- data.frame(
  alfa = grid_alfa,
  Verossimilhança = lik_alfa,
  Priori = dprior_alfa(grid_alfa),
  Posteriori = NA  # será preenchido pela densidade estimada
)
df_beta <- data.frame(
  beta = grid_beta,
  Verossimilhança = lik_beta,
  Priori = dprior_beta(grid_beta),
  Posteriori = NA
)

# Estimar densidades posteriores
dens_alfa <- density(alfa_thin)
dens_beta <- density(beta_thin)

# Interpolar para os grids (para sobreposição nos mesmos pontos)
df_alfa$Posteriori <- approx(dens_alfa$x, dens_alfa$y, xout = grid_alfa, rule = 2)$y
df_beta$Posteriori <- approx(dens_beta$x, dens_beta$y, xout = grid_beta, rule = 2)$y

# Ajustar escalas: a verossimilhança está normalizada (0 a 1), as densidades estão em escalas diferentes.
# Para visualização conjunta, vamos normalizar também a priori e posteriori para terem máximo 1.
df_alfa$Priori_norm <- df_alfa$Priori / max(df_alfa$Priori, na.rm = TRUE)
df_alfa$Posteriori_norm <- df_alfa$Posteriori / max(df_alfa$Posteriori, na.rm = TRUE)
df_beta$Priori_norm <- df_beta$Priori / max(df_beta$Priori, na.rm = TRUE)
df_beta$Posteriori_norm <- df_beta$Posteriori / max(df_beta$Posteriori, na.rm = TRUE)

# 6. Gráfico para alfa
p_alfa <- ggplot(df_alfa, aes(x = alfa)) +
  geom_line(aes(y = Verossimilhança, color = "Verossimilhança"), size = 1) +
  geom_line(aes(y = Priori_norm, color = "Priori"), size = 1, linetype = "dashed") +
  geom_line(aes(y = Posteriori_norm, color = "Posteriori"), size = 1, linetype = "dotted") +
  scale_color_manual(
    name = NULL,
    values = c("Verossimilhança" = "darkgreen", "Priori" = "blue", "Posteriori" = "red")
  ) +
  labs(
    title = expression(alpha ~ " - Verossimilhança, Priori e Posteriori (normalizadas)"),
    x = expression(alpha),
    y = "Densidade / Verossimilhança (normalizada)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

# 7. Gráfico para beta
p_beta <- ggplot(df_beta, aes(x = beta)) +
  geom_line(aes(y = Verossimilhança, color = "Verossimilhança"), size = 1) +
  geom_line(aes(y = Priori_norm, color = "Priori"), size = 1, linetype = "dashed") +
  geom_line(aes(y = Posteriori_norm, color = "Posteriori"), size = 1, linetype = "dotted") +
  scale_color_manual(
    name = NULL,
    values = c("Verossimilhança" = "darkgreen", "Priori" = "blue", "Posteriori" = "red")
  ) +
  labs(
    title = expression(beta ~ " - Verossimilhança, Priori e Posteriori (normalizadas)"),
    x = expression(beta),
    y = "Densidade / Verossimilhança (normalizada)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

# 8. Exibir os gráficos lado a lado
#grid.arrange(p_alfa, p_beta, ncol = 2)

# Exibir os gráficos separadamente
print(p_alfa)
print(p_beta)