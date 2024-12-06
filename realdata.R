## Libraries
library(agop)
library(pracma)
library(expint)
library(hypergeo)
library(rootSolve)
library(ggplot2)
library(dplyr)
library(cubature)
library(plotly)
library(viridisLite)
library(htmlwidgets)
library(HellCor)
library(latex2exp)

## Source code
url <- "https://biostatisticien.eu/Qlm/"
source(paste(url, "qlm.R", sep = ""))

############################
# Application on real data #
############################

# 1) We load the data
#    they are organized in two files whose first
#    column contains the observed values in the
#    3D image and the other columns contain the
#    counts associated to each value (one column
#     per patient)

filebefore <- "df_before.csv" # data collected before treatment
fileafter <- "df_after.csv"   # data collected after treatment

# we read in the data (see Section 2.1 for the description)
df.before <- read.csv(paste(url, filebefore, sep = ""))
df.after <- read.csv(paste(url, fileafter, sep = ""))

# names (codes in fact, to guarantee anonymization) of patients 
patient.names <- substr(colnames(df.before)[-1], 1, 6)


# Figure 2.1
# plotting the two (pre- and post-treatment) empirical quantile curves for a given patient 
# Adjust the graphical parameters
cairo_pdf("fig2.1-right.pdf", width = 6.09, height = 6.33, pointsize = 10)
par(mar = c(3.4, 4.2, 1, 1))
patient.index <- 9
grid.p <- seq(from = 0, to = 1, length = 1000)
# Our qemp() function implements the first formula in Section 2.2.
plot(grid.p, qemp(p = grid.p, sample = df.before[, 1], counts = df.before[, patient.index + 1]), type = "l", ylim = range(df.before[, 1]), xlab = "", ylab = "", main = "", lty = 2, axes = FALSE)
points(grid.p, qemp(p = grid.p, sample = df.after[, 1], counts = df.after[, patient.index + 1]), type = "l", ylim = range(df.before[, 1]), lty = 3)
legend("topleft", lty = 2:3, legend = c("before treatment", "after treatment"), cex = 1.5)
# Add custom x-axis with specific mgp settings for x-axis
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext("HU", side = 2, line = 2.7, adj = 0.5, cex = 2)
mtext("p", side = 1, line = 2.2, adj = 0.5, cex = 2)
box()
dev.off()


# Transform raw data for each patient into a
# quantile function $q_{x,i}=q_{x,i,1}^\star$ in the 
# space of polynomials (d=1) of the standard 
# normal quantile function as described at the
# end of Section 2.2.
nb.patients <- length(patient.names)
qx <- qy <- vector("list", nb.patients)
rawsigmaqx <- rawsigmaqy <- rawmuqx <- rawmuqy <-
    sigmaqx <- sigmaqy <- muqx <- muqy <-
        rep(NA, nb.patients)
values.before <- df.before[, 1]
values.after <- df.after[, 1]
weighted.mean <- function(x, w) sum(na.omit(x * w )) / sum(na.omit(w))
weighted.sd <- function(x, y, w) {
    sqrt(weighted.mean(x ^ 2, w) - y ^ 2)
}
for (i in 1:nb.patients) {
    counts.before <- df.before[, i + 1]
    counts.after <- df.after[, i + 1]
    qx[[i]] <- qopt(sample = values.before, counts = counts.before, d = 1)
    qy[[i]] <- qopt(sample = values.after, counts = counts.after, d = 1)
    muqx[i] <- get("a0", envir = environment(qx[[i]]))
    sigmaqx[i] <- get("a1", envir = environment(qx[[i]]))
    muqy[i] <- get("a0", envir = environment(qy[[i]]))
    sigmaqy[i] <- get("a1", envir = environment(qy[[i]]))
    rawmuqx[i] <- weighted.mean(values.before, w = counts.before)
    rawmuqy[i] <- weighted.mean(values.after, w =counts.after)
    rawsigmaqx[i] <- weighted.sd(values.before, rawmuqx[i], counts.before)
    rawsigmaqy[i] <- weighted.sd(values.after, rawmuqy[i], counts.after)
}

# This graph is not presented in the paper but we believe it is interesting 
# to plot the empirical quantile curve for a given patient, as well as its 
# associated $q_{x,i}$ version (before treatment)
patient.index <- 9
grid.p <- seq(from = 0, to = 1, length = 1000)
plot(grid.p, qemp(p = grid.p, sample = df.before[, 1], counts = df.before[, patient.index + 1]), type = "l", ylim = range(df.before[, 1]), xlab = "p", ylab = "HU", main = paste("Patient ", patient.names[patient.index], " (before treatment)"), lty = 2)
points(grid.p, qx[[patient.index]](grid.p), type = "l", lty = 1)
legend("topleft", lty = 2:1, legend = c("Empirical", expression(q[opt])))

# Figure 3.1, at the end of Section 3.1
cairo_pdf("fig3.1.pdf", width = 10, height = 7)
mu <- 0; sigma <- 1; beta <- 2; delta <- 3
set.seed(1) # For reproduciblity purposes.
n <- 5 # Sample size.
a0 <- rnorm(n, mean = mu, sd = sigma)
a1 <- delta + rexp(n, rate = 1 / beta)
pvec <- seq(from = 0, to = 1, length = 1000)
res <- Vectorize(FUN = function(p) qnorm(p, mean = a0, sd = a1), 
        vectorize.args = "p")(pvec)
# Adjust the graphical parameters
par(mar = c(3.4, 4.2, 1, 1))
matplot(pvec, t(res), type = "l", xlab = "", ylab = "", lty = 2, col = "black", cex.axis = 2, cex.lab = 2, axes = FALSE)
# Add custom x-axis with specific mgp settings for x-axis
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext("q(p)", side = 2, line = 2.7, adj = 0.5, cex = 2)
mtext("p", side = 1, line = 2.2, adj = 0.5, cex = 2)
box()
dev.off()


# Figure 3.2 in Section 3.2
# Left part
cairo_pdf("fig3.2-left.pdf", width = 6.09, height = 6.33, pointsize = 10)
par(mar = c(3.4, 6, 2, 0.3), mfrow = c(2, 1), las = 1)
beta0 <- 5
beta1 <- 1
beta2 <- 1
muqx <- 1
sigmaqx <- 2
muqy <- beta0 + beta1 * muqx
sigmaqy <- beta2 * sigmaqx
curve(qnorm(x, mean = muqx, sd = sigmaqx), xlim = c(0, 1), ylim = c(-10, 15), xlab = "", ylab = "", axes = FALSE, main = "Quantile curves", cex.main = 2)
curve(qnorm(x, mean = muqy, sd = sigmaqy), add = TRUE, lty = 2)
legend1 <- bquote(mu[q[y]] == .(muqy))
legend2 <- bquote(sigma[q[y]] == .(sigmaqy))
legend3 <- bquote(mu[q[x]] == .(muqx))
legend4 <- bquote(sigma[q[x]] == .(sigmaqx))
legend("topleft", legend = c(as.expression(legend1), as.expression(legend3)),
       lty = c(2, 1), cex = 1.8, bty = "n", horiz = TRUE)
legend("topleft", legend = c(as.expression(legend2), as.expression(legend4)), lty = c(2, 1), cex = 1.8, bty = "n", horiz = TRUE, inset = c(0, 0.12), col = "white")
# Add custom x-axis with specific mgp settings for x-axis
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext("q(p)", side = 2, line = 4, adj = 0.5, cex = 2)
mtext("p", side = 1, line = 2.2, adj = 0.5, cex = 2)
box()
curve(dnorm(x, mean = muqx, sd = sigmaqx), xlim = c(-10, 15), ylim = c(0, 0.20), xlab = "", ylab = "", axes = FALSE, main = "Density curves", cex.main = 2)
curve(dnorm(x, mean = muqy, sd = sigmaqy), add = TRUE, lty = 2)
# Add custom x-axis with specific mgp settings for x-axis
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext("Density", side = 2, line = 4.4, adj = 0.5, cex = 2, las = 0)
mtext("Value", side = 1, line = 2.4, adj = 0.5, cex = 2)
box()
dev.off()
# Right part
cairo_pdf("fig3.2-right.pdf", width = 6.09, height = 6.33, pointsize = 10)
par(mar = c(3.4, 6, 2, 0.3), mfrow = c(2, 1), las = 1)
beta0 <- 0
beta1 <- 1
beta2 <- 3
muqx <- 1
sigmaqx <- 2
muqy <- beta0 + beta1 * muqx
sigmaqy <- beta2 * sigmaqx
curve(qnorm(x, mean = muqx, sd = sigmaqx), xlim = c(0, 1), ylim = c(-10, 15), xlab = "", ylab = "", axes = FALSE, main = "Quantile curves", cex.main = 2)
curve(qnorm(x, mean = muqy, sd = sigmaqy), add = TRUE, lty = 2)
legend1 <- bquote(mu[q[y]] == .(muqy))
legend2 <- bquote(sigma[q[y]] == .(sigmaqy))
legend3 <- bquote(mu[q[x]] == .(muqx))
legend4 <- bquote(sigma[q[x]] == .(sigmaqx))
legend("topleft", legend = c(as.expression(legend1), as.expression(legend3)),
       lty = c(2, 1), cex = 1.8, bty = "n", horiz = TRUE)
legend("topleft", legend = c(as.expression(legend2), as.expression(legend4)), lty = c(2, 1), cex = 1.8, bty = "n", horiz = TRUE, inset = c(0, 0.12), col = "white")
# Add custom x-axis with specific mgp settings for x-axis
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext("q(p)", side = 2, line = 4, adj = 0.5, cex = 2)
mtext("p", side = 1, line = 2.2, adj = 0.5, cex = 2)
box()
curve(dnorm(x, mean = muqx, sd = sigmaqx), xlim = c(-10, 15), ylim = c(0, 0.20), xlab = "", ylab = "", axes = FALSE, main = "Density curves", cex.main = 2)
curve(dnorm(x, mean = muqy, sd = sigmaqy), add = TRUE, lty = 2)
# Add custom x-axis with specific mgp settings for x-axis
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext("Density", side = 2, line = 4.4, adj = 0.5, cex = 2, las = 0)
mtext("Value", side = 1, line = 2.4, adj = 0.5, cex = 2)
box()
dev.off()


# We fit the model (3.3) and showcase our new 
# methods coef.qlm(), confint.qlm() and summary.qlm()
res.qlm <- qlm(qy ~ qx, data = list(qx = qx, qy = qy))
coef(res.qlm)
confint(res.qlm)
summary(res.qlm)
# We also show how to extract various other quantities
# used in Proposition 3.2
res.qlm$muqybar
res.qlm$sigmaqybar
res.qlm$muqxbar
res.qlm$sigmaqxbar
res.qlm$w

# Numerical verification of the L_{\alpha} values provided after Proposition 3.7
# This is also needed for Figure 3.3 below
integrand <- function(x, Lalpha, n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w, betahat, beta2hat, sigmaqxnp1, sigmaqxbar) {
    s <- x[1]
    t <- x[2]    
    val <- fQhat(s, t, n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w, betahat, beta2hat, sigmaqxnp1, sigmaqxbar)

    if (val >= Lalpha) return(val) else return(0)
}
muqxnp1 <- -750
sigmaqxnp1 <- 120
n <- res.qlm$n
beta0hat <- res.qlm$coefficients["beta0hat"]
beta1hat <- res.qlm$coefficients["beta1hat"]
sigma2hat <- res.qlm$coefficients["sigma2hat"]
muqxbar <- res.qlm$muqxbar
w <- res.qlm$w
betahat <- res.qlm$coefficients["betahat"]
beta2hat <- res.qlm$coefficients["beta2hat"]
sigmaqxbar <- res.qlm$sigmaqxbar
#Lalpha <- 0.000033 # for alpha = 0.01
Lalpha <- 0.000164 # for alpha = 0.05
#Lalpha <- 0.000328 # for alpha = 0.10
objective_function <- function(par, n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w, betahat, beta2hat, sigmaqxnp1, sigmaqxbar) {
  s <- par[1]
  t <- par[2]
  fQhat(s, t, n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w, betahat, beta2hat, sigmaqxnp1, sigmaqxbar)
}
# Optimization over the given bounds to find the maximum value of the density function
result <- optim(par = c(-750, 130), fn = objective_function, n = n, beta0hat = beta0hat, beta1hat = beta1hat, sigma2hat = sigma2hat, muqxnp1 = muqxnp1, muqxbar = muqxbar, w = w, betahat = betahat, beta2hat = beta2hat, sigmaqxnp1 = sigmaqxnp1, sigmaqxbar = sigmaqxbar, method = "L-BFGS-B",
                lower = c(-770, 90), upper = c(-730, 160), control = list(fnscale = -1))
max_value <- result$value
#lowerlimit <- c(-Inf, beta2hat * sigmaqxnp1)
lowerlimit <- c(-780, beta2hat * sigmaqxnp1)
#upperlimit <- c(Inf, Inf)
upperlimit <- c(-717.5, Inf)
# Gives a value equal (i.e., numerically close) to 1 - alpha 
cubature::adaptIntegrate(f = integrand, lowerLimit = lowerlimit, upperLimit = upperlimit, Lalpha = Lalpha, n = n, beta0hat = beta0hat, beta1hat = beta1hat, sigma2hat = sigma2hat, muqxnp1 = muqxnp1, muqxbar = muqxbar, w = w, betahat = betahat, beta2hat = beta2hat, sigmaqxnp1 = sigmaqxnp1, sigmaqxbar = sigmaqxbar, tol = 1e-16, maxEval = 10000)$integral
# Note that one can also use the following:
# alpha <- 0.05
# optim(par = 0.00016, fn = function(Lalpha) ((1 - alpha) - cubature::adaptIntegrate(f = integrand, lowerLimit = lowerlimit, upperLimit = upperlimit, Lalpha = Lalpha, n = n, beta0hat = beta0hat, beta1hat = beta1hat, sigma2hat = sigma2hat, muqxnp1 = muqxnp1, muqxbar = muqxbar, w = w, betahat = betahat, beta2hat = beta2hat, sigmaqxnp1 = sigmaqxnp1, sigmaqxbar = sigmaqxbar, tol = 1e-16, maxEval = 10000)$integral) ^ 2, lower = 0.0001, upper = 0.0002, method = "L-BFGS-B")


# Figure 3.3
library("ggplot2")
library("dplyr")
library("latex2exp")  # Load the latex2exp package for LaTeX rendering

# Define font sizes and styles
tickSize1 <- 18
tickSize2 <- 13
legendSize1 <- 15
legendKeySize <- 0.5
legendSize2 <- 15
fontSize1 <- 20
fontSize2 <- 15

# Set up the s and t grid for the plot
s_vals <- seq(-770, -727.5, length.out = 100)
t_vals <- seq(90, as.numeric(beta2hat * sigmaqxnp1 * 2 + 10), length.out = 100)
grid <- expand.grid(s = s_vals, t = t_vals)

# Calculate fQhat[s,t] for the grid points
grid$z <- apply(grid, 1, function(row) {
  fQhat(row[1], row[2], n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w, betahat, beta2hat, sigmaqxnp1, sigmaqxbar)
})

# Define breaks to include 0.000 at the bottom of the legend
legend_breaks <- c(0.000, 0.001, 0.002, 0.003)

# Contour plot
contour_plot <- ggplot(grid, aes(x = s, y = t, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", breaks = legend_breaks, limits = c(0, max_value)) +
  guides(fill = guide_colorbar(
    title = "Density",
    title.theme = element_text(size = legendSize1),
    label.theme = element_text(size = legendSize1),
    barwidth = unit(0.8, "lines"),
    barheight = unit(8, "lines")
  )) +
  labs(
    x = TeX(r'($\mu_{\widehat{Q}_{Y,n+1}}$)'),  
    y = TeX(r'($\sigma_{\widehat{Q}_{Y,n+1}}$)')
  ) +
  scale_x_continuous(breaks = seq(min(grid$s), max(grid$s), by = 10), expand = c(0, 0)) +  
  scale_y_continuous(breaks = seq(min(grid$t), max(grid$t), by = 10), expand = c(0, 0)) +  
  theme_minimal(base_size = fontSize1) +
  theme(
    axis.title.x = element_text(angle = 0, vjust = 2, hjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = -0.5, margin = margin(t = 0, r = -18, b = 0, l = 0)), 
    axis.text.x = element_text(size = tickSize1, margin = margin(0)),  
    axis.text.y = element_text(size = tickSize1, margin = margin(0)),  
    axis.ticks.length = unit(0.1, "cm"),  
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = legendSize1),
    legend.key.size = unit(legendKeySize, "cm"),
    plot.margin = margin(1, -8, -9, 2),  # Minimal plot margin
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, -10),  # Bring legend closer
    legend.box.spacing = unit(0.3, "cm"),  # Remove additional horizontal space
    legend.spacing.x = unit(0, "cm"),
    plot.title = element_text(size = 25, face = "bold"),
    panel.grid = element_blank()
  )

# Save contour plot as PDF
cairo_pdf("fig3.3.pdf", width = 8.09, height = 6.33, pointsize = 10)
print(contour_plot)
dev.off()

# Figure 3.3 (3d plot)
library("plotly")
library(viridisLite)
# Generate magma color scale from viridisLite package
magma_colors <- viridisLite::viridis(256, option = "magma")
# Create the colorscale as a list of evenly spaced stops (0 to 1)
colorscale <- list()
for (i in seq_along(magma_colors)) {
  colorscale[[i]] <- list((i - 1) / (length(magma_colors) - 1), magma_colors[i])
}
library("htmlwidgets")
# Interactive 3D plotly widget
s_vals_3d <- seq(-770, -727.5, length.out = 100)
t_vals_3d <- seq(90, as.numeric(beta2hat * sigmaqxnp1 * 2 + 10), length.out = 100)
z_vals_3d <- outer(s_vals_3d, t_vals_3d, Vectorize(function(s, t) fQhat(s, t, n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w, betahat, beta2hat, sigmaqxnp1, sigmaqxbar)))
plot_3d <- plot_ly(x = s_vals_3d,
                   y = t_vals_3d,
                   z = z_vals_3d,
                   type = "surface",
                   colorscale = colorscale) %>%
  layout(scene = list(xaxis = list(title = expression(hat(mu)[Q[Y, n+1]]), titlefont = list(size = fontSize2),
                                   tickfont = list(size = tickSize2)),
                      yaxis = list(title = expression(hat(sigma)[Q[Y, n+1]]), titlefont = list(size = fontSize2),
                                   tickfont = list(size = tickSize2)),
                      zaxis = list(title = 'Density', titlefont = list(size = fontSize2))),
         legend = list(font = list(size = legendSize2)))
saveWidget(plot_3d, file = "fig-3.3.html", selfcontained = TRUE)


# Figure 3.4
sigmahat <- as.numeric(sqrt(res.qlm$coefficients["sigma2hat"]))
ii <- 10
muqxi <- -767.5906
sigmaqxi <- 121.92861
nu2hati <- betahat * sigmaqxi / (n * (n - 1) * sigmaqxbar)
nu1hati <- betahat + nu2hati
factor1i <- (1 - 1/n - (muqxi - muqxbar)^2 / (n * w))
factor2i <- sigmaqxi / (n * sigmaqxbar) 
# Set plot styles
tickSize1 <- 18
tickSize2 <- 13
legendSize1 <- 15
legendKeySize <- 0.5
legendSize2 <- 15
fontSize1 <- 20
fontSize2 <- 15
# Density Plot
library(ggplot2)
library(viridis)

density_plot <- ggplot(data = expand.grid(s = seq(-120, 120, length.out = 100), t = seq(0, 150, length.out = 100)), aes(x = s, y = t)) +
  geom_raster(aes(fill = fEhat(s, t, n, sigmahat, factor1i, factor2i, nu1hati, nu2hati))) +
  scale_fill_viridis(option = "magma", direction = 1, discrete = FALSE) + 
  guides(fill = guide_colorbar(
    title = "Density",
    title.theme = element_text(size = legendSize1),
    label.theme = element_text(size = legendSize1),
    barwidth = unit(0.8, "lines"),
    barheight = unit(8, "lines")
  )) +
  labs(
    x = TeX(r'($mu_{\widehat{E}_i}$)'), 
    y = TeX(r'($\sigma_{\widehat{E}_i}$)')
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = fontSize1) +
  theme(
    axis.title.x = element_text(angle = 0, vjust = 2, hjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = -0.5, margin = margin(t = 0, r = -7, b = 0, l = 0)),
    axis.text.x = element_text(size = tickSize1, margin = margin(0)),  
    axis.text.y = element_text(size = tickSize1, margin = margin(0)),  
    axis.ticks.length = unit(0.1, "cm"),  
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = legendSize1),
    legend.key.size = unit(legendKeySize, "cm"),
    plot.margin = margin(4, -8, -9, 2),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, -10),  # Increase space between plot and legend
    legend.box.spacing = unit(0.3, "cm"),       # Further increase horizontal space
    legend.spacing.x = unit(0, "cm"),
    plot.title = element_text(size = 25, face = "bold"),
    panel.grid = element_blank()
  )

# Save density plot as PDF
cairo_pdf("fig3.4.pdf", width = 8.09, height = 6.33, pointsize = 10)
print(density_plot)
dev.off()

# Figure 3.4 (3d plot)
library(plotly)
library(viridisLite)
# Generate magma color scale from viridisLite package
magma_colors <- viridisLite::viridis(256, option = "magma")
# Create the matrix of z values using fEhat
s_vals <- seq(-120, 120, length.out = 100)
t_vals <- seq(0, 150, length.out = 100)
z_vals <- outer(s_vals, t_vals, fEhat, n = n, sigmahat = sigmahat, factor1i = factor1i, factor2i = factor2i, nu1hati = nu1hati, nu2hati = nu2hati)
# Normalize original z values for color scaling
z_min <- min(z_vals, na.rm = TRUE)
z_max <- max(z_vals, na.rm = TRUE)
# Create the colorscale as a list of evenly spaced stops (0 to 1)
colorscale <- list()
for (i in seq_along(magma_colors)) {
  colorscale[[i]] <- list((i - 1) / (length(magma_colors) - 1), magma_colors[i])
}
# Define tick values manually for colorbar (e.g., 6 ticks) to avoid scientific notation
tick_vals <- seq(z_min, z_max, length.out = 6)  # Evenly spaced tick values from min to max
tick_text <- sprintf("%.6f", tick_vals)         # Format tick values as strings with 6 decimal places
# Create 3D plot with original z-values using plot_ly
plot3d <- plot_ly(
  x = t_vals,  # x-axis values (sigma)
  y = s_vals,  # y-axis values (mu)
  z = z_vals,  # Use original z-values (no rescaling)
  type = "surface", 
  colorscale = colorscale,  # Custom colorscale from viridis magma
  surfacecolor = z_vals,    # Ensuring color scale is applied to z values
  cmin = z_min,             # Set the minimum color value (matches the range of z_vals)
  cmax = z_max              # Set the maximum color value (matches the range of z_vals)
) %>%
    layout(
      scene = list(
        xaxis = list(
          title = list(text = '$\\hat{\\sigma}_{E_i}$', standoff = 30),  # Proper LaTeX for sigma with hat and subscripts
          titlefont = list(size = 12)
        ),
        yaxis = list(
          title = list(text = '$\\hat{\\mu}_{E_i}$', standoff = 30),   # Proper LaTeX for mu with hat and subscripts
          titlefont = list(size = 12)
        ),
        zaxis = list(
          title = "Density",                # Correct label for z-axis
          titlefont = list(size = 12),
          tickformat = ".6f",               # Format for proper numerical representation with 6 decimal places
          tickfont = list(size = 10),       # Ensure proper font size for ticks
          range = c(z_min, z_max)           # Ensure the correct range is applied
        )
      ),
      colorbar = list(
        title = "Density",                  # Proper title for legend
        tickvals = tick_vals,               # Use manual tick values
        ticktext = tick_text,               # Use string format to force exact decimal display
        tickfont = list(size = 10)          # Ensure proper font size for the colorbar ticks
      )
    ) %>%
    config(mathjax = 'cdn')  # Enable MathJax rendering using CDN for LaTeX labels
# Display the 3D plot
#plot3d
saveWidget(plot3d, file = "fig-3.4.html", selfcontained = TRUE)


# Table 4.1
table4.1 <- cbind(res.qlm$muqx, res.qlm$sigmaqx, res.qlm$muqy, res.qlm$sigmaqy) 
rownames(table4.1) <- 1:nrow(table4.1)
colnames(table4.1) <- c("mu_{q_{x,i}}", "sigma_{q_{x,i}}", "mu_{q_{y,i}}", "sigma_{q_{y,i}}")
head(table4.1)


# Table 4.2 in Section 4
resids <- residuals(res.qlm)
n <- length(resids)
muresid <- sigmaresid <- rep(NA, n)
for (i in 1:n) {
  muresid[i] <- get("muresid", envir = environment(resids[[i]]))
  sigmaresid[i] <- get("sdresid", envir = environment(resids[[i]]))
}
table4.2 <- cbind(muresid, sigmaresid)
rownames(table4.2) <- 1:nrow(table4.1)
colnames(table4.2) <- c("mu_{hat{epsilon}_i}", "sigma_{hat{epsilon}_i}")
round(head(table4.2), 2)

set.seed(1); HellCor::HellCor(muresid, sigmaresid, pval.comp = TRUE)$p.value


# Figure 4.1 in Section 4
cairo_pdf("fig4.1.pdf", width = 10, height = 7)
par(mar = c(3.4, 6, 2, 0.3), las = 1)
plot(table4.2, pch = 20, col = "blue", xlab = "", ylab = "", cex.lab = 2, axes = FALSE)
text(x = table4.2[, 1], y = table4.2[, 2], labels = rownames(table4.2), adj = 1.2, cex = 1.5)
abline(v = 0)
abline(h = coef(res.qlm)["betahat"], lty = 2)
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext(expression(sigma[hat(epsilon)[i]]), side = 2, line = 4, adj = 0.5, cex = 2)
mtext(expression(mu[hat(epsilon)[i]]), side = 1, line = 2.2, adj = 0.5, cex = 2)
text(-135, 54, expression(hat(beta)), cex = 2)
box()
dev.off()


# Table 4.3 in Section 4
factor1 <- 1 - (1 / res.qlm$n) - (res.qlm$muqx - mean(res.qlm$muqx)) ^ 2 / (n * res.qlm$w)
factor2 <- res.qlm$sigmaqx / (res.qlm$n * mean(res.qlm$sigmaqx))
sigmahat <- as.numeric(sqrt(res.qlm$coefficients["sigma2hat"]))
betahat <- as.numeric(res.qlm$coefficients["betahat"])
nu2hat <- betahat * factor2 / (res.qlm$n - 1)
nu1hat <- betahat + nu2hat
res.qlm.pvals <- rep(NA, res.qlm$n)
for (i in 1:res.qlm$n) {
    factor1i <- factor1[i]
    factor2i <- factor2[i]
    nu1hati <- nu1hat[i]
    nu2hati <- nu2hat[i]
    muresidi <- muresid[i]
    sigmaresidi <- sigmaresid[i]
    res.qlm.pvals[i] <- compute.residual.pvalue(muresidi, sigmaresidi, n, sigmahat, factor1i, factor2i, nu1hati, nu2hati)
}
table4.3 <- round(matrix(res.qlm.pvals, ncol = 4), 4)


# Figure 4.2 in Section 4
cairo_pdf("fig4.2.pdf", width = 16, height = 7)
par(mar = c(3.2, 5, 1.5, 0.5), xaxs = "i", family = "Arial Narrow")
# Compute the bar heights
bar_heights <- -log(res.qlm.pvals, 10)
num_bars <- length(bar_heights)
# Create the barplot with specified width and space
x_pos <- barplot(
  height = bar_heights,
  width = 0.8,
  space = 0.25,
  col = "skyblue",
  ylim = c(0, 2.5),
  xlim = c(0, num_bars),
  axes = FALSE,
  names.arg = rep("", num_bars)
)
# Add the red dashed horizontal line
abline(h = 2, lty = 2, col = "red", lwd = 2)
# Add y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2.5)
# Rotate x-axis labels at 90 degrees and move them closer to ticks
axis(1, at = x_pos, labels = FALSE)
text(
  x = x_pos,
  y = par("usr")[3] - 0.06,
  labels = 1:num_bars,
  srt = 90,
  adj = 1,
  xpd = TRUE,
  cex = 2.3,
  family = "Arial Narrow"
)
# Add y-axis label closer to the axis
mtext(expression(-Log[10]("p-value")), side = 2, line = 2.5, adj = 0.5, cex = 2.5)
# Draw box around the plot
box()
dev.off()

# Figure 4.3 in Section 4
# Left part
cairo_pdf("fig4.3-left.pdf", width = 10, height = 7)
# Compute the average quantile functions
threshold <- -716
aveqx <- function(x) qnorm(x, mean = mean(res.qlm$muqx[res.qlm$muqx < -716]), sd = mean(res.qlm$sigmaqx[res.qlm$muqx < threshold]))
aveqy <- function(x) qnorm(x, mean = mean(res.qlm$muqy[res.qlm$muqx < -716]), sd = mean(res.qlm$sigmaqy[res.qlm$muqx < threshold]))
# Adjust the graphical parameters
par(mar = c(3.4, 7, 3, 1), las = 1)
curve(aveqx, lty = 1, xlim = c(0, 1), ylim = c(-1023, -200), n = 1000, axes = FALSE, cex.axis = 2, cex.lab = 2, xlab = "", ylab = "", main = "Average quantile functions", cex.main = 2)
curve(aveqy, lty = 4, add = TRUE, n = 1000)
muqyfitted <- sigmaqyfitted <- 0
for (i in which(muqx < threshold)) {
    muqyfitted <- muqyfitted + get("muqyfitted", envir = environment(res.qlm$fitted[[i]]))
    sigmaqyfitted <- sigmaqyfitted + get("sigmaqyfitted", envir = environment(res.qlm$fitted[[i]]))
}
muqyfitted <- muqyfitted / sum(res.qlm$muqx < threshold)
sigmaqyfitted <- sigmaqyfitted / sum(res.qlm$muqx < threshold)
aveqyfitted <- function(x) qnorm(x, mean = muqyfitted, sd = sigmaqyfitted)
curve(aveqyfitted, lty = 5, add = TRUE, n = 1000)
legend("topleft", lty = c(4, 5, 1), legend = c("after treatment", "fitted quantile", "before treatment"), cex = 2)
# Add custom x-axis with specific mgp settings for x-axis
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext("HU", side = 2, line = 5.7, adj = 0.5, cex = 2)
mtext("p", side = 1, line = 2.2, adj = 0.5, cex = 2)
box()
dev.off()
# Right part
cairo_pdf("fig4.3-right.pdf", width = 10, height = 7)
diff <- function(x) aveqy(x) - aveqx(x)
par(mar = c(3.4, 5.5, 3, 1.2), las = 1, xaxs = "i", yaxs = "i")
curve(diff, lty = 1, xlim = c(0, 1), ylim = c(-15, 35), n = 1000, axes = FALSE, cex.axis = 2, cex.lab = 2, xlab = "", ylab = "", main = "Average difference after/before treatment", cex.main = 2)
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2) # Adjust second value to move y-axis labels
mtext("HU", side = 2, line = 4.2, adj = 0.5, cex = 2)
mtext("p", side = 1, line = 2.2, adj = 0.5, cex = 2)
q10 <- uniroot(function(x) diff(x) - 10, lower = 10^-10, upper = 1-10^-10)$root
q20 <- uniroot(function(x) diff(x) - 20, lower = 10^-10, upper = 1-10^-10)$root
segments(x0 = q10, y0 = -15, x1 = q10, y1 = 10, lty = 2)
segments(x0 = q20, y0 = -15, x1 = q20, y1 = 20, lty = 2)
segments(x0 = 0, y0 = 10, x1 = q10, y1 = 10, lty = 2)
segments(x0 = 0, y0 = 20, x1 = q20, y1 = 20, lty = 2)
mtext(round(q10, 2), side = 1, line = -1.1, adj = q10 - 0.04, cex = 2)
mtext(round(q20, 2), side = 1, line = -1.1, adj = q20 - 0.01, cex = 2)
clip(0, 1, -15, 35)
box()
dev.off()


# Figure 4.4
patient.index <- 13
cairo_pdf("fig4.4-top.pdf", width = 8.09, height = 6.33, pointsize = 10)
par(mar = c(4, 8, 3, 1))
#par(mar = c(3.8, 7.2, 2, 1))
curve(qx[[patient.index]](x), axes = FALSE, main = "Quantile curves", xlab = "", ylab = "", lty = 1, ylim = c(-1300, -400), n = 10000, cex.main = 2)
curve(qy[[patient.index]](x), lty = 2, add = TRUE, n = 10000)
curve(res.qlm$fitted[[patient.index]](x), add = TRUE, lty = 3, n = 10000)
legend("topleft", lty = c(3, 2, 1), legend = c("fitted", "after treatment", "before treatment"), cex = 2)
grid(nx = 5, ny = 5, col = "lightgray", lty = "dotted")
# Add custom x-axis with specific mgp settings for x-axis
mtext("HU", side = 2, line = 5, adj = 0.5, cex = 2)
mtext("p", side = 1, line = 2.2, adj = 0.5, cex = 2)
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2, at = seq(from = 0, to = 1, by = 0.1)) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 0.7, 0), cex.axis = 2, las = 1, at = seq(from = -1300, to = -400, by = 100)) # Adjust second value to move y-axis labels
dev.off()

cairo_pdf("fig4.4-bottom.pdf", width = 8.09, height = 6.33, pointsize = 10)
par(mar = c(4, 8, 3, 1))
mu1 <- get("a0", envir = environment(qx[[patient.index]]))
sd1 <- get("a1", envir = environment(qx[[patient.index]]))
curve(dnorm(x, mean = mu1, sd = sd1), xlim = c(-1300, -400), ylim = c(0, 0.0040), ylab = "", xlab = "", lty = 1, axes = FALSE, main = "Density curves", cex.main = 2)
mu2 <- get("a0", envir = environment(qy[[patient.index]]))
sd2 <- get("a1", envir = environment(qy[[patient.index]]))
curve(dnorm(x, mean = mu2, sd = sd2), lty = 2, add = TRUE)
mu3 <- get("muqyfitted", envir = environment(res.qlm$fitted[[patient.index]]))
sd3 <- get("sigmaqyfitted", envir = environment(res.qlm$fitted[[patient.index]]))
curve(dnorm(x, mean = mu3, sd = sd3), lty = 3, add = TRUE)
abline(h = 0)
legend("topleft", lty = 1:3, legend = c("before", "after", "fitted"), cex = 2, horiz = TRUE)
grid(nx = 5, ny = 5, col = "lightgray", lty = "dotted")
mtext("Density", side = 2, line = 6.2, adj = 0.5, cex = 2)
mtext("Value (in HU)", side = 1, line = 2.7, adj = 0.5, cex = 2)
axis(1, mgp = c(3, 1.2, 0), cex.axis = 2, at = seq(from = -1400, to = -400, by = 100)) # Adjust second value to move x-axis labels
# Add custom y-axis with specific mgp settings for y-axis
axis(2, mgp = c(2, 1, 0), cex.axis = 2, las = 1) # Adjust second value to move y-axis labels
segments(x0 = mu1, y0 = 0, x1 = mu1, y1 = dnorm(mu1, mean = mu1, sd = sd1), lty = 1)
segments(x0 = mu2, y0 = 0, x1 = mu2, y1 = dnorm(mu2, mean = mu2, sd = sd2), lty = 2)
segments(x0 = mu3, y0 = 0, x1 = mu3, y1 = dnorm(mu3, mean = mu3, sd = sd3), lty = 3)
dev.off()


# Figure 4.5
patient.index <- 13
muqxnp1 <- res.qlm$muqx[patient.index]
sigmaqxnp1 <- res.qlm$sigmaqx[patient.index]
result <- optim(par = c(-828, 122), fn = objective_function, n = n, beta0hat = beta0hat, beta1hat = beta1hat, sigma2hat = sigma2hat, muqxnp1 = muqxnp1, muqxbar = muqxbar, w = w, betahat = betahat, beta2hat = beta2hat, sigmaqxnp1 = sigmaqxnp1, sigmaqxbar = sigmaqxbar, method = "L-BFGS-B",
                lower = c(-900, 100), upper = c(-750, 157), control = list(fnscale = -1))
max_value <- result$value
library("ggplot2")
library("dplyr")
# Set up the s and t grid for the plot
s_vals <- seq(-860, -800, length.out = 100)
t_vals <- seq(100, 140, length.out = 100)
grid <- expand.grid(s = s_vals, t = t_vals)
# Calculate fQhat[s,t] for the grid points
grid$z <- apply(grid, 1, function(row) {
  fQhat(row[1], row[2], n, beta0hat, beta1hat, sigma2hat, muqxnp1, muqxbar, w, betahat, beta2hat, sigmaqxnp1, sigmaqxbar)
})
# Define breaks to include 0.000 at the bottom of the legend
legend_breaks <- seq(from = 0.000, by = 0.0005, to = max_value + 0.00001)
# Contour plot
# Define font sizes and styles
tickSize1 <- 18
tickSize2 <- 13
legendSize1 <- 15
legendKeySize <- 0.5
legendSize2 <- 15
fontSize1 <- 20
fontSize2 <- 15
contour_plot <- ggplot(grid, aes(x = s, y = t, fill = z)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma", breaks = legend_breaks, limits = c(0, max_value + 0.00001)) +  
  guides(fill = guide_colorbar(
    title = "Density",
    title.theme = element_text(size = legendSize1),
    label.theme = element_text(size = legendSize1),
    barwidth = unit(0.8, "lines"),
    barheight = unit(8, "lines")
  )) +
  labs(
    x = TeX(r'($\mu$)'), 
    y = TeX(r'($\sigma$)')
  ) +
  scale_x_continuous(breaks = seq(min(grid$s), max(grid$s), by = 10), expand = c(0, 0)) +  
  scale_y_continuous(breaks = seq(min(grid$t), max(grid$t), by = 10), expand = c(0, 0)) +  
  theme_minimal(base_size = fontSize1) +
  theme(
    axis.title.x = element_text(angle = 0, vjust = 2, hjust = 0.5),
    axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = -0.5), 
    axis.text.x = element_text(size = tickSize1, margin = margin(0)),  
    axis.text.y = element_text(size = tickSize1, margin = margin(0)),  
    axis.ticks.length = unit(0.1, "cm"),  
    axis.ticks = element_line(color = "black"),
    legend.text = element_text(size = legendSize1),
    legend.key.size = unit(legendKeySize, "cm"),
    plot.margin = margin(4, -8, -6, 2),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, -10),  # Increase space between plot and legend
    legend.box.spacing = unit(0.3, "cm"),       # Further increase horizontal space
    legend.spacing.x = unit(0, "cm"),
    plot.title = element_text(size = 25, face = "bold"),
    panel.grid = element_blank()
  ) +
  geom_point(aes(x = res.qlm$muqx[patient.index], y = res.qlm$sigmaqx[patient.index]), shape = 17, color = "chartreuse2", size = 4) + # filled triangle pointing upwards
  geom_point(aes(x = res.qlm$muqy[patient.index], y = res.qlm$sigmaqy[patient.index]), shape = 15, color = "blue", size = 4) + #  filled square
  geom_point(aes(x = get("muqyfitted", envir = environment(fitted(res.qlm)[[patient.index]])), 
                 y = get("sigmaqyfitted", envir = environment(fitted(res.qlm)[[patient.index]]))), 
             shape = 3, color = "black", size = 5) # plus sign (+)

# Save contour plot as PDF
cairo_pdf("fig4.5.pdf", width = 8.09, height = 6.33, pointsize = 10)
print(contour_plot)
dev.off()


# Figure 4.6
cairo_pdf("fig4.6.pdf", width = 13.33, height = 10, pointsize = 10)  # width and height in inches, matching approximately to 4000x3000 pixels at 300 dpi
par(mar = c(2, 6, 0.4, 0.4),  # Smaller margins for base plot settings
    cex.lab = 3,           # Larger axis labels
    cex.axis = 3)        # Larger tick labels
plot((table4.1[,3] - table4.1[,1])[order(abs(table4.1[,3] - table4.1[,1]))],
     type = "h", 
     xlab = "", 
     ylab = expression(mu[q[y]] - mu[q[x]]), 
     xaxt = "n",
     ylim = c(-150, 150))
text(x = 1:44,
     y = (table4.1[,3] - table4.1[,1])[order(abs(table4.1[,3] - table4.1[,1]))] + 
       3 * sign((table4.1[,3] - table4.1[,1])[order(abs(table4.1[,3] - table4.1[,1]))]), 
     labels = order(abs(table4.1[,3] - table4.1[,1])), 
     cex = 1.6)  # Increase size of labels above plot lines
abline(h = 0)
# Close PDF device
dev.off()
