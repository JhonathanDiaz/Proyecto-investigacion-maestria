library(readxl)

# ##########################################################
# MUESTREO DE GIBBS : #######################
# ##########################################################
# Funciones auxiliares para el paso M-H:

library(readxl)
psi_hj1 <- function(hj1, yj1) {
  b_2 <- exp(-0.5 * (hj1 + exp(-hj1) * yj1^2))
  return(b_2)
}
psi_hjt = function ( hjt , yjt ) {
  b_1= exp ( -1 / 2 * ( hjt + exp ( - hjt ) * yjt ^2) )
  b_1
}

# ##########################################################
# Datos iniciales
n =21000 # n´umero de iteraciones
burn =2000 # tama~no del periodo de calentamiento
n.cadenas =3 # n´umero de cadenas a simular
salto =15 # tomar valor cada " salto " iterciones
#region =1 # simulaci´on para la regi´on n´umero
seed =987 # semilla para recuperar simulaci´on
if ( seed >0) { set.seed ( seed ) }
a1 =0 ; b1 =1 # hiperpar´ametros phi
c1 =3 ; d1 =3 # hiperpar´ametros sigma2
e1 =0 ; f1 =10 # hiperpar´ametros mu
# Generaci´on de valores iniciales
in_mu = rnorm ( n.cadenas ,0 ,1)
in_phi = rnorm ( n.cadenas ,0 ,0.35)
in_sigma2 =1 / rgamma ( n.cadenas ,3 ,2)

# ##########################################################

mide.tiempo = proc.time ()



# Cargar datos
file_path <- "C:/Users/COJEDIAZ/Downloads/Currencies.xlsx"
retornos <- read_excel(file_path, sheet = "Retornos 2012")
precios <- read_excel(file_path, sheet = "Precios 2012")

# Obtener lista de monedas (nombres de columnas)
monedas <- names(retornos)[-1]  # O [-1] si la primera columna es fecha

# Lista para guardar resultados
resultados <- list()


procesar_moneda <- function(nombre_moneda) {
  cat("Procesando:", nombre_moneda, "\n")
  
  y1 <- retornos[[nombre_moneda]][1:2899]
  z <- precios[[nombre_moneda]][1:2899]
  y_real <- diff(log(z))
  N <- length(y1)
  
  espmu = e1 ; espphi = a1 ; espsigma2 = d1 / ( c1 -1)
  # ##########################################################
  # Almacenamiento de las cadenas
  SIMmu = matrix (0 ,n , n.cadenas )
  SIMphi = matrix (0 ,n , n.cadenas )
  SIMsigma2 = matrix (0 ,n , n.cadenas )
  Hjtn = matrix (0 ,n , N * n.cadenas ) # guarda simul . de h.
  Contador = matrix (0 ,N , n.cadenas ) # guarda valores aceptados en MH
  
  
  
  
  
  ###############################################################
  ###############################################################
  for (m in 1:n.cadenas) {
    
    # ##############################################################
    # Simulación Gibbs n=1 latentes iniciales #
    # h_j para t =1 ,2 ,... ,155; N =155 #
    # ##############################################################
    # Paso M-H dentro de GIBBS #
    # Primero para t=1 #
    anterior <- in_mu[m]  # valor inicial para h_j(1)
    u <- runif(1)
    y <- rnorm(1, in_mu[m], sqrt(in_sigma2[m]))  # propuesta
    num <- psi_hj1(y, y1[1])
    den <- psi_hj1(anterior, y1[1])
    rho <- min(1, num / den)  # probabilidad de aceptación
    Hjtn[1, (N * (m - 1) + 1)] <- anterior + (y - anterior) * (u <= rho)
    Contador[1, m] <- (u <= rho)
    
    # Simulación para t = 2, ..., N
    for (t in 2:N) {
      hj_ant <- Hjtn[1, (N * (m - 1) + t - 1)]
      anterior <- rnorm(1, hj_ant, sqrt(in_sigma2[m]))
      u <- runif(1)
      y <- rnorm(1, in_mu[m] + in_phi[m] * (hj_ant - in_mu[m]), sqrt(in_sigma2[m]))
      num <- psi_hjt(y, y1[t])
      den <- psi_hjt(anterior, y1[t])
      rho <- min(1, num / den)
      Hjtn[1, (N * (m - 1) + t)] <- anterior + (y - anterior) * (u <= rho)
      Contador[t, m] <- (u <= rho)
    }
    
    # ##############################################################
    # mu_1 -> n = 1
    h_1 <- Hjtn[1, (N * (m - 1) + 1):(N * (m - 1) + N)]
    
    A <- 1 / f1 + (1 + (1 - in_phi[m])^2 * (N - 1)) / in_sigma2[m]
    diferencias <- numeric(N)
    for (t in 2:N) {
      diferencias[t] <- h_1[t] - in_phi[m] * h_1[t - 1]
    }
    B <- e1 / f1 + (h_1[1] + (1 - in_phi[m]) * sum(diferencias)) / in_sigma2[m]
    C <- B / A
    D <- 1 / A
    SIMmu[1, m] <- rnorm(1, C, sqrt(D))
    
    # phi -> n = 1
    mu_1 <- SIMmu[1, m]
    diferencias1 <- numeric(N)
    for (t in 2:N) {
      diferencias1[t] <- (h_1[t - 1] - mu_1)^2
    }
    A <- 1 / b1 + sum(diferencias1) / in_sigma2[m]
    
    diferencias2 <- numeric(N)
    for (t in 2:N) {
      diferencias2[t] <- (h_1[t] - mu_1) * (h_1[t - 1] - mu_1)
    }
    B <- a1 / b1 + sum(diferencias2) / in_sigma2[m]
    C <- B / A
    D <- 1 / A
    SIMphi[1, m] <- rnorm(1, C, sqrt(D))
    
    # sigma2 -> n = 1
    phi_1 <- SIMphi[1, m]
    A <- c1 + N / 2
    diferencias <- numeric(N)
    for (t in 2:N) {
      diferencias[t] <- (h_1[t] - mu_1 - phi_1 * (h_1[t - 1] - mu_1))^2
    }
    B <- 0.5 * ((h_1[1] - mu_1)^2 + sum(diferencias)) + d1
    SIMsigma2[1, m] <- 1 / rgamma(1, A, B)
    
    # #############################################################
    # Simulación completa de Gibbs para iteraciones k = 2, ..., n
    for (k in 2:n) {
      
      mu_1 <- SIMmu[k - 1, m]
      phi_1 <- SIMphi[k - 1, m]
      sigma2_1 <- SIMsigma2[k - 1, m]
      
      # M-H para h(1)
      anterior <- Hjtn[k - 1, (N * (m - 1) + 1)]
      u <- runif(1)
      y <- rnorm(1, mu_1, sqrt(sigma2_1))
      num <- psi_hj1(y, y1[1])
      den <- psi_hj1(anterior, y1[1])
      rho <- min(1, num / den)
      Hjtn[k, (N * (m - 1) + 1)] <- anterior + (y - anterior) * (u <= rho)
      Contador[1, m] <- Contador[1, m] + (u <= rho)
      
      # M-H para h_j con t = 2, ..., N
      for (t in 2:N) {
        hj_ant <- Hjtn[k, (N * (m - 1) + t - 1)]
        anterior <- Hjtn[k - 1, (N * (m - 1) + t)]
        u <- runif(1)
        y <- rnorm(1, mu_1 + phi_1 * (hj_ant - mu_1), sqrt(sigma2_1))
        num <- psi_hjt(y, y1[t])
        den <- psi_hjt(anterior, y1[t])
        rho <- min(1, num / den)
        Hjtn[k, (N * (m - 1) + t)] <- anterior + (y - anterior) * (u <= rho)
        Contador[t, m] <- Contador[t, m] + (u <= rho)
      }
      
      # mu -> n = k
      h_1 <- Hjtn[k, (N * (m - 1) + 1):(N * (m - 1) + N)]
      A <- 1 / f1 + (1 + (1 - phi_1)^2 * (N - 1)) / sigma2_1
      diferencias <- numeric(N)
      for (t in 2:N) {
        diferencias[t] <- h_1[t] - phi_1 * h_1[t - 1]
      }
      B <- e1 / f1 + (h_1[1] + (1 - phi_1) * sum(diferencias)) / sigma2_1
      C <- B / A
      D <- 1 / A
      SIMmu[k, m] <- rnorm(1, C, sqrt(D))
      
      # phi -> n = k
      mu_1 <- SIMmu[k, m]
      diferencias1 <- numeric(N)
      for (t in 2:N) {
        diferencias1[t] <- (h_1[t - 1] - mu_1)^2
      }
      A <- 1 / b1 + sum(diferencias1) / sigma2_1
      
      diferencias2 <- numeric(N)
      for (t in 2:N) {
        diferencias2[t] <- (h_1[t] - mu_1) * (h_1[t - 1] - mu_1)
      }
      B <- a1 / b1 + sum(diferencias2) / sigma2_1
      C <- B / A
      D <- 1 / A
      SIMphi[k, m] <- rnorm(1, C, sqrt(D))
      
      # sigma2 -> n = k
      phi_1 <- SIMphi[k, m]
      A <- c1 + N / 2
      diferencias <- numeric(N)
      for (t in 2:N) {
        diferencias[t] <- (h_1[t] - mu_1 - phi_1 * (h_1[t - 1] - mu_1))^2
      }
      B <- 0.5 * ((h_1[1] - mu_1)^2 + sum(diferencias)) + d1
      SIMsigma2[k, m] <- 1 / rgamma(1, A, B)
    }
  }
  
  ############################################
  # CONVERGENCIA
  ###########################################
  
  
  library(lattice)  # Paquetes que sirven
  library(coda)     # Para hacer las pruebas
  
  # -------------------------------------
  burning <- burn + 1
  
  las_h1 <- matrix(0, N, n.cadenas)
  las_h2 <- matrix(0, N, n.cadenas)
  SIMmu1 <- matrix(0, n - burn, n.cadenas)
  SIMphi1 <- matrix(0, n - burn, n.cadenas)
  SIMsigma21 <- matrix(0, n - burn, n.cadenas)
  
  nuevo_tno <- (n - burn) / salto
  SIMmu2 <- matrix(0, nuevo_tno, n.cadenas)
  SIMphi2 <- matrix(0, nuevo_tno, n.cadenas)
  SIMsigma22 <- matrix(0, nuevo_tno, n.cadenas)
  
  mmu1 <- numeric(n.cadenas)
  mphi1 <- numeric(n.cadenas)
  msigma21 <- numeric(n.cadenas)
  vmu1 <- numeric(n.cadenas)
  vphi1 <- numeric(n.cadenas)
  vsigma21 <- numeric(n.cadenas)
  
  mmu2 <- numeric(n.cadenas)
  mphi2 <- numeric(n.cadenas)
  msigma22 <- numeric(n.cadenas)
  vmu2 <- numeric(n.cadenas)
  vphi2 <- numeric(n.cadenas)
  vsigma22 <- numeric(n.cadenas)
  
  for (m in 1:n.cadenas) {
    # Estimaciones a posteriori de las h’s
    for (t in 1:N) {
      las_h1[t, m] <- mean(Hjtn[burning:n, (N * (m - 1) + t)])
    }
    
    help.h <- matrix(0, nuevo_tno, N)
    for (t in 1:N) {
      for (i in 1:nuevo_tno) {
        help.h[i, t] <- Hjtn[(burn + salto * i), (N * (m - 1) + t)]
      }
    }
    
    for (t in 1:N) {
      las_h2[t, m] <- mean(help.h[, t])
    }
    
    # Estimaciones para muestra sin minar
    SIMmu1[, m] <- SIMmu[burning:n, m]
    SIMphi1[, m] <- SIMphi[burning:n, m]
    SIMsigma21[, m] <- SIMsigma2[burning:n, m]
    
    mmu1[m] <- mean(SIMmu1[, m])
    mphi1[m] <- mean(SIMphi1[, m])
    msigma21[m] <- mean(SIMsigma21[, m])
    
    vmu1[m] <- var(SIMmu1[, m])
    vphi1[m] <- var(SIMphi1[, m])
    vsigma21[m] <- var(SIMsigma21[, m])
    
    # Estimaciones para muestra minada
    for (i in 1:nuevo_tno) {
      SIMmu2[i, m] <- SIMmu1[salto * i, m]
      SIMphi2[i, m] <- SIMphi1[salto * i, m]
      SIMsigma22[i, m] <- SIMsigma21[salto * i, m]
      
      
      
    }
    
    mmu2[m] <- mean(SIMmu2[, m])
    mphi2[m] <- mean(SIMphi2[, m])
    msigma22[m] <- mean(SIMsigma22[, m])
    
    vmu2[m] <- var(SIMmu2[, m])
    vphi2[m] <- var(SIMphi2[, m])
    vsigma22[m] <- var(SIMsigma22[, m])}
  
  # -------------------------------------
  # Ruta de guardado
  ruta_guardado <- "C:/Users/COJEDIAZ/Downloads/RESULTADOS/Graficos"
  archivo_png <- file.path(ruta_guardado, paste0("convergencia_mcmc_", moneda, ".png"))
  png(archivo_png, width = 1000, height = 1000)
  m=1
  # Gráficos para la cadena m
  #par(mfrow = c(3, 3))
  par(mfrow = c(3, 3), oma = c(0, 0, 4, 0)) 
  plot(SIMmu2[, m], type = "l",
       ylim = c(mmu2[m] - 4 * 1.96 * sqrt(vmu2[m]),
                mmu2[m] + 4 * 1.96 * sqrt(vmu2[m])))
  abline(h = c(mmu2[m] - 1.96 * sqrt(vmu2[m]),
               mmu2[m] + 1.96 * sqrt(vmu2[m])), col = "blue")
  
  plot(SIMphi2[, m], type = "l",
       ylim = c(mphi2[m] - 4 * 1.96 * sqrt(vphi2[m]),
                mphi2[m] + 4 * 1.96 * sqrt(vphi2[m])))
  abline(h = c(mphi2[m] - 1.96 * sqrt(vphi2[m]),
               mphi2[m] + 1.96 * sqrt(vphi2[m])), col = "blue")
  
  plot(SIMsigma22[, m], type = "l",
       ylim = c(msigma22[m] - 4 * 1.96 * sqrt(vsigma22[m]),
                msigma22[m] + 4 * 1.96 * sqrt(vsigma22[m])))
  abline(h = c(msigma22[m] - 1.96 * sqrt(vsigma22[m]),
               msigma22[m] + 1.96 * sqrt(vsigma22[m])), col = "blue")
  
  hist(SIMmu2[, m], main = paste("Hist mu, cadena", m))
  hist(SIMphi2[, m], main = paste("Hist phi, cadena", m))
  hist(SIMsigma22[, m], main = paste("Hist sigma^2, cadena", m))
  
  acf(SIMmu2[, m], lag.max = 100)
  acf(SIMphi2[, m], lag.max = 100)
  acf(SIMsigma22[, m], lag.max = 100)
  mtext(paste("CONVERGENCIA MONEDA", moneda), outer = TRUE, cex = 1.5, font = 2)
  dev.off()  # Cierra el archivo
  
  # -------------------------------------
  # Reporte por consola
  ruta_base <- "C:/Users/COJEDIAZ/Downloads/RESULTADOS/Graficos"
  ruta_txt <- file.path(ruta_base, paste0("resultados_diagnosticos_", moneda, ".txt"))
  sink(ruta_txt)
  
  simtime <- proc.time() - mide.tiempo
  
  rbind ( paste ( " Resultados para muestra usando las " ,n.cadenas , " cadenas. " ) ,
          paste ( " Par´ametros estimados cada " , salto , " datos : " , mean ( mmu2 ) , mean (mphi2 ) , mean ( msigma22 ) ) ,
          paste ( " Par´ametros estimados muestra completa : " , mean ( mmu1 ) , mean (mphi1 ) , mean ( msigma21 ) ) ,
          paste ( " En todos los casos el calentamiento fue de " , burn ) ,
          paste ( " Tama~no de la muestra cada cadena (n - burn ) : " ,n - burn ) ,
          paste ( " Tama~no de la muestra cada cadena usando cada " , salto , " datos(( n - burn ) / salto ) : " , nuevo_tno ) ,
          paste ( " Resultados para cadena " ,m ) ,
          paste ( " Valores iniciales : " ," mu0 = " , in_mu [ m ] , " phi0 = " , in_phi [ m ] , "sigma20 = " , in_sigma2 [ m ]) ,
          paste ( " Par´ametros estimados cada " , salto , " datos : " , mmu2 [ m ] , mphi2 [ m ] ,msigma22 [ m ]) ,
          paste ( " Varianzas muestrales cada " , salto , " datos : " , vmu2 [ m ] , vphi2 [ m ] ,vsigma22 [ m ]) ,
          paste ( " Par´ametros estimados muestra completa : " , mmu1 [ m ] , mphi1 [ m ] ,msigma21 [ m ]) ,
          paste ( " Varianzas muestrales con muestra completa : " , vmu1 [ m ] , vphi1 [ m
          ] , vsigma21 [ m ]) ,
          paste ( " Hiperp ( mu , phi , sigma2 ) = " ,e1 , f1 , a1 , b1 , c1 , d1 ) ,
          paste ( " Esperanzas a priori ( mu , phi , sigma2 ) = " , espmu , espphi , espsigma2
          ) ,
          paste ( " Varianzas a priori ( mu , phi , sigma2 ) = " ,f1 , b1 ,( espsigma2 ^2 / ( c1
                                                                                              -2) ) ) ,
          paste ( " Tiempo simulaci´on : " , simtime [3]) )
  
  ######### Pruebas con coda de R
  ##1: Geweke
  Gibbs_mu1 = mcmc ( SIMmu1 [ , m ])
  Gibbs_phi1 = mcmc ( SIMphi1 [ , m ])
  Gibbs_sigma21 = mcmc ( SIMsigma21 [ , m ])
  geweke.diag( Gibbs_mu1 )
  geweke.diag( Gibbs_phi1 )
  geweke.diag( Gibbs_sigma21 )
  ##
  Gibbs_mu2 = mcmc ( SIMmu2 [ , m ])
  Gibbs_phi2 = mcmc ( SIMphi2 [ , m ])
  Gibbs_sigma22 = mcmc ( SIMsigma22 [ , m ])
  geweke.diag( Gibbs_mu2 )
  geweke.diag( Gibbs_phi2 )
  geweke.diag( Gibbs_sigma22 )
  
  ##2: Gelman and Rubin
  GibbsGR_mu11 = mcmc ( SIMmu1 [ ,1])
  GibbsGR_phi11 = mcmc ( SIMphi1 [ ,1])
  GibbsGR_sigma211 = mcmc ( SIMsigma21 [ ,1])
  GibbsGR_mu12 = mcmc ( SIMmu1 [ ,2])
  GibbsGR_phi12 = mcmc ( SIMphi1 [ ,2])
  GibbsGR_sigma212 = mcmc ( SIMsigma21 [ ,2])
  GibbsGR_mu13 = mcmc ( SIMmu1 [ ,3])
  GibbsGR_phi13 = mcmc ( SIMphi1 [ ,3])
  GibbsGR_sigma213 = mcmc ( SIMsigma21 [ ,3])
  GibbsGR_mu1 = mcmc.list ( list ( GibbsGR_mu11 , GibbsGR_mu12 , GibbsGR_mu13 ) )
  GibbsGR_phi1 = mcmc.list ( list ( GibbsGR_phi11 , GibbsGR_phi12 , GibbsGR_phi13 ) )
  GibbsGR_sigma21 = mcmc.list ( list ( GibbsGR_sigma211 , GibbsGR_sigma212 ,
                                       GibbsGR_sigma213 ) )
  gelman.diag( GibbsGR_mu1 )
  gelman.diag( GibbsGR_phi1 )
  gelman.diag( GibbsGR_sigma21 )
  ##
  GibbsGR_mu21 = mcmc ( SIMmu2 [ ,1])
  GibbsGR_phi21 = mcmc ( SIMphi2 [ ,1])
  GibbsGR_sigma221 = mcmc ( SIMsigma22 [ ,1])
  GibbsGR_mu22 = mcmc ( SIMmu2 [ ,2])
  GibbsGR_phi22 = mcmc ( SIMphi2 [ ,2])
  GibbsGR_sigma222 = mcmc ( SIMsigma22 [ ,2])
  GibbsGR_mu23 = mcmc ( SIMmu2 [ ,3])
  GibbsGR_phi23 = mcmc ( SIMphi2 [ ,3])
  GibbsGR_sigma223 = mcmc ( SIMsigma22 [ ,3])
  GibbsGR_mu2 = mcmc.list( list( GibbsGR_mu21 , GibbsGR_mu22 , GibbsGR_mu23 ) )
  GibbsGR_phi2 = mcmc.list( list( GibbsGR_phi21 , GibbsGR_phi22 , GibbsGR_phi23 ) )
  GibbsGR_sigma22 = mcmc.list ( list ( GibbsGR_sigma221 , GibbsGR_sigma222 ,
                                       GibbsGR_sigma223 ) )
  gelman.diag ( GibbsGR_mu2 )
  gelman.diag ( GibbsGR_phi2 )
  gelman.diag ( GibbsGR_sigma22 )
  # gelman . plot ( GibbsGR _mu1)
  
  ##3: Raftery and Lewis
  raftery.diag( Gibbs_mu1 )
  raftery.diag( Gibbs_phi1 )
  raftery.diag( Gibbs_sigma21 )
  ##
  raftery.diag( Gibbs_mu2 )
  raftery.diag( Gibbs_phi2 )
  raftery.diag( Gibbs_sigma22 )
  
  ##4: Heidelberger and Welch
  heidel.diag ( Gibbs_mu1 )
  heidel.diag( Gibbs_phi1 )
  heidel.diag( Gibbs_sigma21 )
  ##
  heidel.diag( Gibbs_mu2 )
  heidel.diag( Gibbs_phi2 )
  heidel.diag( Gibbs_sigma22 )
  sink()
  cat("✅ Resultados guardados en:\n", ruta_txt, "\n")
  # Ruta base donde guardar
  ruta_base <- "C:/Users/COJEDIAZ/Downloads/RESULTADOS/Graficos"
  
  # Construcción dinámica del nombre del archivo
  ruta_txt <- file.path(ruta_base, paste0("resultados_diagnosticos_", moneda, ".txt"))
  
  # Captura de salida de consola en el archivo
  sink(ruta_txt)
  
  cat("###############################\n")
  cat("### RESUMEN DE RESULTADOS\n")
  cat("### Moneda:", moneda, "\n")
  cat("### Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("###############################\n\n")
  
  simtime <- proc.time() - mide.tiempo
  print(
    rbind(
      paste("Resultados para muestra usando las", n.cadenas, "cadenas."),
      paste("Parámetros estimados cada", salto, "datos:", mean(mmu2), mean(mphi2), mean(msigma22)),
      paste("Parámetros estimados muestra completa:", mean(mmu1), mean(mphi1), mean(msigma21)),
      paste("Calentamiento (burn-in):", burn),
      paste("Tamaño muestra post-burn-in:", n - burn),
      paste("Muestra usada cada", salto, "datos:", nuevo_tno),
      paste("Resultados para cadena", m),
      paste("Valores iniciales: mu0 =", in_mu[m], "phi0 =", in_phi[m], "sigma20 =", in_sigma2[m]),
      paste("Estimaciones cada", salto, "datos:", mmu2[m], mphi2[m], msigma22[m]),
      paste("Varianzas muestrales:", vmu2[m], vphi2[m], vsigma22[m]),
      paste("Estimaciones muestra completa:", mmu1[m], mphi1[m], msigma21[m]),
      paste("Varianzas completas:", vmu1[m], vphi1[m], vsigma21[m]),
      paste("Hiperparámetros:", e1, f1, a1, b1, c1, d1),
      paste("Esperanzas a priori:", espmu, espphi, espsigma2),
      paste("Varianzas a priori:", f1, b1, (espsigma2^2 / (c1 - 2))),
      paste("Tiempo de simulación:", simtime[3])
    )
  )
  
  cat("\n###############################\n")
  cat("### DIAGNÓSTICOS CODA\n")
  cat("###############################\n")
  
  # -- Geweke
  cat("\n--- Geweke\n")
  print(geweke.diag(Gibbs_mu1))
  print(geweke.diag(Gibbs_phi1))
  print(geweke.diag(Gibbs_sigma21))
  print(geweke.diag(Gibbs_mu2))
  print(geweke.diag(Gibbs_phi2))
  print(geweke.diag(Gibbs_sigma22))
  
  # -- Gelman-Rubin
  cat("\n--- Gelman-Rubin\n")
  print(gelman.diag(GibbsGR_mu1))
  print(gelman.diag(GibbsGR_phi1))
  print(gelman.diag(GibbsGR_sigma21))
  print(gelman.diag(GibbsGR_mu2))
  print(gelman.diag(GibbsGR_phi2))
  print(gelman.diag(GibbsGR_sigma22))
  
  # -- Raftery y Lewis
  cat("\n--- Raftery y Lewis\n")
  print(raftery.diag(Gibbs_mu1))
  print(raftery.diag(Gibbs_phi1))
  print(raftery.diag(Gibbs_sigma21))
  print(raftery.diag(Gibbs_mu2))
  print(raftery.diag(Gibbs_phi2))
  print(raftery.diag(Gibbs_sigma22))
  
  # -- Heidelberger y Welch
  cat("\n--- Heidelberger y Welch\n")
  print(heidel.diag(Gibbs_mu1))
  print(heidel.diag(Gibbs_phi1))
  print(heidel.diag(Gibbs_sigma21))
  print(heidel.diag(Gibbs_mu2))
  print(heidel.diag(Gibbs_phi2))
  print(heidel.diag(Gibbs_sigma22))
  
  # Finalizar la redirección
  sink()
  
  cat("✅ Resultados guardados en:\n", ruta_txt, "\n")
  ##########################
  # Ruta base donde guardar
  ruta_base <- "C:/Users/COJEDIAZ/Downloads/RESULTADOS/Graficos"
  
  # Construcción dinámica del nombre del archivo
  ruta_txt <- file.path(ruta_base, paste0("resultados_diagnosticos_", moneda, ".txt"))
  
  # Captura de salida de consola en el archivo
  sink(ruta_txt)
  
  cat("###############################\n")
  cat("### RESUMEN DE RESULTADOS\n")
  cat("### Moneda:", moneda, "\n")
  cat("### Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("###############################\n\n")
  
  simtime <- proc.time() - mide.tiempo
  print(
    rbind(
      paste("Resultados para muestra usando las", n.cadenas, "cadenas."),
      paste("Parámetros estimados cada", salto, "datos:", mean(mmu2), mean(mphi2), mean(msigma22)),
      paste("Parámetros estimados muestra completa:", mean(mmu1), mean(mphi1), mean(msigma21)),
      paste("Calentamiento (burn-in):", burn),
      paste("Tamaño muestra post-burn-in:", n - burn),
      paste("Muestra usada cada", salto, "datos:", nuevo_tno),
      paste("Resultados para cadena", m),
      paste("Valores iniciales: mu0 =", in_mu[m], "phi0 =", in_phi[m], "sigma20 =", in_sigma2[m]),
      paste("Estimaciones cada", salto, "datos:", mmu2[m], mphi2[m], msigma22[m]),
      paste("Varianzas muestrales:", vmu2[m], vphi2[m], vsigma22[m]),
      paste("Estimaciones muestra completa:", mmu1[m], mphi1[m], msigma21[m]),
      paste("Varianzas completas:", vmu1[m], vphi1[m], vsigma21[m]),
      paste("Hiperparámetros:", e1, f1, a1, b1, c1, d1),
      paste("Esperanzas a priori:", espmu, espphi, espsigma2),
      paste("Varianzas a priori:", f1, b1, (espsigma2^2 / (c1 - 2))),
      paste("Tiempo de simulación:", simtime[3])
    )
  )
  
  cat("\n###############################\n")
  cat("### DIAGNÓSTICOS CODA\n")
  cat("###############################\n")
  
  # -- Geweke
  cat("\n--- Geweke\n")
  print(geweke.diag(Gibbs_mu1))
  print(geweke.diag(Gibbs_phi1))
  print(geweke.diag(Gibbs_sigma21))
  print(geweke.diag(Gibbs_mu2))
  print(geweke.diag(Gibbs_phi2))
  print(geweke.diag(Gibbs_sigma22))
  
  # -- Gelman-Rubin
  cat("\n--- Gelman-Rubin\n")
  print(gelman.diag(GibbsGR_mu1))
  print(gelman.diag(GibbsGR_phi1))
  print(gelman.diag(GibbsGR_sigma21))
  print(gelman.diag(GibbsGR_mu2))
  print(gelman.diag(GibbsGR_phi2))
  print(gelman.diag(GibbsGR_sigma22))
  
  # -- Raftery y Lewis
  cat("\n--- Raftery y Lewis\n")
  print(raftery.diag(Gibbs_mu1))
  print(raftery.diag(Gibbs_phi1))
  print(raftery.diag(Gibbs_sigma21))
  print(raftery.diag(Gibbs_mu2))
  print(raftery.diag(Gibbs_phi2))
  print(raftery.diag(Gibbs_sigma22))
  
  # -- Heidelberger y Welch
  cat("\n--- Heidelberger y Welch\n")
  print(heidel.diag(Gibbs_mu1))
  print(heidel.diag(Gibbs_phi1))
  print(heidel.diag(Gibbs_sigma21))
  print(heidel.diag(Gibbs_mu2))
  print(heidel.diag(Gibbs_phi2))
  print(heidel.diag(Gibbs_sigma22))
  
  # Finalizar la redirección
  sink()
  
  cat("✅ Resultados guardados en:\n", ruta_txt, "\n")
  
  
  
    
  ###############Pronostico 
  grafo <- function(seed, m, p, s2, z, y_real, region = "Moneda", linf = 2700, horizon = 4) {
    set.seed(seed)
    
    n_real <- length(z)
    n_total <- n_real + horizon  # total de puntos simulados
    
    h <- numeric(n_total)
    z_est <- numeric(n_total)
    y_est <- numeric(n_total)
    volatilidad <- numeric(n_total)
    
    # Epsilon truncado (opción 1)
    epsilon <- pmax(pmin(rnorm(n_total), 2), -2)
    
    # Inicialización
    h[1] <- rnorm(1, m, sqrt(s2))
    volatilidad[1] <- exp(h[1])
    z_est[1] <- z[1]
    
    # Simulación de datos observados
    for (i in 2:n_real) {
      h[i] <- rnorm(1, m + p * (h[i - 1] - m), sqrt(s2))
      volatilidad[i] <- exp(h[i])
      z_est[i] <- z[i - 1] * exp(0.5 * sqrt(volatilidad[i]) * epsilon[i])
      y_est[i] <- log(z_est[i] / z_est[i - 1])
    }
    
    # Simulación futura (proyección)
    for (i in (n_real + 1):n_total) {
      h[i] <- rnorm(1, m + p * (h[i - 1] - m), sqrt(s2))
      volatilidad[i] <- exp(h[i])
      z_est[i] <- z_est[i - 1] * exp(0.5 * sqrt(volatilidad[i]) * epsilon[i])
      y_est[i] <- log(z_est[i] / z_est[i - 1])
    }
    
    # Matrices para graficar
    graf <- cbind(c(z, rep(NA, horizon)), z_est)
    graf_y <- cbind(c(NA, y_real, rep(NA, n_total - length(y_real) - 1)), y_est)
    
    errores <- abs(graf[, 1] - graf[, 2])
    error_total <- sum(errores[1:n_real], na.rm = TRUE)
    vol_total <- sum(volatilidad, na.rm = TRUE)
    
    # Gráfico
    par(mfrow = c(1, 1))
    ylim_range <- range(graf, na.rm = TRUE) * c(0.95, 1.05)
    
    plot(1:n_total, graf[, 1], type = 'n', lwd = 2.5,
         ylim = ylim_range,
         ylab = "Precio diario",
         xlab = "Día",
         main = paste("Simulación para", region))
    
    abline(h = seq(ylim_range[1], ylim_range[2], length.out = 10), col = 'gray', lwd = 0.5)
    abline(v = seq(0, n_total, by = 250), col = 'gray', lwd = 0.5)
    abline(h = 0, lwd = 1.5)
    
    legend("topright", c("Observada", "Estimada"), col = c(1, 2), lty = 1:2, lwd = 2)
    
    lines(1:n_real, graf[1:n_real, 1], lwd = 2.5)
    lines((n_real + 1):n_total, graf[(n_real + 1):n_total, 1], type = 'o', col = "blue", lwd = 2.5)
    lines(1:n_total, graf[, 2], type = 'o', lwd = 2.5, col = 2, lty = 2)
    
    return(list(
      error_total = error_total,
      volatilidad_total = vol_total,
      z_observado = z,
      z_est = z_est,
      y_est = y_est,
      volatilidad = volatilidad,
      h = h,
      epsilon = epsilon,
      graf = graf,
      graf_y = graf_y
    ))
  }
  
  # PARÁMETROS FINAL ESTIMADOS
  m <- mean(SIMmu2[, 1])
  p <- mean(SIMphi2[, 1])
  s2 <- mean(SIMsigma22[, 1])
  # Crear archivo de salida PNG
  ruta_guardado <- "C:/Users/COJEDIAZ/Downloads/RESULTADOS/Graficos"
  archivo_png <- file.path(ruta_guardado, paste0("pronostico_", moneda, ".png"))
  
  png(filename = archivo_png, width = 1000, height = 600)
  
  # Simulación con grafo
  res <- grafo(seed = 123, m = m, p = p, s2 = s2,
               z = z, y_real = y_real, region = nombre_moneda, linf = 2700)

  dev.off()
  
  cat("✅ Gráfico de pronóstico guardado en:\n", archivo_png, "\n")


  
  return(res)
}
  
for (moneda in monedas) {
  res <- procesar_moneda(moneda)
  assign(paste0("res_", moneda), res, envir = .GlobalEnv)
}



# Define la ruta completa con el nombre del archivo
ruta_guardado <- "C:/Users/COJEDIAZ/Downloads/RESULTADOS/resultados_monedas.RData"

# Guarda todos los objetos que comienzan con "res_"
save(list = ls(pattern = "^res_"), file = ruta_guardado)





cat("✅ Objetos res_* guardados en:\n", ruta_guardado, "\n")

    
    
library(LaplacesDemon)

monedas <- c("AUD", "EUR", "GBP", "NZD", "CAD", "CHF", "CNH", "JPY", "MXN", "NOK", "SEK", "TRY", "ZAR", "XAU")
coins <- length(monedas)
days <- length(res_AUD$h)    
    
    
# Crear matriz de desviaciones estándar diarias
h_t2 <- matrix(nrow = days, ncol = coins)
for (i in 1:coins) {
  res <- get(paste0("res_", monedas[i]))
  h_t2[, i] <- res$h
}

    
# Crear distribuciones normales para cada moneda y día
x <- seq(-0.0005, 0.0005, by = 0.0000005)
norm <- matrix(nrow = length(x), ncol = coins * days)

for (i in 1:coins) {
  for (j in 1:days) {
    sigma <- sqrt(exp(h_t2[j, i]))
    norm[, (i - 1) * days + j] <- dnorm(x, mean = 0, sd = sigma)
  }
}

#### Calcular KL de todas las monedas respecto a EUR
idx_base <- which(monedas == "EUR")
KLlist <- numeric(coins - 1)
names(KLlist) <- monedas[-idx_base]  # asignar nombres
for (i in 1:coins) {
  if (i == idx_base) next
  sumaKL <- 0
  for (j in 1:days) {
    col_base <- (idx_base - 1) * days + j
    col_comp <- (i - 1) * days + j
    kl <- KLD(norm[, col_base], norm[, col_comp])
    sumaKL <- sumaKL + kl$sum.KLD.px.py
  }
  KLlist[i - (i > idx_base)] <- sumaKL
}



############# TODAS LAS MONEDAS##########################################################################


library(LaplacesDemon)

monedas <- c("AUD", "EUR", "GBP", "NZD", "CAD", "CHF", "CNH", "JPY", "MXN", "NOK", "SEK", "TRY", "ZAR", "XAU")
coins <- length(monedas)
days <- length(res_AUD$h)    


# Crear matriz de desviaciones estándar diarias
h_t2 <- matrix(nrow = days, ncol = coins)
for (i in 1:coins) {
  res <- get(paste0("res_", monedas[i]))
  h_t2[, i] <- res$h
}


# Crear distribuciones normales para cada moneda y día
x <- seq(-0.0005, 0.0005, by = 0.0000005)




norm <- matrix(nrow = length(x), ncol = coins * (days - 1))


 for (i in 1:coins) {
     for (j in 1:(days - 1)) {
           
             h_val <- h_t2[j, i]
            
               # Cálculo seguro de sigma
               sigma <- sqrt(exp(h_val))
               sigma <- max(sigma, 1e-6)  #esto me evita sd = 0 o valores muy pequeños
               
                 # Densidad normal + normalización
                 pdf <- dnorm(x, mean = 0, sd = sigma)
                 pdf <- pdf / sum(pdf)  # esto me asegura que la suma sea 1 (recomendado para KL)
                 
                   # Asignación
                   norm[, (i - 1) * (days - 1) + j] <- pdf
                   }
 }

for (i in 1:coins) {
  for (j in 1:(days - 1)) {
    
    h_val <- h_t2[j, i]
    
    # Cálculo seguro de sigma
    sigma <- sqrt(exp(h_val))
    sigma <- max(sigma, 1e-6)  # evita sd = 0 o valores muy pequeños
    
    # Densidad normal SIN normalización manual
    pdf <- dnorm(x, mean = 0, sd = sigma)
    
    # Asignación directa
    norm[, (i - 1) * (days - 1) + j] <- pdf
  }
}



i=1
j=1

library(LaplacesDemon)
KL_matrix <- matrix(0, nrow = coins, ncol = coins)
for (i in 1:coins) {
     for (j in 1:coins) {
           if (i == j) next
           auxsum <- 0
             for (k in 1:(days-1)) {
                   col_i <- (i - 1) * (days - 1) + k
                   col_j <- (j - 1) * (days - 1) + k
                   KLdist <- KLD(norm[, col_i], norm[, col_j])
                   auxsum <- auxsum + KLdist$sum.KLD.px.py
               }
          KL_matrix[i, j] <- auxsum
     }
  }


rownames(KL_matrix) <- monedas
colnames(KL_matrix) <- monedas

library(writexl)
write_xlsx(KL_df, path = "C:/Users/COJEDIAZ/Downloads/RESULTADOS/KL_matrix2.xlsx")

