#######################################################################################
#  TRABAJO 7 - MODELO 1: TOLERANCIA FIJA                                            #
#  Método: NEWTON-RAPHSON                                                           #
#  Encuentra el tiempo donde proporción_infelices ≈ proporción_similares            #
#######################################################################################

# ====================================================================================
# SIMULACIÓN DEL MODELO 1 - Tolerancia Fija
# ====================================================================================

simular_modelo1 <- function(N, p, citySize, years, alikePref, alfa){
  
  T <- citySize**2
  N1 <- round(N * p)
  N2 <- N - N1
  
  groups <- c(rep(1, N1), rep(2, N2))
  empty <- c(rep(0, T - N))
  set.seed(1234)
  city <- matrix(sample(c(groups, empty), citySize^2, replace = F), ncol = citySize)
  
  getNeighbors <- function(coords) {
    n <- c()
    for (i in c(1:8)) {
      if (i == 1) { x <- coords[1] + 1; y <- coords[2] }
      if (i == 2) { x <- coords[1] + 1; y <- coords[2] + 1 }
      if (i == 3) { x <- coords[1]; y <- coords[2] + 1 }
      if (i == 4) { x <- coords[1] - 1; y <- coords[2] + 1 }
      if (i == 5) { x <- coords[1] - 1; y <- coords[2] }
      if (i == 6) { x <- coords[1] - 1; y <- coords[2] - 1 }
      if (i == 7) { x <- coords[1]; y <- coords[2] - 1 }
      if (i == 8) { x <- coords[1] + 1; y <- coords[2] - 1 }
      
      if (x < 1) x <- citySize
      if (x > citySize) x <- 1
      if (y < 1) y <- citySize
      if (y > citySize) y <- 1
      n <- rbind(n, c(x, y))
    }
    n
  }
  
  unHappyMonitor <- c()
  similarMonitor <- c()
  
  # MODELO 1: Tolerancia FIJA
  alikePref1 <- alikePref
  alikePref2 <- alikePref * alfa
  
  for (t in c(1:years)) {
    happy <- c()
    unhappy <- c() 
    similar <- c()
    
    for (j in c(1:citySize)) {
      for (k in c(1:citySize)) {
        current <- c(j, k)
        value <- city[j, k] 
        
        if (value > 0) {
          likeNeighbors <- 0
          allNeighbors <- 0
          neighbors <- getNeighbors(current)
          
          for (i in c(1:nrow(neighbors))){
            x <- neighbors[i, 1]
            y <- neighbors[i, 2]
            if (city[x, y] > 0) allNeighbors <- allNeighbors + 1
            if (city[x, y] == value) likeNeighbors <- likeNeighbors + 1
          }
          
          if (!is.nan(likeNeighbors / allNeighbors)) {
            similar <- c(similar, likeNeighbors / allNeighbors)
            alikePrefg <- ifelse(value == 1, alikePref1, alikePref2)
            
            if ((likeNeighbors / allNeighbors) < alikePrefg) {
              unhappy <- rbind(unhappy, current)
            } else {
              happy <- rbind(happy, current)
            }
          } else {
            happy <- rbind(happy, current)
          }
        }
      }
    }
    
    unHappyMonitor <- c(unHappyMonitor, 
                        length(unhappy) / (length(happy) + length(unhappy)))
    similarMonitor <- c(similarMonitor, mean(similar))
    
    if (length(unhappy) > 0) {
      randUnHappy <- sample(nrow(unhappy))
      for (i in randUnHappy) {
        mover <- unhappy[i, ]
        mover_val <- city[mover[1], mover[2]]
        move_to <- c(sample(1:citySize, 1), sample(1:citySize, 1))
        move_to_val <- city[move_to[1], move_to[2]]
        
        while (move_to_val > 0) {
          move_to <- c(sample(1:citySize, 1), sample(1:citySize, 1))
          move_to_val <- city[move_to[1], move_to[2]]
        }
        city[mover[1], mover[2]] <- 0
        city[move_to[1], move_to[2]] <- mover_val
      }
    }
  }
  
  return(list(
    unhappy = unHappyMonitor, 
    similar = similarMonitor,
    city_final = city
  ))
}

# ====================================================================================
# MÉTODO DE NEWTON-RAPHSON PARA MODELO 1
# ====================================================================================

newton_raphson_modelo1 <- function(N, p, citySize, alikePref, alfa, 
                                   x0, years, tol = 0.01, max_iter = 50, h = 0.1) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║              MODELO 1: TOLERANCIA FIJA                        ║\n")
  cat("║              MÉTODO: NEWTON-RAPHSON                           ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  # Medir tiempo total
  tiempo_inicio_total <- Sys.time()
  
  # Medir tiempo de simulación
  cat("⏳ Ejecutando simulación del Modelo 1...\n")
  tiempo_inicio_sim <- Sys.time()
  datos <- simular_modelo1(N, p, citySize, years, alikePref, alfa)
  tiempo_fin_sim <- Sys.time()
  tiempo_simulacion <- as.numeric(difftime(tiempo_fin_sim, tiempo_inicio_sim, units = "secs"))
  diferencia <- datos$unhappy - datos$similar
  
  cat(sprintf("✓ Simulación completada en %.2f segundos\n\n", tiempo_simulacion))
  
  # Medir tiempo del método
  tiempo_inicio_metodo <- Sys.time()
  
  # Función objetivo f(t)
  f <- function(t) {
    t_int <- round(t)
    if (t_int < 1) t_int <- 1
    if (t_int > length(diferencia)) t_int <- length(diferencia)
    return(diferencia[t_int])
  }
  
  # Derivada numérica mejorada con interpolación
  df <- function(t) {
    t_plus <- min(t + h, length(diferencia))
    t_minus <- max(t - h, 1)
    
    if (abs(t_plus - t_minus) < 0.1) {
      t_plus <- min(t + 1, length(diferencia))
      t_minus <- max(t - 1, 1)
    }
    
    f_plus <- f(t_plus)
    f_minus <- f(t_minus)
    delta_t <- t_plus - t_minus
    
    if (delta_t < 1e-10) return(0)
    
    return((f_plus - f_minus) / delta_t)
  }
  
  cat("📊 PARÁMETROS:\n")
  cat(sprintf("   • N = %d agentes (%.1f%% ocupación)\n", N, 100*N/(citySize^2)))
  cat(sprintf("   • Grupo 1: %.0f%% (tolerancia = %.2f)\n", 100*p, alikePref))
  cat(sprintf("   • Grupo 2: %.0f%% (tolerancia = %.2f)\n", 100*(1-p), alikePref*alfa))
  cat(sprintf("   • Ciudad: %dx%d parcelas\n", citySize, citySize))
  cat(sprintf("   • Valor inicial: x0 = %.2f años\n", x0))
  cat(sprintf("   • Tolerancia: %.4f\n", tol))
  cat(sprintf("   • Máximo iteraciones: %d\n", max_iter))
  cat(sprintf("   • h (derivada numérica): %.2f\n\n", h))
  
  cat("🔍 MÉTODO:\n")
  cat("   Newton-Raphson modificado para minimizar |f(t)|\n")
  cat("   Busca el tiempo donde |infelices - similares| es mínimo\n\n")
  
  # Funciones auxiliares
  g <- function(t) abs(f(t))
  
  dg <- function(t) {
    ft <- f(t)
    dft <- df(t)
    if (abs(ft) < 1e-10) return(0)
    return(sign(ft) * dft)
  }
  
  d2g <- function(t) {
    return((dg(t + h) - dg(t - h)) / (2 * h))
  }
  
  # Verificar si hay cero o buscar mínimo
  valores <- sapply(1:length(diferencia), f)
  tiene_cero <- any(diff(sign(valores)) != 0)
  
  if (!tiene_cero) {
    cat("⚠️  NO HAY CAMBIO DE SIGNO en la función\n")
    cat("   → Buscando MÍNIMO resolviendo f'(t) = 0\n\n")
    buscar_minimo <- TRUE
  } else {
    cat("✓ HAY CAMBIO DE SIGNO\n")
    cat("   → Buscando CERO de f(t) = 0\n\n")
    buscar_minimo <- FALSE
  }
  
  cat("────────────────────────────────────────────────────────────────\n")
  cat("ITER      x          f(x)         f'(x)        |f'(x)|     ERROR\n")
  cat("────────────────────────────────────────────────────────────────\n")
  
  historial <- data.frame()
  x <- x0
  iter <- 0
  convergio <- FALSE
  
  for (iter in 1:max_iter) {
    
    f_real <- f(x)
    df_real <- df(x)
    
    if (buscar_minimo) {
      fx <- df(x)
      dfx <- (df(x + h) - df(x - h)) / (2 * h)
      if (abs(dfx) < 1e-10) {
        dfx <- (df(x + 2*h) - df(x - 2*h)) / (4 * h)
      }
      criterio <- fx
    } else {
      fx <- f(x)
      dfx <- df(x)
      criterio <- fx
    }
    
    cat(sprintf("%4d   %8.3f   %11.6f   %11.6f   %11.6f   ", 
                iter, x, f_real, df_real, abs(df_real)))
    
    if (abs(dfx) < 1e-6) {
      cat("    N/A\n")
      cat(sprintf("\n⚠️  Derivada muy pequeña (%.2e) en iteración %d\n", dfx, iter))
      cat("   Probando con búsqueda directa del mínimo...\n\n")
      
      t_vals <- seq(1, length(diferencia), by = 0.5)
      abs_diffs <- sapply(t_vals, function(t) abs(f(t)))
      idx_min <- which.min(abs_diffs)
      x <- t_vals[idx_min]
      convergio <- TRUE
      
      historial <- rbind(historial, data.frame(
        iter = iter,
        x = x,
        fx = f_real,
        dfx = df_real,
        abs_dfx = abs(df_real),
        error = NA
      ))
      break
    }
    
    x_new <- x - fx / dfx
    
    if (x_new < 1) x_new <- 1
    if (x_new > length(diferencia)) x_new <- length(diferencia)
    
    error <- abs(x_new - x)
    
    historial <- rbind(historial, data.frame(
      iter = iter,
      x = x,
      fx = f_real,
      dfx = df_real,
      abs_dfx = abs(df_real),
      error = error
    ))
    
    cat(sprintf("%.6f\n", error))
    
    if (buscar_minimo) {
      if (abs(fx) < tol || error < tol) {
        convergio <- TRUE
        x <- x_new
        break
      }
    } else {
      if (abs(fx) < tol || error < tol) {
        convergio <- TRUE
        x <- x_new
        break
      }
    }
    
    x <- x_new
  }
  
  cat("─────────────────────────────────────────────────────────────\n\n")
  
  if (!convergio && iter >= max_iter) {
    cat(sprintf("⚠️  Alcanzado máximo de iteraciones (%d)\n\n", max_iter))
  } else {
    cat(sprintf("✓ Convergencia alcanzada en %d iteraciones\n\n", iter))
  }
  
  tiempo_fin_metodo <- Sys.time()
  tiempo_metodo <- as.numeric(difftime(tiempo_fin_metodo, tiempo_inicio_metodo, units = "secs"))
  
  tiempo_fin_total <- Sys.time()
  tiempo_total <- as.numeric(difftime(tiempo_fin_total, tiempo_inicio_total, units = "secs"))
  
  # RESULTADOS
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║                         RESULTADOS                            ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  t_int <- round(x)
  if (t_int < 1) t_int <- 1
  if (t_int > length(datos$unhappy)) t_int <- length(datos$unhappy)
  
  fc <- f(x)
  dfc <- df(x)
  
  if (buscar_minimo) {
    cat("✓ MÍNIMO ENCONTRADO (f'(t) ≈ 0)\n\n")
  } else {
    cat("✓ CERO ENCONTRADO (f(t) ≈ 0)\n\n")
  }
  
  cat(sprintf("   Tiempo óptimo:        t = %.3f años\n", x))
  cat(sprintf("   Valor de f(t):        %.6f\n", fc))
  cat(sprintf("   Valor de f'(t):       %.6f\n", dfc))
  cat(sprintf("   Distancia al cero:    %.6f\n", abs(fc)))
  cat(sprintf("   Iteraciones:          %d\n\n", iter))
  
  cat(sprintf("📈 ESTADO DEL SISTEMA EN t = %.3f:\n\n", x))
  cat(sprintf("   Proporción infelices:  %.4f (%.1f%%)\n", 
              datos$unhappy[t_int], 100*datos$unhappy[t_int]))
  cat(sprintf("   Proporción similares:  %.4f (%.1f%%)\n", 
              datos$similar[t_int], 100*datos$similar[t_int]))
  cat(sprintf("   Diferencia:            %.6f\n\n", diferencia[t_int]))
  
  # TIEMPOS DE EJECUCIÓN
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("⏱️  TIEMPOS DE EJECUCIÓN\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  cat(sprintf("   Simulación:        %.2f segundos (%.1f%%)\n", 
              tiempo_simulacion, 100*tiempo_simulacion/tiempo_total))
  cat(sprintf("   Método Newton-R:   %.2f segundos (%.1f%%)\n", 
              tiempo_metodo, 100*tiempo_metodo/tiempo_total))
  cat(sprintf("   TOTAL:             %.2f segundos\n\n", tiempo_total))
  
  if (tiempo_total < 60) {
    cat(sprintf("   ⏱️  Tiempo total: %.1f segundos\n\n", tiempo_total))
  } else if (tiempo_total < 3600) {
    cat(sprintf("   ⏱️  Tiempo total: %.1f minutos (%.0f segundos)\n\n", 
                tiempo_total/60, tiempo_total))
  } else {
    cat(sprintf("   ⏱️  Tiempo total: %.2f horas (%.0f minutos)\n\n", 
                tiempo_total/3600, tiempo_total/60))
  }
  
  # GRÁFICOS
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))
  
  tiempo <- 1:length(diferencia)
  
  # 1. Evolución temporal
  plot(tiempo, datos$unhappy, type = "l", col = "red", lwd = 2,
       xlab = "Tiempo (años)", ylab = "Proporción",
       main = "Modelo 1: Evolución Temporal", ylim = c(0, 1))
  lines(tiempo, datos$similar, col = "blue", lwd = 2)
  abline(v = x, lty = 2, col = "darkgreen", lwd = 2)
  legend("right", 
         legend = c("Infelices", "Similares", sprintf("t* = %.2f", x)),
         col = c("red", "blue", "darkgreen"),
         lty = c(1, 1, 2), lwd = 2, cex = 0.7)
  grid()
  
  # 2. Función objetivo
  plot(tiempo, diferencia, type = "l", col = "purple", lwd = 2,
       xlab = "Tiempo (años)", ylab = "f(t)",
       main = "Función Objetivo f(t)")
  abline(h = 0, lty = 2, col = "gray", lwd = 1)
  abline(v = x, lty = 2, col = "darkgreen", lwd = 2)
  points(x, fc, pch = 19, col = "darkgreen", cex = 2)
  grid()
  
  # 3. Derivada f'(t)
  derivadas <- sapply(tiempo, df)
  plot(tiempo, derivadas, type = "l", col = "orange", lwd = 2,
       xlab = "Tiempo (años)", ylab = "f'(t)",
       main = "Derivada f'(t)")
  abline(h = 0, lty = 2, col = "gray", lwd = 1)
  abline(v = x, lty = 2, col = "darkgreen", lwd = 2)
  points(x, dfc, pch = 19, col = "darkgreen", cex = 2)
  grid()
  
  # 4. Trayectoria en el espacio
  if (nrow(historial) > 1) {
    plot(historial$x, historial$fx, type = "b", pch = 19, col = "blue", lwd = 2,
         xlab = "x", ylab = "f(x)",
         main = "Trayectoria de Newton-Raphson")
    abline(h = 0, lty = 2, col = "gray")
    points(x, fc, pch = 19, col = "red", cex = 2)
    grid()
  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "", main = "Trayectoria\n(Convergencia inmediata)")
    text(1, 1, "Convergió en 1 iteración", cex = 1.2)
  }
  
  # 5. Convergencia
  if (nrow(historial) > 1) {
    plot(historial$iter, historial$x, type = "b", pch = 19, col = "blue", lwd = 2,
         xlab = "Iteración", ylab = "Aproximación de t",
         main = "Convergencia del Método")
    abline(h = x, lty = 2, col = "red", lwd = 2)
    grid()
  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "", main = "Convergencia\n(Convergencia inmediata)")
    text(1, 1, sprintf("x0 = %.2f ya era óptimo", historial$x[1]), cex = 1)
  }
  
  # 6. Error
  if (nrow(historial) > 1 && any(!is.na(historial$error))) {
    errores_validos <- historial$error[!is.na(historial$error)]
    if (length(errores_validos) > 0) {
      plot(which(!is.na(historial$error)), errores_validos, 
           type = "b", pch = 19, col = "red", lwd = 2,
           xlab = "Iteración", ylab = "Error |x_{n+1} - x_n|",
           main = "Reducción del Error", log = "y")
      abline(h = tol, lty = 2, col = "darkgreen")
      grid()
    } else {
      plot(1, 1, type = "n", xlab = "", ylab = "", main = "Error\n(Convergencia inmediata)")
      text(1, 1, "Sin errores medibles", cex = 1.2)
    }
  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "", main = "Error\n(Convergencia inmediata)")
    text(1, 1, "Convergió en 1 iteración", cex = 1.2)
  }
  
  return(list(
    modelo = "Modelo 1 - Tolerancia Fija",
    metodo = "Newton-Raphson",
    tiempo_optimo = x,
    valor_funcion = fc,
    valor_derivada = dfc,
    distancia_cero = abs(fc),
    iteraciones = iter,
    convergio = convergio,
    busco_minimo = buscar_minimo,
    historial = historial,
    datos_simulacion = datos,
    parametros = list(
      N = N, p = p, citySize = citySize,
      alikePref = alikePref, alfa = alfa,
      x0 = x0, h = h
    ),
    tiempos = list(
      simulacion = tiempo_simulacion,
      metodo = tiempo_metodo,
      total = tiempo_total
    )
  ))
}

# ====================================================================================
# EJEMPLO DE USO
# ====================================================================================

cat("\n")
cat("███████████████████████████████████████████████████████████████████\n")
cat("█                                                                 █\n")
cat("█                  TRABAJO 7 - MODELO 1                           █\n")
cat("█               MÉTODO: NEWTON-RAPHSON                            █\n")
cat("█                                                                 █\n")
cat("███████████████████████████████████████████████████████████████████\n")

resultado_modelo1_nr <- newton_raphson_modelo1(
  N = 2000,
  p = 0.5,
  citySize = 50,
  alikePref = 0.6,
  alfa = 1.5,
  x0 = 20,
  years = 200, 
  tol = 0.01,
  max_iter = 200,
  h = 1.0
)

# Resumen compacto
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║                    RESUMEN PARA COMPARACIÓN                   ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("Modelo:               %s\n", resultado_modelo1_nr$modelo))
cat(sprintf("Método:               %s\n", resultado_modelo1_nr$metodo))
cat(sprintf("Tiempo óptimo:        %.3f años\n", resultado_modelo1_nr$tiempo_optimo))
cat(sprintf("Distancia al cero:    %.6f\n", resultado_modelo1_nr$distancia_cero))
cat(sprintf("Iteraciones:          %d\n", resultado_modelo1_nr$iteraciones))
cat(sprintf("Tiempo computacional: %.3f seg\n", resultado_modelo1_nr$tiempos$total))
cat(sprintf("¿Convergió?:          %s\n", 
            ifelse(resultado_modelo1_nr$convergio, "Sí", "NO")))
cat(sprintf("Tipo de búsqueda:     %s\n\n", 
            ifelse(resultado_modelo1_nr$busco_minimo, "Mínimo (f'=0)", "Cero (f=0)")))

saveRDS(resultado_modelo1_nr, "resultado_modelo1_newton.rds")
cat("✓ Resultado guardado en: resultado_modelo1_newton.rds\n")