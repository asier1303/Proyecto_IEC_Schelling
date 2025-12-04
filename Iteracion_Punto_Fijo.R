#######################################################################################
#  TRABAJO 7 - MODELO 1: TOLERANCIA FIJA                                            #
#  Método: PUNTO FIJO (Versión Clásica Corregida)                                   #
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
# FUNCIÓN DE PUNTO FIJO CLÁSICA (Versión General)
# ====================================================================================

pf <- function(FUN, x0, tol = 1e-5, maxiter = 1000) {
  x_new <- FUN(x0)
  iter <- 1
  
  historial <- data.frame(iter = 0, x = x0, x_new = x_new, error = abs(x_new - x0))
  
  while (abs(x_new - x0) > tol & iter < maxiter) {
    x0 <- x_new
    x_new <- FUN(x0)
    iter <- iter + 1
    
    historial <- rbind(historial, data.frame(
      iter = iter - 1,
      x = x0,
      x_new = x_new,
      error = abs(x_new - x0)
    ))
    
    cat(paste('x_', iter, sep = ''), '=', x_new, '\n')
  }
  
  if (iter == maxiter) {
    cat('No hubo convergencia luego de ', iter, ' iteraciones ', '\n')
  } else {
    cat('Si hubo convergencia luego de ', iter, ' iteraciones ', '\n')
  }
  
  salida <- list(pf = x_new, iter = iter, historial = historial)
  return(salida)
}

# ====================================================================================
# MÉTODO DE PUNTO FIJO PARA MODELO 1 - VERSIÓN CORREGIDA
# ====================================================================================

punto_fijo_modelo1 <- function(N, p, citySize, alikePref, alfa, 
                               x0, years, tol = 0.01, maxiter = 50, lambda = 0.1) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║              MODELO 1: TOLERANCIA FIJA                        ║\n")
  cat("║              MÉTODO: PUNTO FIJO (VERSIÓN CLÁSICA)             ║\n")
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
  
  # Función objetivo f(t) con interpolación lineal para suavizar
  f <- function(t) {
    if (t <= 1) return(diferencia[1])
    if (t >= length(diferencia)) return(diferencia[length(diferencia)])
    
    t_floor <- floor(t)
    t_ceil <- ceiling(t)
    
    if (t_floor == t_ceil) {
      return(diferencia[t_floor])
    } else {
      peso <- t - t_floor
      return(diferencia[t_floor] * (1 - peso) + diferencia[t_ceil] * peso)
    }
  }
  
  cat("📊 PARÁMETROS:\n")
  cat(sprintf("   • N = %d agentes (%.1f%% ocupación)\n", N, 100*N/(citySize^2)))
  cat(sprintf("   • Grupo 1: %.0f%% (tolerancia = %.2f)\n", 100*p, alikePref))
  cat(sprintf("   • Grupo 2: %.0f%% (tolerancia = %.2f)\n", 100*(1-p), alikePref*alfa))
  cat(sprintf("   • Ciudad: %dx%d parcelas\n", citySize, citySize))
  cat(sprintf("   • Valor inicial: x0 = %.2f años\n", x0))
  cat(sprintf("   • Parámetro λ: %.3f\n", lambda))
  cat(sprintf("   • Tolerancia: %.4f\n", tol))
  cat(sprintf("   • Máximo iteraciones: %d\n\n", maxiter))
  
  cat("🔍 MÉTODO:\n")
  cat("   Iteración de Punto Fijo Clásica: x_{n+1} = g(x_n)\n")
  cat("   Función g(x) transformada para garantizar convergencia\n")
  cat("   Criterio: |x_{n+1} - x_n| < tol\n\n")
  
  # Verificar si hay cero real
  valores <- sapply(1:length(diferencia), f)
  tiene_cero <- any(diff(sign(valores)) != 0)
  
  if (!tiene_cero) {
    cat("⚠️  NO HAY CAMBIO DE SIGNO en la función\n")
    cat("   → Buscando MÍNIMO de |f(t)|\n\n")
  } else {
    cat("✓ HAY CAMBIO DE SIGNO\n")
    cat("   → Buscando CERO de f(t) = 0\n\n")
  }
  
  # Encontrar intervalo aproximado del cero
  if (tiene_cero) {
    for (i in 1:(length(diferencia)-1)) {
      if (diferencia[i] * diferencia[i+1] < 0) {
        x0 <- (i + i + 1) / 2
        cat(sprintf("   Ajustando x0 al intervalo con cambio de signo: x0 = %.2f\n\n", x0))
        break
      }
    }
  }
  
  # Definir función g(x)
  g <- function(x) {
    fx <- f(x)
    
    h <- 0.5
    x_plus <- min(x + h, length(diferencia))
    x_minus <- max(x - h, 1)
    
    dfx <- (f(x_plus) - f(x_minus)) / (x_plus - x_minus)
    
    if (abs(dfx) < 1e-6) {
      if (fx > 0) {
        x_new <- x - lambda
      } else {
        x_new <- x + lambda
      }
    } else {
      x_new <- x - lambda * fx / dfx
    }
    
    if (x_new < 1) x_new <- 1
    if (x_new > length(diferencia)) x_new <- length(diferencia)
    
    return(x_new)
  }
  
  cat("─────────────────────────────────────────────────────────────\n")
  cat("ITERACIONES DEL MÉTODO:\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat(paste('x_0 = ', x0, '\n', sep = ''))
  
  resultado_pf <- pf(FUN = g, x0 = x0, tol = tol, maxiter = maxiter)
  
  cat("─────────────────────────────────────────────────────────────\n\n")
  
  tiempo_fin_metodo <- Sys.time()
  tiempo_metodo <- as.numeric(difftime(tiempo_fin_metodo, tiempo_inicio_metodo, units = "secs"))
  
  tiempo_fin_total <- Sys.time()
  tiempo_total <- as.numeric(difftime(tiempo_fin_total, tiempo_inicio_total, units = "secs"))
  
  # Extraer solución
  x_final <- resultado_pf$pf
  iter_final <- resultado_pf$iter
  convergio <- (iter_final < maxiter)
  
  # RESULTADOS
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║                         RESULTADOS                            ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  t_int <- round(x_final)
  if (t_int < 1) t_int <- 1
  if (t_int > length(datos$unhappy)) t_int <- length(datos$unhappy)
  
  fc <- f(x_final)
  
  if (!tiene_cero) {
    cat("✓ MÍNIMO ENCONTRADO\n\n")
  } else {
    cat("✓ PUNTO FIJO ENCONTRADO\n\n")
  }
  
  cat(sprintf("   Tiempo óptimo:        t = %.3f años\n", x_final))
  cat(sprintf("   Valor de f(t):        %.6f\n", fc))
  cat(sprintf("   Distancia al cero:    %.6f\n", abs(fc)))
  cat(sprintf("   Iteraciones:          %d\n\n", iter_final))
  
  cat(sprintf("📈 ESTADO DEL SISTEMA EN t = %.3f:\n\n", x_final))
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
  cat(sprintf("   Método punto fijo: %.2f segundos (%.1f%%)\n", 
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
  
  # Crear historial extendido con valores de f(x)
  historial <- resultado_pf$historial
  historial$fx <- sapply(historial$x_new, f)
  historial$abs_fx <- abs(historial$fx)
  
  # GRÁFICOS
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))
  
  tiempo <- 1:length(diferencia)
  
  # 1. Evolución temporal
  plot(tiempo, datos$unhappy, type = "l", col = "red", lwd = 2,
       xlab = "Tiempo (años)", ylab = "Proporción",
       main = "Modelo 1: Evolución Temporal", ylim = c(0, 1))
  lines(tiempo, datos$similar, col = "blue", lwd = 2)
  abline(v = x_final, lty = 2, col = "darkgreen", lwd = 2)
  legend("right", 
         legend = c("Infelices", "Similares", sprintf("t* = %.2f", x_final)),
         col = c("red", "blue", "darkgreen"),
         lty = c(1, 1, 2), lwd = 2, cex = 0.7)
  grid()
  
  # 2. Función objetivo
  plot(tiempo, diferencia, type = "l", col = "purple", lwd = 2,
       xlab = "Tiempo (años)", ylab = "f(t)",
       main = "Función Objetivo f(t)")
  abline(h = 0, lty = 2, col = "gray", lwd = 1)
  abline(v = x_final, lty = 2, col = "darkgreen", lwd = 2)
  points(x_final, fc, pch = 19, col = "darkgreen", cex = 2)
  grid()
  
  # 3. Trayectoria de iteraciones
  if (nrow(historial) > 1) {
    plot(historial$x_new, historial$fx, type = "b", pch = 19, col = "blue", lwd = 2,
         xlab = "x_n", ylab = "f(x_n)",
         main = "Trayectoria de Punto Fijo")
    abline(h = 0, lty = 2, col = "gray")
    points(x_final, fc, pch = 19, col = "red", cex = 2)
    grid()
  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "", 
         main = "Trayectoria\n(Convergencia inmediata)")
    text(1, 1, "Convergió en 1 iteración", cex = 1.2)
  }
  
  # 4. Convergencia de x_n
  if (nrow(historial) > 1) {
    plot(historial$iter, historial$x_new, type = "b", pch = 19, col = "blue", lwd = 2,
         xlab = "Iteración", ylab = "x_n",
         main = "Convergencia del Método")
    abline(h = x_final, lty = 2, col = "red", lwd = 2)
    grid()
  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "", 
         main = "Convergencia\n(Convergencia inmediata)")
    text(1, 1, sprintf("x0 = %.2f ya era óptimo", historial$x[1]), cex = 1)
  }
  
  # 5. Error absoluto
  if (nrow(historial) > 1) {
    errores_validos <- historial$error[historial$error > 0]
    if (length(errores_validos) > 0) {
      iters_validos <- historial$iter[historial$error > 0]
      plot(iters_validos, errores_validos, 
           type = "b", pch = 19, col = "red", lwd = 2,
           xlab = "Iteración", ylab = "Error |x_{n+1} - x_n|",
           main = "Reducción del Error", log = "y")
      abline(h = tol, lty = 2, col = "darkgreen")
      grid()
    } else {
      plot(1, 1, type = "n", xlab = "", ylab = "", 
           main = "Error\n(Convergencia inmediata)")
      text(1, 1, "Sin errores medibles", cex = 1.2)
    }
  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "", 
         main = "Error\n(Convergencia inmediata)")
    text(1, 1, "Convergió en 1 iteración", cex = 1.2)
  }
  
  # 6. Distancia al cero
  if (nrow(historial) > 1) {
    plot(historial$iter, historial$abs_fx, type = "b", pch = 19, col = "purple", lwd = 2,
         xlab = "Iteración", ylab = "|f(x_n)|",
         main = "Distancia al Cero", log = "y")
    abline(h = tol, lty = 2, col = "darkgreen")
    grid()
  } else {
    plot(1, 1, type = "n", xlab = "", ylab = "", 
         main = "Distancia al Cero\n(Convergencia inmediata)")
    text(1, 1, "Convergió en 1 iteración", cex = 1.2)
  }
  
  return(list(
    modelo = "Modelo 1 - Tolerancia Fija",
    metodo = "Punto Fijo (Clásico)",
    tiempo_optimo = x_final,
    valor_funcion = fc,
    distancia_cero = abs(fc),
    iteraciones = iter_final,
    convergio = convergio,
    tiene_cero = tiene_cero,
    historial = historial,
    datos_simulacion = datos,
    parametros = list(
      N = N, p = p, citySize = citySize,
      alikePref = alikePref, alfa = alfa,
      x0 = x0, lambda = lambda
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
cat("█                                                                   █\n")
cat("█                  TRABAJO 7 - MODELO 1                            █\n")
cat("█               MÉTODO: PUNTO FIJO (CLÁSICO)                       █\n")
cat("█                                                                   █\n")
cat("███████████████████████████████████████████████████████████████████\n")

resultado_modelo1_pf <- punto_fijo_modelo1(
  N = 2000,
  p = 0.5,
  citySize = 50,
  alikePref = 0.6,
  alfa = 1.5,
  x0 = 20,
  years = 250,
  tol = 0.01,
  maxiter = 2,
  lambda = 0.1
)

# Resumen compacto
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║                    RESUMEN PARA COMPARACIÓN                    ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("Modelo:               %s\n", resultado_modelo1_pf$modelo))
cat(sprintf("Método:               %s\n", resultado_modelo1_pf$metodo))
cat(sprintf("Tiempo óptimo:        %.3f años\n", resultado_modelo1_pf$tiempo_optimo))
cat(sprintf("Distancia al cero:    %.6f\n", resultado_modelo1_pf$distancia_cero))
cat(sprintf("Iteraciones:          %d\n", resultado_modelo1_pf$iteraciones))
cat(sprintf("Tiempo computacional: %.3f seg\n", resultado_modelo1_pf$tiempos$total))
cat(sprintf("¿Convergió?:          %s\n", 
            ifelse(resultado_modelo1_pf$convergio, "Sí", "NO")))
cat(sprintf("¿Tiene cero real?:    %s\n\n", 
            ifelse(resultado_modelo1_pf$tiene_cero, "Sí", "NO")))

saveRDS(resultado_modelo1_pf, "resultado_modelo1_puntofijo.rds")
cat("✓ Resultado guardado en: resultado_modelo1_puntofijo.rds\n")