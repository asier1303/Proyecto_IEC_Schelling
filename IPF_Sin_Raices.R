#######################################################################################
#  TRABAJO 7 - MODELO 1: TOLERANCIA FIJA                                             #
#  Método: PUNTO FIJO PARA CASOS SIN CEROS                                           #
#  Busca el mínimo de |f(t)| cuando las curvas no se cruzan                          #
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
# FUNCIÓN DE PUNTO FIJO PARA MINIMIZACIÓN
# ====================================================================================

pf_minimo <- function(FUN, x0, tol = 1e-5, maxiter = 1000) {
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
    
    cat(sprintf("x_%d = %.3f\n", iter, x_new))
  }
  
  if (iter == maxiter) {
    cat('⚠️  No hubo convergencia luego de ', iter, ' iteraciones\n')
  } else {
    cat('✓ Convergencia alcanzada luego de ', iter, ' iteraciones\n')
  }
  
  salida <- list(pf = x_new, iter = iter, historial = historial)
  return(salida)
}

# ====================================================================================
# MÉTODO DE PUNTO FIJO ESPECIALIZADO PARA CASOS SIN CEROS
# ====================================================================================

punto_fijo_sin_ceros <- function(N, p, citySize, alikePref, alfa, 
                                 x0, years, tol = 0.01, maxiter = 50, lambda = 0.1) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║              MODELO 1: TOLERANCIA FIJA                        ║\n")
  cat("║     MÉTODO: PUNTO FIJO PARA MINIMIZACIÓN (SIN CEROS)          ║\n")
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
  
  # Función objetivo f(t) con interpolación lineal
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
  
  # VERIFICAR QUE NO HAY CEROS
  valores <- sapply(1:length(diferencia), f)
  tiene_cero <- any(diff(sign(valores)) != 0)
  
  if (tiene_cero) {
    cat("❌ ADVERTENCIA: Se detectaron cambios de signo\n")
    cat("   Este método está diseñado para casos SIN ceros\n")
    cat("   Considera usar otro método\n\n")
  } else {
    cat("✓ Confirmado: No hay cambios de signo\n")
    cat("   → Buscando MÍNIMO de |f(t)|\n\n")
  }
  
  cat("🔍 ESTRATEGIA DE MINIMIZACIÓN:\n")
  cat("   Función g(x) = x - λ·sign(f)·f'(x)\n")
  cat("   Minimiza |f(t)| usando descenso por gradiente\n")
  cat("   Criterio: |x_{n+1} - x_n| < tol\n\n")
  
  # BUSCAR PUNTO INICIAL ÓPTIMO
  cat("🔍 Buscando punto inicial óptimo...\n")
  abs_vals <- abs(diferencia)
  x0_optimo <- which.min(abs_vals)
  f_min_global <- diferencia[x0_optimo]
  
  cat(sprintf("   Mínimo global aproximado: t = %d, |f(t)| = %.6f\n", 
              x0_optimo, abs(f_min_global)))
  cat(sprintf("   Usando x0 = %.2f\n\n", x0_optimo))
  
  x0 <- x0_optimo
  
  # FUNCIÓN g(x) ESPECIALIZADA PARA MINIMIZACIÓN
  g <- function(x) {
    fx <- f(x)
    
    h <- 0.5
    x_plus <- min(x + h, length(diferencia))
    x_minus <- max(x - h, 1)
    dfx <- (f(x_plus) - f(x_minus)) / (x_plus - x_minus)
    
    if (abs(dfx) < 1e-8) {
      x_new <- x
    } else {
      paso <- lambda * sign(fx) * dfx
      
      paso_max <- 2.0
      if (abs(paso) > paso_max) {
        paso <- sign(paso) * paso_max
      }
      
      x_new <- x - paso
    }
    
    if (x_new < 1) x_new <- 1
    if (x_new > length(diferencia)) x_new <- length(diferencia)
    
    return(x_new)
  }
  
  cat("─────────────────────────────────────────────────────────────\n")
  cat("ITERACIONES DE MINIMIZACIÓN:\n")
  cat("─────────────────────────────────────────────────────────────\n")
  cat(sprintf("x_0 = %.3f, f(x_0) = %.6f\n", x0, f(x0)))
  
  resultado_pf <- pf_minimo(FUN = g, x0 = x0, tol = tol, maxiter = maxiter)
  
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
  
  cat("✓ MÍNIMO ENCONTRADO (minimizando |f(t)|)\n\n")
  
  cat(sprintf("   Tiempo óptimo:        t = %.3f años\n", x_final))
  cat(sprintf("   Valor de f(t):        %.6f\n", fc))
  cat(sprintf("   Distancia al cero:    %.6f\n", abs(fc)))
  cat(sprintf("   Iteraciones:          %d\n\n", iter_final))
  
  # COMPARACIÓN CON MÍNIMO GLOBAL
  idx_min_global <- which.min(abs(diferencia))
  f_min_global <- diferencia[idx_min_global]
  
  cat("📈 COMPARACIÓN CON BÚSQUEDA EXHAUSTIVA:\n")
  cat(sprintf("   • Punto fijo:    t = %.3f, |f(t)| = %.6f\n", x_final, abs(fc)))
  cat(sprintf("   • Mínimo global: t = %d, |f(t)| = %.6f\n", idx_min_global, abs(f_min_global)))
  cat(sprintf("   • Diferencia:    Δt = %.1f años, Δ|f| = %.6f\n\n", 
              abs(x_final - idx_min_global), abs(abs(fc) - abs(f_min_global))))
  
  cat(sprintf("📊 ESTADO DEL SISTEMA EN t = %.3f:\n\n", x_final))
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
  
  # Crear historial extendido
  historial <- resultado_pf$historial
  historial$fx <- sapply(historial$x_new, f)
  historial$abs_fx <- abs(historial$fx)
  
  # GRÁFICOS ESPECIALIZADOS
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  tiempo <- 1:length(diferencia)
  
  # 1. Función objetivo con mínimo
  plot(tiempo, diferencia, type = "l", col = "purple", lwd = 2,
       xlab = "Tiempo (años)", ylab = "f(t)",
       main = "Búsqueda del Mínimo |f(t)|")
  abline(h = 0, lty = 2, col = "gray", lwd = 1)
  abline(v = x_final, lty = 2, col = "red", lwd = 2)
  abline(v = idx_min_global, lty = 2, col = "blue", lwd = 2)
  points(x_final, fc, pch = 19, col = "red", cex = 2)
  points(idx_min_global, f_min_global, pch = 19, col = "blue", cex = 2)
  legend("topright", 
         legend = c("f(t)", "Punto Fijo", "Mínimo Global"),
         col = c("purple", "red", "blue"),
         lty = c(1, 2, 2), lwd = 2)
  grid()
  
  # 2. Trayectoria de minimización
  if (nrow(historial) > 1) {
    plot(historial$x_new, historial$abs_fx, type = "b", pch = 19, col = "blue", lwd = 2,
         xlab = "x_n", ylab = "|f(x_n)|",
         main = "Trayectoria de Minimización")
    points(x_final, abs(fc), pch = 19, col = "red", cex = 2)
    points(idx_min_global, abs(f_min_global), pch = 19, col = "green", cex = 2)
    grid()
  }
  
  # 3. Convergencia
  if (nrow(historial) > 1) {
    plot(historial$iter, historial$abs_fx, type = "b", pch = 19, col = "red", lwd = 2,
         xlab = "Iteración", ylab = "|f(x_n)|",
         main = "Convergencia del Mínimo", log = "y")
    abline(h = abs(f_min_global), lty = 2, col = "blue")
    grid()
  }
  
  # 4. Error de posición
  if (nrow(historial) > 1) {
    plot(historial$iter, historial$error, type = "b", pch = 19, col = "green", lwd = 2,
         xlab = "Iteración", ylab = "|x_{n+1} - x_n|",
         main = "Error de Posición", log = "y")
    abline(h = tol, lty = 2, col = "red")
    grid()
  }
  
  return(list(
    modelo = "Modelo 1 - Tolerancia Fija",
    metodo = "Punto Fijo para Minimización",
    tiempo_optimo = x_final,
    valor_funcion = fc,
    distancia_cero = abs(fc),
    iteraciones = iter_final,
    convergio = convergio,
    minimo_global = list(
      tiempo = idx_min_global,
      valor = f_min_global,
      distancia = abs(f_min_global)
    ),
    diferencia_con_global = abs(x_final - idx_min_global),
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
# EJEMPLO DE USO - CASO SIN CEROS
# ====================================================================================

cat("\n")
cat("███████████████████████████████████████████████████████████████████\n")
cat("█                                                                 █\n")
cat("█        PUNTO FIJO PARA MINIMIZACIÓN (CASOS SIN CEROS)           █\n")
cat("█                                                                 █\n")
cat("███████████████████████████████████████████████████████████████████\n")

resultado_sin_ceros <- punto_fijo_sin_ceros(
  N = 2000,
  p = 0.5,
  citySize = 50,
  alikePref = 0.6,
  alfa = 1.5,
  x0 = 20,
  years =100,
  tol = 0.01,
  maxiter = 5,
  lambda = 0.05
)

# Resumen final
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║                      RESUMEN FINAL                            ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")

cat(sprintf("Método:               %s\n", resultado_sin_ceros$metodo))
cat(sprintf("Tiempo óptimo:        %.3f años\n", resultado_sin_ceros$tiempo_optimo))
cat(sprintf("Valor de f(t):        %.6f\n", resultado_sin_ceros$valor_funcion))
cat(sprintf("Distancia al cero:    %.6f\n", resultado_sin_ceros$distancia_cero))
cat(sprintf("Iteraciones:          %d\n", resultado_sin_ceros$iteraciones))
cat(sprintf("¿Convergió?:          %s\n", ifelse(resultado_sin_ceros$convergio, "SÍ", "NO")))
cat(sprintf("Diferencia con global: %.1f años\n", resultado_sin_ceros$diferencia_con_global))
cat(sprintf("Tiempo computacional: %.3f seg\n", resultado_sin_ceros$tiempos$total))

saveRDS(resultado_sin_ceros, "resultado_puntofijo_sin_ceros.rds")
cat("\n✓ Resultado guardado en: resultado_puntofijo_sin_ceros.rds\n")