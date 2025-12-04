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
  
  # MODELO 1: Tolerancia FIJA durante toda la simulación
  alikePref1 <- alikePref      # Grupo 1
  alikePref2 <- alikePref * alfa  # Grupo 2
  
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
    
    # Mover agentes insatisfechos
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
# MÉTODO DE BISECCIÓN PARA MODELO 1
# ====================================================================================

biseccion_modelo1 <- function(N, p, citySize, alikePref, alfa, 
                              a, b, tol = 0.01, max_iter) {
  
  cat("\n")
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║              MODELO 1: TOLERANCIA FIJA                        ║\n")
  cat("║              MÉTODO: BISECCIÓN                                ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  # Medir tiempo total
  tiempo_inicio_total <- Sys.time()
  
  # Medir tiempo de simulación
  cat("⏳ Ejecutando simulación del Modelo 1...\n")
  tiempo_inicio_sim <- Sys.time()
  datos <- simular_modelo1(N, p, citySize, b, alikePref, alfa)
  tiempo_fin_sim <- Sys.time()
  tiempo_simulacion <- as.numeric(difftime(tiempo_fin_sim, tiempo_inicio_sim, units = "secs"))
  diferencia <- datos$unhappy - datos$similar
  
  cat(sprintf("✓ Simulación completada en %.2f segundos\n\n", tiempo_simulacion))
  
  # Medir tiempo del método de bisección
  tiempo_inicio_metodo <- Sys.time()
  
  # Función objetivo
  f <- function(t) {
    t_int <- round(t)
    if (t_int < 1) t_int <- 1
    if (t_int > length(diferencia)) t_int <- length(diferencia)
    return(diferencia[t_int])
  }
  
  cat("📊 PARÁMETROS:\n")
  cat(sprintf("   • N = %d agentes (%.1f%% ocupación)\n", N, 100*N/(citySize^2)))
  cat(sprintf("   • Grupo 1: %.0f%% (tolerancia = %.2f)\n", 100*p, alikePref))
  cat(sprintf("   • Grupo 2: %.0f%% (tolerancia = %.2f)\n", 100*(1-p), alikePref*alfa))
  cat(sprintf("   • Ciudad: %dx%d parcelas\n", citySize, citySize))
  cat(sprintf("   • Intervalo: [%d, %d] años\n\n", a, b))
  
  # Verificar cambio de signo
  fa <- f(a)
  fb <- f(b)
  
  cat("🔍 ANÁLISIS INICIAL:\n")
  cat(sprintf("   f(%d) = %.6f\n", a, fa))
  cat(sprintf("   f(%d) = %.6f\n", b, fb))
  
  historial <- data.frame()
  iter <- 0
  hay_cero <- (fa * fb < 0)
  
  if (!hay_cero) {
    cat("\n⚠️  NO HAY CAMBIO DE SIGNO\n")
    cat("   → Buscando punto más cercano a cero\n\n")
    cat("─────────────────────────────────────────────────────────────\n")
    cat("ITER      t         f(t)        |f(t)|      INTERVALO\n")
    cat("─────────────────────────────────────────────────────────────\n")
    
    # Búsqueda ternaria
    while (iter < max_iter && (b - a) > tol) {
      iter <- iter + 1
      
      m1 <- a + (b - a) / 3
      m2 <- b - (b - a) / 3
      
      fm1 <- abs(f(m1))
      fm2 <- abs(f(m2))
      
      if (fm1 < fm2) {
        b <- m2
        c <- m1
        fc <- f(m1)
      } else {
        a <- m1
        c <- m2
        fc <- f(m2)
      }
      
      error <- b - a
      
      historial <- rbind(historial, data.frame(
        iter = iter, 
        t = c, 
        ft = fc, 
        abs_ft = abs(fc), 
        intervalo = error
      ))
      
      cat(sprintf("%4d   %7.2f   %10.6f   %10.6f   %10.2f\n", 
                  iter, c, fc, abs(fc), error))
    }
    
    if (iter >= max_iter) {
      cat("\n⚠️  Alcanzado máximo de iteraciones\n")
    }
    
  } else {
    cat("\n✓ HAY CAMBIO DE SIGNO\n")
    cat("   → Buscando cero exacto\n\n")
    cat("─────────────────────────────────────────────────────────────\n")
    cat("ITER     a        b         c        f(c)        ERROR\n")
    cat("─────────────────────────────────────────────────────────────\n")
    
    # Bisección 
    while (iter < max_iter) {
      iter <- iter + 1
      c <- (a + b) / 2
      fc <- f(c)
      error <- abs(b - a) / 2
      
      historial <- rbind(historial, data.frame(
        iter = iter,
        a = a,
        b = b,
        t = c,
        ft = fc,
        abs_ft = abs(fc),
        error = error
      ))
      
      cat(sprintf("%4d  %7.2f  %7.2f  %7.2f  %10.6f  %10.4f\n", 
                  iter, a, b, c, fc, error))
      
      if (abs(fc) < 0.001 || error < tol) break
      
      if (f(a) * fc < 0) {
        b <- c
      } else {
        a <- c
      }
    }
    
    if (iter >= max_iter) {
      cat("\n⚠️  Alcanzado máximo de iteraciones\n")
    }
  }
  
  cat("─────────────────────────────────────────────────────────────\n\n")
  
  tiempo_fin_metodo <- Sys.time()
  tiempo_metodo <- as.numeric(difftime(tiempo_fin_metodo, tiempo_inicio_metodo, units = "secs"))
  
  tiempo_fin_total <- Sys.time()
  tiempo_total <- as.numeric(difftime(tiempo_fin_total, tiempo_inicio_total, units = "secs"))
  
  # RESULTADOS
  cat("╔═══════════════════════════════════════════════════════════════╗\n")
  cat("║                         RESULTADOS                            ║\n")
  cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
  
  t_int <- round(c)
  if (t_int < 1) t_int <- 1
  if (t_int > length(datos$unhappy)) t_int <- length(datos$unhappy)
  
  if (hay_cero) {
    cat("✓ CERO ENCONTRADO\n\n")
  } else {
    cat("✓ PUNTO MÁS CERCANO A CERO\n\n")
  }
  
  cat(sprintf("   Tiempo óptimo:        t = %.2f años\n", c))
  cat(sprintf("   Valor de f(t):        %.6f\n", fc))
  cat(sprintf("   Distancia al cero:    %.6f\n", abs(fc)))
  cat(sprintf("   Iteraciones:          %d\n\n", iter))
  
  cat(sprintf("📈 ESTADO DEL SISTEMA EN t = %.2f:\n", c))
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
  cat(sprintf("   Método bisección:  %.2f segundos (%.1f%%)\n", 
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
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  
  tiempo <- 1:length(diferencia)
  
  # 1. Evolución temporal
  plot(tiempo, datos$unhappy, type = "l", col = "red", lwd = 2,
       xlab = "Tiempo (años)", ylab = "Proporción",
       main = "Modelo 1: Evolución Temporal", ylim = c(0, 1))
  lines(tiempo, datos$similar, col = "blue", lwd = 2)
  abline(v = c, lty = 2, col = "darkgreen", lwd = 2)
  legend("right", 
         legend = c("Infelices", "Similares", sprintf("t* = %.1f", c)),
         col = c("red", "blue", "darkgreen"),
         lty = c(1, 1, 2), lwd = 2, cex = 0.8)
  grid()
  
  # 2. Función objetivo
  plot(tiempo, diferencia, type = "l", col = "purple", lwd = 2,
       xlab = "Tiempo (años)", ylab = "f(t) = Infelices - Similares",
       main = "Función Objetivo")
  abline(h = 0, lty = 2, col = "gray", lwd = 1)
  abline(v = c, lty = 2, col = "darkgreen", lwd = 2)
  points(c, fc, pch = 19, col = "darkgreen", cex = 2)
  text(c, fc, sprintf("  (%.1f, %.4f)", c, fc), pos = 4, cex = 0.8)
  grid()
  
  # 3. Convergencia
  plot(historial$iter, historial$t, type = "b", pch = 19, col = "blue", lwd = 2,
       xlab = "Iteración", ylab = "Aproximación de t",
       main = "Convergencia del Método")
  abline(h = c, lty = 2, col = "red", lwd = 2)
  grid()
  
  # 4. Distancia al cero
  plot(historial$iter, historial$abs_ft, type = "b", pch = 19, col = "red", lwd = 2,
       xlab = "Iteración", ylab = "|f(t)|",
       main = "Distancia al Cero", log = "y")
  grid()
  
  # 5. Visualización espacial final
  par(mfrow = c(1, 1), mar = c(2, 2, 3, 2))
  image(datos$city_final, col = c("white", "red", "blue"), axes = FALSE,
        main = sprintf("Modelo 1: Distribución Espacial Final (t=%f)", b))
  legend("topright", legend = c("Vacío", "Grupo 1", "Grupo 2"),
         fill = c("white", "red", "blue"), cex = 0.8)
  
  return(list(
    modelo = "Modelo 1 - Tolerancia Fija",
    metodo = "Bisección",
    tiempo_optimo = c,
    valor_funcion = fc,
    distancia_cero = abs(fc),
    iteraciones = iter,
    hay_cero_real = hay_cero,
    historial = historial,
    datos_simulacion = datos,
    parametros = list(
      N = N, p = p, citySize = citySize,
      alikePref = alikePref, alfa = alfa
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
cat("█                    TRABAJO 7 - MODELO 1                         █\n")
cat("█                  MÉTODO: BISECCIÓN                              █\n")
cat("█                                                                 █\n")
cat("███████████████████████████████████████████████████████████████████\n")

# Ejecutar bisección para Modelo 1
resultado_modelo1 <- biseccion_modelo1(
  N = 2000,
  p = 0.5,
  citySize = 50,
  alikePref = 0.6,
  alfa = 1.5,
  a = 1,
  b = 200,
  tol = 0.005,
  max_iter = 100
)

# Resumen compacto
cat("\n")
cat("╔═══════════════════════════════════════════════════════════════╗\n")
cat("║                    RESUMEN PARA COMPARACIÓN                   ║\n")
cat("╚═══════════════════════════════════════════════════════════════╝\n\n")
cat(sprintf("Modelo:               %s\n", resultado_modelo1$modelo))
cat(sprintf("Método:               %s\n", resultado_modelo1$metodo))
cat(sprintf("Tiempo óptimo:        %.2f años\n", resultado_modelo1$tiempo_optimo))
cat(sprintf("Distancia al cero:    %.6f\n", resultado_modelo1$distancia_cero))
cat(sprintf("Iteraciones:          %d\n", resultado_modelo1$iteraciones))
cat(sprintf("Tiempo computacional: %.3f seg\n", resultado_modelo1$tiempos$total))
cat(sprintf("¿Cero real?:          %s\n\n", 
            ifelse(resultado_modelo1$hay_cero_real, "SÍ", "NO")))

# Guardar resultado
saveRDS(resultado_modelo1, "resultado_modelo1_biseccion.rds")
cat("✓ Resultado guardado en: resultado_modelo1_biseccion.rds\n")