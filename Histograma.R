# Función para extraer todas las raíces (tiempos t*) de una sola simulación.
obtener_raices_run <- function(resultado, years) {
  tiempo <- 1:years
  # f(t) = Infelices - Similares
  diferencia <- resultado$unhappy - resultado$similar
  
  raices_encontradas <- c()
  
  # 1. Encontrar el Mínimo Global Discreto (su t_referencia*)
  idx_cero <- which.min(abs(diferencia))
  t_minimo <- tiempo[idx_cero]
  raices_encontradas <- c(raices_encontradas, t_minimo)
  
  # 2. Encontrar y calcular los cruces interpolados (raíces continuas)
  cruces <- which(diff(sign(diferencia)) != 0)
  
  if (length(cruces) > 0) {
    for (idx in cruces) {
      t1 <- tiempo[idx]
      t2 <- tiempo[idx + 1]
      y1 <- diferencia[idx]
      y2 <- diferencia[idx + 1]
      
      # Interpolación lineal para encontrar t_cruce (donde f(t)=0)
      t_cruce <- t1 - y1 * (t2 - t1) / (y2 - y1)
      raices_encontradas <- c(raices_encontradas, t_cruce)
    }
  }
  
  # Devolver todas las raíces (únicas) encontradas
  return(unique(raices_encontradas))
}


histograma_raices <- function(n = 100, N=2000, p=.5, citySize=50, years=30, alikePref=.6, alfa=1.5) {
  
  cat(sprintf("\n=== INICIANDO SIMULACIÓN (%d simulaciones) ===\n", n))
  
  # Vector para almacenar todas las raices encontradas
  todas_las_raices <- c()
  
  # Bucle Monte Carlo
  for (i in 1:n) {
    # Ejecutar la simulación (usamos la función que ya tiene)
    resultado_i <- simular_modelo1(N, p, citySize, years, alikePref, alfa)
    
    # Extraer las raíces de esta ejecución
    raices_i <- obtener_raices_run(resultado_i, years)
    
    # Acumular las raíces en el vector principal
    todas_las_raices <- c(todas_las_raices, raices_i)
    
    if (i %% 5 == 0) {
      cat(sprintf("Procesado %d de %d simulaciones...\n", i, n))
    }
  }
  
  cat("\n✓ Extracción de raíces completada.\n")
  
  # ====================================================================================
  # GENERAR HISTOGRAMA
  # ====================================================================================
  
  df_raices <- data.frame(Tiempo = todas_las_raices)
  
  # AÑADIDO: Creación de la columna de tiempo entero
  df_raices$Tiempo_entero <- floor(df_raices$Tiempo)
  
  media_t <- mean(todas_las_raices)
  desvio_estandar_t <- sd(todas_las_raices)
  
  # CORREGIDO: Se usa Tiempo_entero en el aes(x) para el histograma discreto
  grafico_hist <- ggplot(df_raices, aes(x = Tiempo_entero)) +
    geom_histogram(binwidth = 1, fill = "lightblue", color = "white", alpha = 0.8) +
    geom_vline(xintercept = media_t, color = "red",
               linetype = "dashed", linewidth = 1.2) +
    labs(
      title = "Distribución del Tiempo de Equilibrio",
      subtitle = sprintf("Basado en %d simulaciones | Media de t: %.2f años | Desv. Est.: %.2f", n, media_t, desvio_estandar_t),
      x = "Tiempo de Equilibrio (años)",
      y = "Frecuencia de Raíces"
    ) +
    scale_x_continuous(breaks = unique(df_raices$Tiempo_entero)) +
    theme_minimal(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  
  print(grafico_hist)
  
  cat(sprintf("\n=== RESULTADOS DE LA SIMULACION ===\n"))
  cat(sprintf("Total de raíces encontradas: %d\n", length(todas_las_raices)))
  cat(sprintf("Tiempo medio de equilibrio (Media t): %.3f años\n", media_t))
  cat(sprintf("Desviación estándar de t*: %.3f\n", desvio_estandar_t))
  cat(sprintf("Moda: %d\n", as.numeric(names(sort(table(df_raices$Tiempo_entero), decreasing = TRUE)[1]))))
  
}